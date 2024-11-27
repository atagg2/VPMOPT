import FLOWUnsteady as uns
import PreallocationTools as PT
vpm = uns.vpm
vlm = uns.vlm

struct VPMModel{TF, TI, TB}
    sim::uns.Simulation
    pfield::vpm.ParticleField
    cache
    dt::TF
    mass::TF
    I::TF
    g::TF
    rho::TF
    cog::Array{TF, 1}
    r_rotors::Array{TF, 2}
    dir_rotors::Array{TF, 2}
    A::Array{TF, 3}
    B::Array{TF, 3}
    x_dots::Array{TF, 2}
    i::TI
    t_sparse::Array{TF, 1}

    vlm_rlx::TF
    sigma_vlm_surf::TF
    sigma_rotor_surf::TF
    hubtiploss_correction
    vlm_init::TB
    sigmafactor_vpmonvlm::TF
    sigmafactor_vpm::TF
    sigma_vpm_overwrite::TF
    omit_shedding

end

function (m::VPMModel)(s)
    
    sim_dual, pfield_dual, static_pfield_dual, static_particles_function = m.cache[s]

    if m.i < 2
        vlm.setVinf(sim_dual.vehicle.system, (x,t)->zeros(eltype(s), 3))
    end

    sim_dual.vehicle = update_vehicle!(sim_dual.vehicle, m.sim.vehicle, s)

    update_pfield!(pfield_dual, m.pfield)

    dt = m.dt
    rho = m.rho
    mass = m.mass
    I = m.I
    g = m.g
    cog = m.cog
    r_rotors = m.r_rotors
    dir_rotors = m.dir_rotors

    Vinf = (x,t)->zeros(eltype(s), 3)
    wake_coupled = true
    sound_spd = 343

    x = s[1:6]
    u = s[7:end]

    vx, vz, theta_dot, rx, rz, theta = x

    V = [vx, 0.0, vz]
    W = [0.0, theta_dot, 0.0]

    sim_dual.nt = m.sim.nt
    sim_dual.t = m.sim.t

    RPM = (t->u[1]*30/pi,
            t->u[2]*30/pi,
            t->u[3]*30/pi)

    sim_dual.maneuver = uns.KinematicManeuver(sim_dual.maneuver.angle, RPM, sim_dual.maneuver.Vvehicle, sim_dual.maneuver.anglevehicle)

    org_np = vpm.get_np(pfield_dual)

    relax = m.pfield.relaxation != vpm.relaxation_none &&
    m.pfield.relaxation.nsteps_relax >= 1 &&
    m.sim.nt>0 && (m.sim.nt%m.pfield.relaxation.nsteps_relax == 0)

    # # Add static particles
    # remove = static_particles_function(pfield_dual, pfield.t, dt)

    # # Step in time solving governing equations
    # vpm.nextstep(pfield_dual, dt; relax=relax, custom_UJ=nothing)

    # # Remove static particles (assumes particles remained sorted)
    # if remove==nothing || remove
    #     for pi in vpm.get_np(pfield_dual):-1:(org_np+1)
    #         vpm.remove_particle(pfield_dual, pi)
    #     end
    # end

    prev_vlm_system = deepcopy(sim_dual.vehicle.vlm_system)

    uns.nextstep_kinematic(sim_dual.vehicle, m.dt)

    # Solver-specific pre-calculations
    uns.precalculations(sim_dual, Vinf, pfield_dual, m.sim.t, m.dt)

    uns.solve(sim_dual, Vinf, pfield_dual, wake_coupled, m.dt, m.vlm_rlx,
                m.sigma_vlm_surf, m.sigma_rotor_surf, m.rho, sound_spd,
                static_pfield_dual, m.hubtiploss_correction;
                init_sol=m.vlm_init, sigmafactor_vpmonvlm=m.sigmafactor_vpmonvlm,
                debug=false)

    # calculate wing forces 
    F, M = calc_aerodynamicforce_kuttajoukowski(sim_dual.vehicle.vlm_system, prev_vlm_system,
                                                pfield_dual, (x,t)->zeros(eltype(s),3), dt, rho, cog)

    # calculate rotor forces
    indices = [1,3,5,7,2,4,6,8,9]
    i = 1
    for rotor_system in sim_dual.vehicle.rotor_systems
        for rotor in rotor_system
            T, Q = uns.vlm.calc_thrust_torque(rotor)
            T_vec = T*dir_rotors[:,indices[i]]
            F += sim_dual.vehicle.system.Oaxis' * T_vec
            M += LA.cross(r_rotors[:,indices[i]] - m.cog, T_vec)   
            i += 1
        end
    end

    # add gravitational force
    F += [0.0, 0.0, -mass*g]

    # assemble state derivative
    x_dot = [F[1]/mass
             F[3]/mass
             M[2]/I
             vx
             vz
             theta_dot]

    # return state derivative
    return x_dot
end

function (m::VPMModel)(sim, pfield, T, DT; vprintln=nothing)

    if T > m.t_sparse[m.i]

        V = sim.vehicle.V
        W = sim.vehicle.W
        Oaxis = sim.vehicle.system.Oaxis

        theta = acosd(LA.dot(Oaxis[:,1], [1,0,0]))*-sign(Oaxis[3,1])

        x = [V[1], V[3], W[2], 0.0, 0.0, theta]
        RPM_front = sim.vehicle.rotor_systems[1][1].RPM
        RPM_back = sim.vehicle.rotor_systems[2][1].RPM
        RPM_pusher = sim.vehicle.rotor_systems[3][1].RPM
        u = [RPM_front, RPM_back, RPM_pusher]*pi/30

        @show x u

        update_pfield!(m.pfield, pfield)

        # A = FD.jacobian(m, [x; u], FD.JacobianConfig(m, [x; u], FD.Chunk(length([x; u]))), Val{false}())

        A = FD.jacobian(m, [x; u])

        m.A[:,:,m.i] = A[:,1:length(x)]
        m.B[:,:,m.i] = A[:,length(x)+1:end]
        m.x_dots[:,m.i] = m([x;u])
        m.i += 1

        f = open("A.txt", "w")
        print(f, "")
        for i in eachindex(A[:,1])
            for j in eachindex(A[1,:])
                print(f, A[i,j], " ")
            end
            print(f, "\n")
        end
        close(f)
        f = open("x_dot.txt", "w")
        print(f, "")
        for i in eachindex(m.x_dots[:,m.i])
            print(f, m.x_dots[i,m.i], " ")
        end
        print(f, "\n")
        close(f)
        

    end

    if T > DT*30
        stop
    end

    return false
end

function update_vehicle!(vehicle_dual::uns.UVLMVehicle{N, M, TD}, vehicle::uns.UVLMVehicle{N, M, TF}, s) where {N, M, TD, TF}

    x = s[1:6]
    u = s[7:end]
    vx, vz, theta_dot, rx, rz, theta = x

    # update vehicle velocity
    V = [vx, 0.0, vz]
    W = [0.0, theta_dot, 0.0]

    Oaxis = [cosd(theta) 0.0 sind(theta); 0.0 1.0 0.0; -sind(theta) 0.0 cosd(theta)]

    vlm.setcoordsystem(vehicle_dual.system, TD.(vehicle.system.O), Oaxis)



    update_wing_system!(vehicle_dual.system, vehicle.system, TD)

    update_wing_system!(vehicle_dual.wake_system, vehicle.wake_system, TD)
    
    update_wing_system!(vehicle_dual.vlm_system, vehicle.vlm_system, TD)

    for i in eachindex(vehicle.rotor_systems)
        for j in eachindex(vehicle.rotor_systems[i])
            update_wing!(vehicle_dual.rotor_systems[i][j], vehicle.rotor_systems[i][j], TD)
        end
    end

    vehicle_dual.V .= V
    vehicle_dual.W .= W

    # # update vehicle RPMs
    # for i in eachindex(u)
    #     for rotor in vehicle_dual.rotor_systems[i]
    #         rotor.RPM = u[i]*30/pi
    #     end
    # end

    return vehicle_dual
end

function update_wing_system!(system_dual::vlm.WingSystem, system::vlm.WingSystem, TD)

    system_dual.sol = system.sol

    for i in eachindex(system.wings)
        update_wing!(system_dual.wings[i], system.wings[i], TD)
    end
end

function update_wing!(wing_dual, wing, TD)

    wing_dual.sol = wing.sol


    if !(typeof(wing) <: vlm.Rotor)
        
        wing_dual._HSs = wing._HSs

        for i in eachindex(wing._HSs)
            for j in eachindex(wing._HSs[i])
                wing_dual._HSs[i][j] = TD.(wing._HSs[i][j])
            end
        end
    else
        wing_dual.RPM = TD.(wing.RPM)
        update_wing_system!(wing_dual._wingsystem, wing._wingsystem, TD)
        wing_dual.sol = wing.sol
    end
end

function update_pfield!(pfield_dual::uns.vpm.ParticleField{TD, FD, VD, TUinfD, SD, TkernelD, TUJD, TintegrationD, TRelaxationD, TGPUD}, 
                        pfield::uns.vpm.ParticleField{TF, F, V, TUinf, S, Tkernel, TUJ, Tintegration, TRelaxation, TGPU}) where {TD, FD, VD, TUinfD, SD, TkernelD, TUJD, TintegrationD, TRelaxationD, TGPUD, 
                                                                                                                                TF, F, V, TUinf, S, Tkernel, TUJ, Tintegration, TRelaxation, TGPU}
    pfield_dual.particles = TD.(pfield.particles)
    pfield_dual.t = pfield.t
    pfield_dual.nt = pfield.nt
    pfield_dual.np = pfield.np

end

function calc_aerodynamicforce_kuttajoukowski(vlm_system::vlm.AbstractWing{TF_design,TF_trajectory},
                                prev_vlm_system, pfield, Vinf, dt, rho, cog; t=0.0,
                                per_unit_span=false, spandir=[0, 1, 0],
                                include_trailingboundvortex=false,
                                Vout=nothing, lenout=nothing,
                                lencrit=-1, debug=false) where {TF_design,TF_trajectory}

    m = vlm.get_m(vlm_system)    # Number of horseshoes

    # Nodes of every horseshoe
    Ap = uns._get_Xs(vlm_system, "Ap")
    A = uns._get_Xs(vlm_system, "A")
    B = uns._get_Xs(vlm_system, "B")
    Bp = uns._get_Xs(vlm_system, "Bp")

    # Midpoints of bound vortices
    ApA = (Ap .+ A)/2
    AB = (A .+ B)/2
    BBp = (B .+ Bp)/2

    # Evaluate VPM on each midpoint
    Xs = vcat(ApA, AB, BBp)

    # Evaluate VLM on each midpoint
    Vvlm = vlm.Vind.(Ref(vlm_system), Xs; t=t, ign_col=true, ign_infvortex=true)

    # Evaluate Vinf on each midpoint
    Vinfs = Vinf.(Xs, t)

    # Evaluate kinematic velocity on each node
    Vtran = uns._Vkinematic(vlm_system, prev_vlm_system, dt; t=t,
                                                targetX=["Ap", "A", "B", "Bp"])

    # Pre-allocate memory
    TF = promote_type(TF_design, TF_trajectory)
    Ftot = [zeros(TF,3) for i in 1:m]
    if debug
        Vvpms, Vvlms, Vtranmids = ([zeros(TF,3) for i in 1:m] for j in 1:3)
    end

    # Bound vortices to include in force calculation
    vortices = include_trailingboundvortex ? (1:3) : (2:2)  # 1==ApA, 2==AB, 3==BBp



    if !vlm_vortexsheet || KJforce_type=="regular"

        ## NOTE: Instead of calling the VPM, we use what was calculated
        ## by `solve()`, which includes Rotor-on-VLM velocities
        # Vvpm = Vvpm_on_Xs(pfield, Xs; dt=dt)
        Vvpm = vcat(vlm_system.sol["Vvpm_ApA"], vlm_system.sol["Vvpm_AB"], vlm_system.sol["Vvpm_BBp"])

    elseif vehicle == nothing

        error("KJ force calculation based on vortex sheet has been
                requested, but vehicle has not been provided.")

    else

        # ----- Estimate VPM velocity (including rotors) on vortex sheet

        # Generate function of static particles
        static_particles_function = uns.generate_static_particle_fun(pfield, pfield,
                                        vehicle,
                                        sigma_vlm_surf, sigma_rotor_surf;
                                        vlm_vortexsheet=vlm_vortexsheet,
                                        vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
                                        vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
                                        vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv)

        # Pre-calculate direction of lifting bound vortices
        BVdir = [zeros(TF,3) for i in 1:m]

        for i in 1:m
            BVdir[i] .= B[i]
            BVdir[i] .-= A[i]
            BVdir[i] ./= norm(BVdir[i])
        end

        # Add static particles representing the vortex sheet
        org_np = vpm.get_np(pfield)
        static_particles_function(pfield)

        # Nullify strength of vortex sheet particles to avoid double accounting for VLM velocity
        for P in vpm.iterate(pfield; include_static=true)

            if P.static[1] && 1 <= P.index[1] <= m

                # Save strength for future reference
                for i in 1:3; P.M[i] = P.Gamma[i]; end;

                # Save index before it gets overwritten by the FMM
                P.vol[1] = P.index[1]

                # Zero it out
                P.Gamma .*= 1e-14
            end
        end

        # Evaluate velocity in entire particle field
        vpm._reset_particles(pfield)
        pfield.UJ(pfield)

        # Collect velocities
        weights = zeros(TF,m)
        Vvpm = [zeros(TF,3) for i in 1:m]

        for P in vpm.iterate(pfield; include_static=true)

            # if P.static[1] && 1 < P.index[1] <= m
            if P.static[1] && 1 <= P.vol[1] <= m

                # ind = P.index[1]
                ind = Int(P.vol[1])

                # Calculate weight of this static particle as projection to lifting BV
                weight = KJforce_type=="weighted" ? abs(dot(view(P.M, 1:3), BVdir[ind])) :
                            KJforce_type=="averaged" ? 1.0 :
                            error("Invalid KJforce_type. Options are `\"regular\"`, `\"weighted\"`, or `\"averaged\"`; got $(KJforce_type).")

                # Add weighted velocity
                for i in 1:3
                    Vvpm[ind][i] += weight*P.U[i]
                end
                weights[ind] += weight

            end
        end

        # Divide by total weight to get the final weighted velocity
        Vvpm ./= weights

        Vvpm = vcat(Vvpm, Vvpm, Vvpm)
        # ------
    end


    # Calculate KJ force
    M = zeros(TF, 3)
    for i in 1:m                 # Iterate over horseshoes
        for j in vortices        # Iterate over bound vortices (1==ApA, 2==AB, 3==BBp)

            # Bound vortex' length vector
            if j==1
                l = (A[i]-Ap[i])
            elseif j==2
                l = (B[i]-A[i])
            else
                l = (Bp[i]-B[i])
            end

            # Bound vortex' midpoint
            # X = Xs[i + m*(j-1)]

            # Kinematic velocity at midpoint
            Vtranmid = (Vtran[i + m*(j-1)] + Vtran[i + m*j])/2

            # Effective velocity at midpoint
            V = Vvpm[i + m*(j-1)] + Vvlm[i + m*(j-1)] + Vinfs[i + m*(j-1)] + Vtranmid
            if Vout != nothing
                push!(Vout, [Vvpm[i + m*(j-1)], Vvlm[i + m*(j-1)], Vinfs[i + m*(j-1)], Vtranmid])
            end

            # Circulation
            Gamma = vlm_system.sol["Gamma"][i]

            # Normalization factor
            if per_unit_span
                len = abs((B[i][1]-A[i][1])*spandir[1] + (B[i][2]-A[i][2])*spandir[2] + (B[i][3]-A[i][3])*spandir[3])
            else
                len = 1
            end
            if lenout != nothing
                push!(lenout, per_unit_span==false || len>lencrit ? len : -10*lencrit)
            end

            # Kutta–Joukowski theorem: F = rho * V × vecGamma
            if per_unit_span==false || len > lencrit # NOTE: Use lencrit to avoid dividing by zero
                Ftot[i][1] += rho * Gamma * (V[2]*l[3] - V[3]*l[2]) / len
                Ftot[i][2] += rho * Gamma * (V[3]*l[1] - V[1]*l[3]) / len
                Ftot[i][3] += rho * Gamma * (V[1]*l[2] - V[2]*l[1]) / len

                M += LA.cross(Xs[i] .- cog, Ftot[i])

                if debug && j==2
                    Vtranmids[i] .+= Vtranmid
                    Vvpms[i] .+= Vvpm[i + m*(j-1)]
                    Vvlms[i] .+= Vvlm[i + m*(j-1)]
                end
            end

        end
    end


    if vlm_vortexsheet && KJforce_type!="regular"
        # Remove static particles
        for pi in vpm.get_np(pfield):-1:org_np+1
            vpm.remove_particle(pfield, pi)
        end
    end

    if debug
        # Save Kutta-Joukowski force as a solution field
        vlm._addsolution(vlm_system, (per_unit_span ? "f" : "F")*"rk-vector", deepcopy(Ftot); t=t)
        vlm._addsolution(vlm_system, "Vtranmid-vector", Vtranmids; t=t)
        vlm._addsolution(vlm_system, "Vvpm-vector", Vvpms; t=t)
        vlm._addsolution(vlm_system, "Vvlm-vector", Vvlms; t=t)
    end

    # Add Kutta-Joukowski force to any existing force calculation
    fieldname = per_unit_span ? "ftot" : "Ftot"
    if fieldname in keys(vlm_system.sol)
        Ftot .+= vlm_system.sol[fieldname]
    end

    # Save total force as a solution field
    vlm._addsolution(vlm_system, fieldname, Ftot; t=t)

    return sum(Ftot), M
end

