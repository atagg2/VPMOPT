import FLOWUnsteady as uns
import PreallocationTools as PT
vpm = uns.vpm
vlm = uns.vlm


file = "_5_10_25.txt"

A_file = "A"*file
x_dot_file = "x_dot"*file
aero_file = "aero"*file
thrust_file = "thrust"*file
force_file = "force"*file

function get_value(x)
    if eltype(x) == Float64
        return x
    else
        return x.value
    end
end

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
    r_surfs::Array{TF, 2}
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

    x_scaling::Array{TF, 1}
    u_scaling::Array{TF, 1}

end

function (m::VPMModel)(s)

    s[1:6] .*= m.x_scaling
    s[7:end] .*= m.u_scaling
    
    sim_dual, pfield_dual, static_pfield_dual, static_particles_function = m.cache[s];

    if m.i[1] < 2
        vlm.setVinf(sim_dual.vehicle.system, (x,t)->[eltype(s)(1e-5), 0.0, 0.0])
    end

    sim_dual.vehicle = update_vehicle!(sim_dual.vehicle, m.sim.vehicle, s);

    update_pfield!(pfield_dual, m.pfield);

    dt = m.dt
    rho = m.rho
    mass = m.mass
    I = m.I
    g = m.g
    cog = m.cog
    r_rotors = m.r_rotors
    dir_rotors = m.dir_rotors

    wake_coupled = true
    sound_spd = 343

    x = s[1:6]
    u = s[7:end]

    vx, vz, theta_dot, rx, rz, theta = x

    Vdiff = [vx, 0.0, -vz] #- sim_dual.vehicle.V*(-1)
    Wdiff = [0.0, theta_dot*pi/180, 0.0] #- sim_dual.vehicle.W

    Vinf = (x,t)->Vdiff + LA.cross(-Wdiff, x-sim_dual.vehicle.system.O)

    sim_dual.nt = m.sim.nt
    sim_dual.t = m.sim.t

    RPM = (t->u[1]*30/pi,
            t->u[2]*30/pi,
            t->u[3]*30/pi)


    sim_dual.maneuver = uns.KinematicManeuver(sim_dual.maneuver.angle, RPM, sim_dual.maneuver.Vvehicle, sim_dual.maneuver.anglevehicle)
    # sim_dual.maneuver = uns.KinematicManeuver(angle, RPM, v_vehicle, angle_vehicle)

    # relax = m.pfield.relaxation != vpm.relaxation_none &&
    # m.pfield.relaxation.nsteps_relax >= 1 &&
    # m.sim.nt>0 && (m.sim.nt%m.pfield.relaxation.nsteps_relax == 0)

    # org_np = vpm.get_np(pfield_dual)

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


    # prev_vlm_system = uns._get_prev_vlm_system(sim_dual.vehicle)



    # # uns.nextstep_kinematic(sim_dual, dt)
    # uns.nextstep_kinematic(sim_dual.vehicle, dt)



    # Solver-specific pre-calculations
    # uns.precalculations(sim_dual, Vinf, pfield_dual, m.sim.t, dt)


    # # Shed semi-infinite wake
    # uns.shed_wake(sim_dual.vehicle, Vinf, pfield_dual, dt, m.sim.nt; t=m.sim.t,
    #             unsteady_shedcrit=-1,
    #             p_per_step=p_per_step, sigmafactor=m.sigmafactor_vpm,
    #             overwrite_sigma=m.sigma_vpm_overwrite,
    #             omit_shedding=m.omit_shedding)

    solve_system(sim_dual, Vinf, pfield_dual, wake_coupled, dt, m.vlm_rlx,
                m.sigma_vlm_surf, m.sigma_rotor_surf, m.rho, sound_spd,
                static_pfield_dual, m.hubtiploss_correction;
                init_sol=m.vlm_init, sigmafactor_vpmonvlm=m.sigmafactor_vpmonvlm,
                debug=false)


    # calculate wing forces 
    prev_vlm_system = deepcopy(sim_dual.vehicle.vlm_system)
    F, M = calc_aerodynamicforce_kuttajoukowski(sim_dual.vehicle.vlm_system, prev_vlm_system, r_surfs,
                                                pfield_dual, Vinf, dt, rho, cog)

    # @show F M
    # @show sim_dual.vehicle.vlm_system.sol["Gamma"]

    # f = open(aero_file, "a")
    # for i in eachindex(F)
    #     print(f, get_value(F[i]), " ")
    # end
    # for i in eachindex(M)
    #     print(f, get_value(M[i]), " ")
    # end
    # print(f, "\n")
    # close(f)
    
    # calculate rotor forces
    # indices = [1,3,5,7,2,4,6,8,9]
    indices = [1,2,3]
    i = 1
    # f = open(thrust_file, "a")
    for (i, rotor_system) in enumerate(sim_dual.vehicle.rotor_systems)
        for rotor in rotor_system
            T, Q = uns.vlm.calc_thrust_torque(rotor)
            if i < 3 T *= 4 end
            # @show T
            # print(f, get_value(T), " ")
            T_vec = T*dir_rotors[:,indices[i]]
            F += sim_dual.vehicle.system.Oaxis' * T_vec
            M += LA.cross(r_rotors[:,indices[i]] - m.cog, T_vec)   
            i += 1
        end
    end
    # print(f, "\n")
    # close(f)

    F[1] *= -1

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
    return x_dot ./ m.x_scaling
end

function (m::VPMModel)(sim, pfield, T, DT, args...; optargs...)
    # if T > 100*DT
    #     stop
    # end
    # global save_forces

    # if save_forces && sim.t > 10*DT
    #     V = sim.vehicle.V
    #     W = sim.vehicle.W
    #     Oaxis = sim.vehicle.system.Oaxis
    #     O = sim.vehicle.system.O

    #     theta = acosd(LA.dot(Oaxis[:,1], [1,0,0]))*sign(Oaxis[3,1])

    #     x = [-V[1], V[3], W[2], -O[1], O[3], theta]
    #     RPM_front = sim.vehicle.rotor_systems[1][1].RPM
    #     RPM_back = sim.vehicle.rotor_systems[2][1].RPM
    #     RPM_pusher = sim.vehicle.rotor_systems[3][1].RPM
    #     u = [RPM_front, RPM_back, RPM_pusher]*pi/30

    #     x_dot = m([x;u])

    #     f = open(force_file, "a")
    #     for i in eachindex(x_dot)
    #         print(f, x_dot[i], " ")
    #     end
    #     print(f, "\n")
    #     close(f)

    #     if sim.t >= sim.ttot
    #         save_forces = false
    #     end
    # end

    if T > m.t_sparse[m.i[1]]

        V = sim.vehicle.V
        W = sim.vehicle.W
        Oaxis = sim.vehicle.system.Oaxis
        O = sim.vehicle.system.O

        theta = acosd(LA.dot(Oaxis[:,1], [1,0,0]))*sign(Oaxis[3,1])

        x = [-V[1], V[3], W[2]*180/pi, -O[1], O[3], theta]
        RPM_front = sim.vehicle.rotor_systems[1][1].RPM
        RPM_back = sim.vehicle.rotor_systems[2][1].RPM
        RPM_pusher = sim.vehicle.rotor_systems[3][1].RPM
        u = [RPM_front, RPM_back, RPM_pusher]*pi/30


        @show x u

        x ./= m.x_scaling
        u ./= m.u_scaling

        update_pfield!(m.pfield, pfield)

        # A = FD.jacobian(m, [x; u], FD.JacobianConfig(m, [x; u], FD.Chunk(length([x; u]))), Val{false}())

        m.x_dots[:,m.i[1]] = m([x;u])

        A = FD.jacobian(m, [x; u])

        m.A[:,:,m.i[1]] = A[:,1:length(x)]
        m.B[:,:,m.i[1]] = A[:,length(x)+1:end]

        # if m.i[1] < 2
        #     @show "clear dynamic history"
        #     f = open(A_file, "w")
        #     print(f, "")
        #     close(f)
        #     f = open(x_dot_file, "w")
        #     print(f, "")
        #     close(f)
        #     f = open(aero_file, "w")
        #     print(f, "")
        #     close(f)
        #     f = open(thrust_file, "w")
        #     print(f, "")
        #     close(f)
        # end

        # f = open(A_file, "a")
        # for i in eachindex(A[:,1])
        #     for j in eachindex(A[1,:])
        #         print(f, A[i,j], " ")
        #     end
        #     print(f, "\n")
        # end
        # close(f)
        # f = open(x_dot_file, "a")
        # for i in eachindex(m.x_dots[:,m.i[1]])
        #     print(f, m.x_dots[i,m.i[1]], " ")
        # end
        # print(f, "\n")
        # close(f)

        m.i .+= 1

    end

    return false
end

function update_vehicle!(vehicle_dual::uns.UVLMVehicle{N, M, TD}, vehicle::uns.UVLMVehicle{N, M, TF}, s) where {N, M, TD, TF}

    x = s[1:6]
    u = s[7:end]
    vx, vz, theta_dot, rx, rz, theta = x

    # update vehicle velocity
    V = [-vx, 0.0, vz]
    W = [0.0, theta_dot, 0.0]

    Oaxis = [cosd(-theta) 0.0 sind(-theta); 0.0 1.0 0.0; -sind(-theta) 0.0 cosd(-theta)]
    O = [-rx, zero(TD), rz]

    update_wing_system!(vehicle_dual.system, vehicle.system, u)

    update_wing_system!(vehicle_dual.vlm_system, vehicle.vlm_system, u)

    update_wing_system!(vehicle_dual.wake_system, vehicle.wake_system, u)

    for i in eachindex(vehicle.rotor_systems)
        for j in eachindex(vehicle.rotor_systems[i])
            update_rotor!(vehicle_dual.rotor_systems[i][j], vehicle.rotor_systems[i][j], u[i])
        end
    end


    # @show vehicle_dual.rotor_systems[1][1]._wingsystem.wings[1]._HSs

    # vehicle_dual.grids .= vehicle.grids
    vehicle_dual.prev_data .= vehicle.prev_data
    # vehicle_dual.grid_O .= vehicle.grid_O

    vlm.setcoordsystem(vehicle_dual.system, TD.(vehicle.system.O), Oaxis; reset=false)

    vehicle_dual.V .= 0.0# vehicle.V
    vehicle_dual.W .= 0.0#vehicle.W


    return vehicle_dual
end


function update_wing_system!(system_dual::vlm.WingSystem, system::vlm.WingSystem, u)

    TD = eltype(u)

    k = 1
    for i in eachindex(system.wings)
        if !(typeof(system.wings[i]) <: vlm.Rotor)
            update_wing!(system_dual.wings[i], system.wings[i], TD)
        else
            update_rotor!(system_dual.wings[i], system.wings[i], u[k])
            k += 1
        end
    end

    # vlm.setcoordsystem(system_dual, TD.(system.O), TD.(system.Oaxis); reset=false)


    # system_dual.Vinf = deepcopy(system.Vinf)
    system_dual.sol = deepcopy(system.sol)
    system_dual.O = TD.(system.O)
    system_dual.Oaxis = TD.(system.Oaxis)
    system_dual.invOaxis = TD.(system.invOaxis)


end

function update_wing!(wing_dual, wing, TD)

    # wing_dual.leftxl = deepcopy(wing.leftxl)
    # wing_dual.lefty = deepcopy(wing.lefty)
    # wing_dual.leftzl = deepcopy(wing.leftzl)
    # wing_dual.leftchord = deepcopy(wing.leftchord)
    # wing_dual.leftchordtwist = deepcopy(wing.leftchordtwist)
    
    wing_dual.m = deepcopy(wing.m)
    wing_dual.O = TD.(wing.O)
    wing_dual.Oaxis = TD.(wing.Oaxis)
    wing_dual.invOaxis = TD.(wing.invOaxis)

    # wing_dual.Vinf = deepcopy(wing.Vinf)
    wing_dual.sol = deepcopy(wing.sol)
    # wing_dual._xlwingdcr = deepcopy(wing._xlwingdcr)
    # wing_dual._xtwingdcr = deepcopy(wing._xtwingdcr)
    # wing_dual._ywingdcr = deepcopy(wing._ywingdcr)
    # wing_dual._zlwingdcr = deepcopy(wing._zlwingdcr)
    # wing_dual._ztwingdcr = deepcopy(wing._ztwingdcr)
    # wing_dual._xm = deepcopy(wing._xm)
    # wing_dual._ym = deepcopy(wing._ym)
    # wing_dual._zm = deepcopy(wing._zm)
    # wing_dual._xn = deepcopy(wing._xn)
    # wing_dual._yn = deepcopy(wing._yn)
    # wing_dual._zn = deepcopy(wing._zn)
    wing_dual._HSs = deepcopy(wing._HSs)

end

function update_rotor!(rotor_dual, rotor, u)

    # rotor_dual.CW = rotor.CW
    # rotor_dual.r = rotor.r
    # rotor_dual.chord = rotor.chord
    # rotor_dual.theta = rotor.theta
    # rotor_dual.LE_x = rotor.LE_x
    # rotor_dual.LE_z = rotor.LE_z
    # rotor_dual.B = rotor.B
    # rotor_dual.airfoils = rotor.airfoils
    # rotor_dual.turbine_flag = rotor.turbine_flag
    rotor_dual.RPM = u*30/pi
    # rotor_dual.hubR = rotor.hubR
    # rotor_dual.rotorR = rotor.rotorR
    # rotor_dual.m = rotor.m
    rotor_dual.sol = deepcopy(rotor.sol)
    # rotor_dual._r = rotor._r
    # rotor_dual._chord = rotor._chord
    # rotor_dual._theta = rotor._theta
    # rotor_dual._LE_x = rotor._LE_x
    # rotor_dual._LE_z = rotor._LE_z
    # rotor_dual._polars = rotor._polars
    # rotor_dual._polarroot = rotor._polarroot
    # rotor_dual._polartip = rotor._polartip


    update_wing_system!(rotor_dual._wingsystem, rotor._wingsystem, u)

end

function update_pfield!(pfield_dual::uns.vpm.ParticleField{TD, FD, VD, TUinfD, SD, TkernelD, TUJD, TintegrationD, TRelaxationD, TGPUD}, 
                        pfield::uns.vpm.ParticleField{TF, F, V, TUinf, S, Tkernel, TUJ, Tintegration, TRelaxation, TGPU}) where {TD, FD, VD, TUinfD, SD, TkernelD, TUJD, TintegrationD, TRelaxationD, TGPUD, 
                                                                                                                                TF, F, V, TUinf, S, Tkernel, TUJ, Tintegration, TRelaxation, TGPU}
    pfield_dual.particles = TD.(pfield.particles)
    pfield_dual.t = deepcopy(pfield.t)
    pfield_dual.nt = deepcopy(pfield.nt)
    pfield_dual.np = deepcopy(pfield.np)

end

function calc_aerodynamicforce_kuttajoukowski(vlm_system::vlm.AbstractWing{TF_design,TF_trajectory},
                                prev_vlm_system, r_surfs, pfield, Vinf, dt, rho, cog; t=0.0,
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

                # if i == 26
                #     @show Gamma Vvpm[i + m*(j-1)] Vvlm[i + m*(j-1)] Vinfs[i + m*(j-1)] Vtranmid l V
                # end

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

    # calculate moments
    M = zeros(TF, 3)
    for i in eachindex(vlm_system.wings)
        surface = vlm_system.wings[i]
        F_surf = surface.sol["Ftot"]
        r_surf = r_surfs[:,i]
        xn, yn, zn = surface._xn, surface._yn, surface._zn
        for j in 1:length(xn)-1
            r = [(xn[j] + xn[j+1])/2,
                 (yn[j] + yn[j+1])/2,
                 (zn[j] + zn[j+1])/2] - cog + r_surf
            M += LA.cross(r, surface.Oaxis*F_surf[j])
        end
    end

    return sum(Ftot), M
end


function solve_system(self::uns.Simulation{V, M, R}, Vinf::Function,
                pfield::uns.vpm.ParticleField, wake_coupled::Bool,
                dt::Real, rlx::Real, sigma_vlm::Real, sigma_rotor::Real,
                rho::Real, speedofsound, staticpfield::uns.vpm.ParticleField,
                hubtiploss_correction;
                mirror=false, mirror_X=nothing, mirror_normal=nothing,
                init_sol::Bool=false, sigmafactor_vpmonvlm=1,
                debug=false
                ) where {V<:uns.UVLMVehicle, M<:uns.AbstractManeuver, R}


    # Bring the rotors back to their initial positions at beginning of this
    # time step before `precalculations`. This is needed to avoid double
    # accounting for the RPM when the solver is run

    # uns.rotate_rotors(self, -dt)

    vhcl = self.vehicle
    mnvr = self.maneuver
    t = self.t

    # TIME-STEPPING PROCEDURE:
    # -1) Solve one particle field time step
    # 0) Translate and rotate systems
    # 1) Recalculate horseshoes with kinematic velocity
    # 2) Paste previous Gamma solution to new system position after translation
    # 3) Shed semi-infinite wake after translation
    # 4) Calculate wake-induced velocity on VLM and Rotor system
    # 5) Solve VLM and Rotor system
    # 6) Shed unsteady-loading wake after new solution
    # 7) Save new solution as prev solution
    # Iterate
    #
    # On the first time step (pfield.nt==0), it only does steps (5) and (7),
    # meaning that the unsteady wake of the first time step is never shed.
   
        # ---------- 4) Calculate VPM velocity on VLM and Rotor system -----
        m = vlm.get_m(vhcl.vlm_system)

        # Control points of VLM
        Xs_cp_vlm = uns._get_Xs(vhcl.vlm_system, "CP"; t=t)

        # Mid-points along VLM filaments for calc_aerodynamicforce()
        ## Nodes of every horseshoe
        Ap = uns._get_Xs(vhcl.vlm_system, "Ap")
        A = uns._get_Xs(vhcl.vlm_system, "A")
        B = uns._get_Xs(vhcl.vlm_system, "B")
        Bp = uns._get_Xs(vhcl.vlm_system, "Bp")
        ## Midpoints of bound vortices
        ApA = (Ap .+ A)/2
        AB = (A .+ B)/2
        BBp = (B .+ Bp)/2
        ## Organize them
        Xs_ApA_AB_BBp_vlm = vcat(ApA, AB, BBp)

        # Points at rotors where to calculate induced velocities
        Xs_rotors = uns._get_midXs(vhcl.rotor_systems)

        # Calculate VPM velocity on all points (VLM and rotors)
        Vvpm = uns.Vvpm_on_Xs(pfield, vcat(Xs_cp_vlm, Xs_ApA_AB_BBp_vlm, Xs_rotors); dt=dt, fsgm=sigmafactor_vpmonvlm,
                         mirror, mirror_X, mirror_normal)

        Vvpm_cp_vlm = Vvpm[1:m]
        Vvpm_ApA_AB_BBp_vlm = Vvpm[m+1:4*m]
        Vvpm_rotors = Vvpm[4*m+1:end]

        # Stitches previous rotor solutions for blade induced velocity
        prev_rotor_systems = uns._get_prev_rotor_systems(vhcl)
        allrotors = vlm.WingSystem()

        for (si, rotors) in enumerate(vhcl.rotor_systems)
            for (ri, rotor) in enumerate(rotors)

                vlm.addwing(allrotors, "S$(si)R$(ri)", rotor; reset=false)
                vlm.getHorseshoe(rotor, 1)          # Force HS calculation
                vlm._addsolution(rotor, "Gamma",    # Give previous solution
                                prev_rotor_systems[si][ri]._wingsystem.sol["Gamma"])
            end
        end

        # Particles for Rotor-on-VLM induced velocity (and Rotor-on-Rotor)
        static_particles_fun2(pfield, args...) = uns._static_particles(pfield, allrotors, sigma_rotor)

        # Evaluate Rotor-on-VLM induced velocity
        Vrotor_on_wing = uns.Vvpm_on_Xs(staticpfield, vcat(Xs_cp_vlm, Xs_ApA_AB_BBp_vlm);
                                    static_particles_fun=static_particles_fun2, dt=dt, fsgm=sigmafactor_vpmonvlm,
                                    mirror, mirror_X, mirror_normal)

        # Add and save VPM-on-VLM and Rotor-on-VLM induced velocity
        vlm._addsolution(vhcl.vlm_system, "Vvpm",
                                         Vvpm_cp_vlm + Vrotor_on_wing[1:m]; t=t)
        vlm._addsolution(vhcl.vlm_system, "Vvpm_ApA",
                        Vvpm_ApA_AB_BBp_vlm[1:m] + Vrotor_on_wing[m+1:2*m]; t=t)
        vlm._addsolution(vhcl.vlm_system, "Vvpm_AB",
                        Vvpm_ApA_AB_BBp_vlm[m+1:2*m] + Vrotor_on_wing[2*m+1:3*m]; t=t)
        vlm._addsolution(vhcl.vlm_system, "Vvpm_BBp",
                        Vvpm_ApA_AB_BBp_vlm[2*m+1:3*m] + Vrotor_on_wing[3*m+1:4*m]; t=t)

        # Calculate induced velocities to use in rotor solver
        ## Particles for VLM-on-Rotor induced velocity
        # NOTE: If I keep the rotor particles wouldn't I be double-counting
        # the blade induced velocity and make the solution unstable?
        # ANSWER: We are never double-accounting for the blade induced velocity
        function static_particles_fun3(pfield, args...)
            # Rotor static particles
            static_particles_fun2(pfield, args...)
            # VLM static particles
            uns._static_particles(pfield, vhcl.vlm_system, sigma_vlm)
        end

        ## Evaluate VLM-on-Rotor and Rotor-on-Rotor induced velocity
        Vvlmrotor_on_rotor = uns.Vvpm_on_Xs(staticpfield, Xs_rotors;
                              static_particles_fun=static_particles_fun3, dt=dt, fsgm=sigmafactor_vpmonvlm,
                              mirror, mirror_X, mirror_normal)

        # Add VPM-on-Rotor, VLM-on-Rotor, and Rotor-on-Rotor induced velocity
        Vinds = Vvpm_rotors #+ Vvlmrotor_on_rotor

        # ---------- 5) Solve VLM system -----------------------------------
        # Wake-coupled solution
        if wake_coupled
            vlm.solve(vhcl.vlm_system, Vinf; t=t, extraVinf=_extraVinf4,
                                    keep_sol=true, vortexsheet=(X,t)->zeros(3))
        # Wake-decoupled solution
        else
            vlm.solve(vhcl.vlm_system, Vinf; t=t, extraVinf=_extraVinf1,
                                                                keep_sol=true)
        end

        # Relaxes (rlx->1) or stiffens (rlx->0) the VLM solution
        if rlx > 0
            rlxd_Gamma = rlx*vhcl.vlm_system.sol["Gamma"] +
                                (1-rlx)*uns._get_prev_vlm_system(vhcl).sol["Gamma"]
            vlm._addsolution(vhcl.vlm_system, "Gamma", rlxd_Gamma)
        end

        # ---------- 5) Solve Rotor system ---------------------------------

        # Solve Rotors
        for (si, rotors) in enumerate(vhcl.rotor_systems)

            # RPM of this rotor system
            RPM = self.RPMref*uns.get_RPM(mnvr, si, t/self.ttot)
            for (ri, rotor) in enumerate(rotors)

                # Calculate kinematic velocities
                Vkin = uns._Vkinematic_rotor(vhcl.rotor_systems,
                                         prev_rotor_systems, si, ri, dt)

                # Get velocities induced at every blade of this rotor
                this_Vinds = uns._parse_midXs(vhcl.rotor_systems, Vinds, si, ri)

                VindVkin = uns._format_blades(this_Vinds, vhcl.rotor_systems,
                                                                        si, ri)
                if wake_coupled
                    vlm.solvefromV(rotor, VindVkin, Vinf, RPM, rho; t=t,
                                    sound_spd=speedofsound,
                                    hubtiploss_correction=hubtiploss_correction,
                                    debug=debug, verbosewarn=false)
                else
                    vlm.solvefromCCBlade(rotor, Vinf, RPM, rho; t=t,
                                            sound_spd=speedofsound,
                                            debug=debug, verbosewarn=false)
                end

                # Relaxes (rlx->1) or stiffens (rlx->0) the rotor solution
                if rlx > 0
                    rlxd_Gamma = rlx*vhcl.rotor_systems[si][ri]._wingsystem.sol["Gamma"] +
                                        (1-rlx)*prev_rotor_systems[si][ri]._wingsystem.sol["Gamma"]
                    vlm._addsolution(rotor, "Gamma", rlxd_Gamma)
                end
            end
        end

        # Rotate rotors
        # NOTE: this means that the rotational velocity hadn't been included in
        # the kinematic velocities up to this point
        # uns.rotate_rotors(self, dt)
end

function _extraVinf4(i, t; wing=nothing)
    if wing==nothing; error("Logic error!"); end;
    return wing.sol["Vvpm"][i]
end