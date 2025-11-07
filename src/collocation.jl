import ForwardDiff as FD
# import FastDifferentiation as FsD
import FLOWMath as FM
import LinearAlgebra as LA
using Plots
using OSQP
using SparseArrays


file = "_5_10_25.txt"

A_file = "A"*file
x_dot_file = "x_dot"*file
x_file = "x_hist"*file
u_file = "u_hist"*file
force_file = "force"*file
step_file = "step"*file

function interpolate_matrix(ts, As, ts_fine, As_fine, order;
                            lambda =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    nx = length(As[:,1,1])
    nu = length(As[1,:,1])
    nt = length(ts)

    As_r = zeros(nx, nu, nt)
    A = zeros(nx, nu)

    psi = zeros(nt, order+1)
    for i in 1:nt
        for j in 1:order+1
            psi[i,j] = ts[i]^(j-1)
        end
    end

    for i in 1:nx
        for j in 1:nu
            w = (psi' * psi + lambda[i]*LA.I) \ psi' * As[i,j,:]
            As_r[i,j,:] = psi * w
        end
    end

    A_spline = (t)->begin
        for i in eachindex(A[:,1])
            for j in eachindex(A[1,:])
                A[i,j] = FM.akima(ts, As_r[i,j,:],t)
            end
        end
        return A
    end

    j = 1
    for i in 1:nx:nx*length(ts_fine)
        As_fine[i:i+nx-1,:] = A_spline(ts_fine[j])
        j += 1
    end

    return As_fine
end


function interpolate_xdots(ts, xdots, ts_fine, xdots_fine, order)

    nx = length(xdots[:,1])
    nt = length(ts)

    xdots_r = zeros(size(xdots))

    psi = zeros(nt, order+1)
    for i in 1:nt
        for j in 1:order+1
            psi[i,j] = ts[i]^(j-1)
        end
    end

    for i in 1:nx
        w = (psi' * psi) \ psi' * xdots[i,:]
        xdots_r[i,:] = psi * w
    end

    xdot = zeros(nx)
    xdot_spline = (t)->begin
        for i in 1:nx
            xdot[i] = FM.akima(ts, xdots_r[i,:], t)
        end
        return xdot
    end

    for i in eachindex(ts_fine)
        xdots_fine[:,i] = xdot_spline(ts_fine[i])
    end

    return xdots_fine
end



# function interpolate_B(ts, Bs, nx, ts_fine, Bs_fine)

#     nu = length(Bs[1,:])
    # nt = length(ts)

    # Bs_r = zeros(nx, nu, nt)
    # B = zeros(nx,nu)

    # j = 1
    # for i in 1:nt
    #     Bs_r[:,:,i] = Bs[j:j+nx-1,:]
    #     j += nx
    # end

#     B_spline = (t)->begin
#         for i in eachindex(B[:,1])
#             for j in eachindex(B[1,:])
#                 B[i,j] = FM.akima(ts, Bs_r[i,j,:],t)
#             end
#         end
#         return B
#     end

#     j = 1
#     for i in 1:nx:nx*length(ts_fine)
#         Bs_fine[i:i+nx-1,:] = B_spline(ts_fine[j])
#         j += 1
#     end

#     return Bs_fine
# end


function collocation(model, sim, A0, B0, Cf, x0, u0, xf, ts;
    max_iter = 100,
    xs_0 = nothing,
    us_0 = nothing,
    step_indices = [4,5,6],
    verbose = true,
    relax_factor = 0.01,
    relax_step =  1.0,
    relax_index = [5],
    u_lims = [8e-5 600.0; 8e-5 600.0; 8e-5 600.0],
    x_lims = [8e-5 200.0; 8e-5 200.0],
    x_scaling = ones(length(x0)),
    u_scaling = ones(length(u0)),
    g_scaling = ones(length(x0)),
    o_scaling = 1.0,
    na = length(ts),
    n_steps_vpm = 5000,
    rho = 1.225,
    mu = 1.81e-5,
    p_per_step=5,
    max_particles=Int(1e7),
    sigma_vlm_surf=1,
    sigma_rotor_surf=1,
    sigma_vpm_overwrite=1,
    sigmafactor_vpmonvlm=1,
    vlm_vortexsheet=1,
    vlm_vortexsheet_overlap=1,
    vlm_vortexsheet_distribution=1,
    vlm_vortexsheet_sigma_tbv=1,
    vlm_rlx=1,
    hubtiploss_correction=1,
    shed_starting=1,
    shed_unsteady=1,
    unsteady_shedcrit=1,
    vpm_SFS = vpm.SFS_none,
    # ----- OUTPUT OPTIONS ------------------
    save_path="runs/example_run_3/",
    run_name="example_run_3")

    f = open(x_file, "w")
    print(f, "")
    close(f)

    f = open(u_file, "w")
    print(f, "")
    close(f)

    f = open(A_file, "w")
    print(f, "")
    close(f)

    f = open(x_dot_file, "w")
    print(f, "")
    close(f)

    f = open(step_file, "w")
    print(f, "")
    close(f)

    Cx = [1.0  0.0  0.0  0.0  0.0  0.0
          0.0  0.0  0.0  0.0  1.0  0.0]

    nt = length(ts)
    nx = length(x0)
    nu = length(u0)
    xs = zeros(nx, nt+1)
    us = zeros(nu, nt)
    xdots = zeros(nx, nt+1)
    xdots_sparse = zeros(nx, na)
    As = zeros(nx*nt, nx)
    Bs = zeros(nx*nt, nu)
    As_sparse = zeros(nx*na, nx)
    Bs_sparse = zeros(nx*na, nu)
    xs_hist = zeros(nx, nt+1, max_iter)
    us_hist = zeros(nu, nt, max_iter)
    xdots_hist = zeros(nx, nt+1, max_iter)
    xs[:,1] = x0

    x0_unscaled = x0 .* model.x_scaling

    As_d1 = zeros(size(model.A))
    Bs_d1 = zeros((size(model.B)))
    xdots_d1 = zeros(size(model.x_dots))

    px_d1 = zeros(size(xs_hist[:,:,1]))
    pu_d1 = zeros(size(us_hist[:,:,1]))

    px_unscaled = zeros(size(xs_hist[:,:,1]))
    pu_unscaled = zeros(size(us_hist[:,:,1]))

    # indices = Int.(round.(range(1, nt, na-1))) 
    # ts_sparse = ts[indices]

    ts_sparse = model.t_sparse

    j = 1
    for i in 2:nt+1
        xs[:,i] = x0
        us[:,i-1] = u0

        # xdots[:,i] = zeros(length(x0))
        # xdots[:,i-1] = (xs[:,i] - xs[:,i-1])/dt
        As[j:j+nx-1, :] = A0 #A(xs[:,i-1], us[:,i-1])
        Bs[j:j+nx-1, :] = B0 #B(xs[:,i-1], us[:,i-1])
        j += nx
    end

    # As = interpolate_matrix(ts_sparse, As_sparse, ts, As, 3)
    # Bs = interpolate_matrix(ts_sparse, Bs_sparse, ts, Bs, 3)
    # xdots = interpolate_xdots(ts_sparse, xdots_sparse, ts, xdots, 3)

    xs_hist[:,:,1] = xs
    us_hist[:,:,1] = us
    xdots_hist[:,:,1] = xdots

    iter = 2
    step = 100
    step_d1 = step
    step_size = 0.1
    l_relax = 0.01
    u_relax = 0.3
    beta = 0.0
    while iter < max_iter && step >= 0.3 #&& step_size < 1.0

        l = ones(3)*400
        q = [1.0, 1.0, 100.0, 0.0, 0.0, 100.0]*0

        Q = [1.0 0.0 0.0 0.0 0.0 0.0 
             0.0 1.0 0.0 0.0 0.0 0.0 
             0.0 0.0 1.0 0.0 0.0 0.0
             0.0 0.0 0.0 1.0 0.0 0.0 
             0.0 0.0 0.0 0.0 1.0 0.0 
             0.0 0.0 0.0 0.0 0.0 1.0]*0
        R = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 1.0]*0

        ts_c, xs_c, us_c = solve_QP(As, Bs, Q, R, Cf, x0, u0, xf, ts[end], xs_hist[:,:,iter-1]', us_hist[:,:,iter-1]', xdots_hist[:,:,iter-1]'; lambda = l, q=q, u_lims = u_lims, x_lims = x_lims, Cx = Cx, x_scaling=x_scaling, u_scaling=u_scaling, g_scaling=g_scaling, o_scaling=o_scaling)
        # ts_c, xs_c, us_c = solve_QP(As, Bs, R, Cf, x0, u0, xf, ts[end], xs_hist[:,:,iter-1]', us_hist[:,:,iter-1]', xdots_hist[:,:,iter-1]'; lambda = l, u_lims = u_lims, x_scaling=ones(7), u_scaling=ones(3), g_scaling=ones(7), o_scaling=1.0)

        px = xs_c .* x_scaling - xs_hist[:,:,iter-1]
        pu = us_c .* u_scaling - us_hist[:,:,iter-1]
        # step, _ = findmax(abs.(px[step_indices,:]))

        # x_scaling = [33.0, 3.5, 2.7, 247.9, 20.0, 1.2]
        # u_scaling = [90.6, 76.7, 180.2]

        # px[1,:] /= 33.0
        # px[2,:] /= 3.5
        # px[3,:] /= 2.7
        # px[4,:] /= 247.9
        # px[5,:] /= 20.0
        # px[6,:] /= 1.2
        # pu[1,:] /= 90.6
        # pu[2,:] /= 76.7
        # pu[3,:] /= 180.2


        f = open(step_file, "a")
        for i in eachindex(px)
            print(f, px[i], " ")
        end
        for i in eachindex(pu)
            print(f, pu[i], " ")
        end
        print(f, "\n")
        close(f)

        step_size = relax_step/findmax(abs.(px[relax_index,:]))[1]

        # px_d1 = beta*px_d1 + (1-beta)*px
        # pu_d1 = beta*pu_d1 + (1-beta)*pu

        # px = px_d1/(1-beta^(iter-1))
        # pu = pu_d1/(1-beta^(iter-1))

        step = LA.norm([px[:,1:end-1]; pu])


        # px[1,:] *= 33.0
        # px[2,:] *= 3.5
        # px[3,:] *= 2.7
        # px[4,:] *= 247.9
        # px[5,:] *= 20.0
        # px[6,:] *= 1.2
        # pu[1,:] *= 90.6
        # pu[2,:] *= 76.7
        # pu[3,:] *= 180.2


        # f_relax = (u_relax - l_relax)*exp(-0.2*(iter-1)) + l_relax
        f_relax = 0.01


        if verbose @show iter, step, f_relax end


        if f_relax > u_relax f_relax = u_relax end
        if f_relax < l_relax f_relax = l_relax end



        xs_step = xs_hist[:,:,iter-1] + f_relax*px
        us_step = us_hist[:,:,iter-1] + f_relax*pu

        if iter < 3 && xs_0 !== nothing
            xs_step = xs_0
            us_step = us_0
            # px_d1 = zeros(size(xs_hist[:,:,1]))
            # pu_d1 = zeros(size(us_hist[:,:,1]))
        end

        p1 = plot(xs_step[4,:], xs_step[5,:], xlabel = "X (m)", ylabel = "Y (m)", label = false)
        p2 = plot(ts_c, us_step[1,:], xlabel = "Time (s)", ylabel = "Control Inputs", label = "Front", legend = :outertopright)
        plot!(ts_c, us_step[2,:], label = "Back")
        plot!(ts_c, us_step[3,:], label = "Pusher")

        p3 = plot(ts_c, xs_step[6,2:end], xlabel = "Time (s)", ylabel = "Pitch (deg)", label = false)

        p_step = plot(p1,p2,p3, layout = [1; 2])

        display(p_step)

        xs_hist[:,:,iter] = xs_step
        us_hist[:,:,iter] = us_step

        f = open(x_file, "a")
        for i in eachindex(xs_step[:,1])
            for j in eachindex(xs_step[1,:])
                print(f, xs_step[i,j], " ")
            end
            print(f, "\n")
        end
        close(f)

        f = open(u_file, "a")
        for i in eachindex(us_step[:,1])
            for j in eachindex(us_step[1,:])
                print(f, us_step[i,j], " ")
            end
            print(f, "\n")
        end
        close(f)

        angle = ()
        RPM = (t->FM.linear(ts, us_step[1,:]*model.u_scaling[1], t*ts[end])*30/pi, 
               t->FM.linear(ts, us_step[2,:]*model.u_scaling[2], t*ts[end])*30/pi, 
               t->FM.linear(ts, us_step[3,:]*model.u_scaling[3], t*ts[end])*30/pi)

        v_vehicle = t->[-FM.linear(ts, xs_step[1,1:end-1]*model.x_scaling[1], t*ts[end])
                        0.0
                        FM.linear(ts, xs_step[2,1:end-1]*model.x_scaling[2], t*ts[end])]

        angle_vehicle = t->[0.0
                            FM.linear(ts, xs_step[6,1:end-1]*model.x_scaling[6], t*ts[end])
                            0.0]

        # RPM_hover = (t->u0[1]*30/pi,
        #              t->u0[2]*30/pi,
        #              t->u0[3]*30/pi)

        # v_vehicle_hover = t->[-x0[1], 0.0, 0.0]
        # angle_vehicle_hover = t->[0.0, x0[6], 0.0]

        sim.maneuver = uns.KinematicManeuver(angle, RPM, v_vehicle, angle_vehicle)

        sim.nt = -1
        sim.t = 0.0

        sim.vehicle, _, _, _ = generate_uli_vehicle(Float64)

        min_gamma = 0.001
        max_gamma = 2.5
        wake_treatment_strength = uns.remove_particles_strength(min_gamma^2, max_gamma^2; every_nsteps=1)

        minsigma = sigma_vpm_overwrite*0.1
        maxsigma = sigma_vpm_overwrite*5
        wake_treatment_sigma = uns.remove_particles_sigma(minsigma, maxsigma; every_nsteps=1)

        # wake_treatment_box = uns.remove_particles_box([-1.5, -10.0, -0.3], [13.64, 10.0, 5.0], 1)
        # wake_treatment_box = uns.remove_particles_box([-1.0, -3.0, -0.7], [3.0, 3.0, 0.7], 1)
        wake_treatment_box = uns.remove_particles_box([-4.0, -10.0, -2.0], [6.5, 10.0, 2.0], 1)

        # extra_runtime_function = uns.concatenate(wake_treatment_box, wake_treatment_strength, wake_treatment_sigma, model)
        extra_runtime_function = uns.concatenate(wake_treatment_box, wake_treatment_strength, model)

        # extra_runtime_function(args...; optargs...) = wake_treatment(args...; optargs...)||model(args...; optargs...)

        uns.set_V(sim.vehicle, [-x0_unscaled[1], 0.0, 0.0])
        uns.set_W(sim.vehicle, [0.0, x0_unscaled[3], 0.0])

        Oaxis = [cosd(-x0_unscaled[6]) 0.0 sind(-x0_unscaled[6]); 0.0 1.0 0.0; -sind(-x0_unscaled[6]) 0.0 cosd(-x0_unscaled[6])]


        uns.vlm.setcoordsystem(sim.vehicle.system, sim.vehicle.system.O, Oaxis)

        model.i .= 1

        uns.run_simulation(sim, n_steps_vpm;
                        Vinf=(x,t)->0*[5.0, 0.0, -5.0]*exp(-100*t),
                        rho=rho, mu=mu,
                        # ----- SOLVERS OPTIONS ----------------
                        p_per_step=p_per_step,
                        max_particles=max_particles,
                        max_static_particles=nothing,
                        vpm_integration=vpm.euler,
                        vpm_viscous=vpm.Inviscid(),
                        vpm_SFS=vpm_SFS,
                        sigma_vlm_surf=sigma_vlm_surf,
                        sigma_rotor_surf=sigma_rotor_surf,
                        sigma_vpm_overwrite=sigma_vpm_overwrite,
                        sigmafactor_vpmonvlm=sigmafactor_vpmonvlm,
                        vlm_vortexsheet=vlm_vortexsheet,
                        vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
                        vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
                        vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv,
                        vlm_rlx=vlm_rlx,
                        hubtiploss_correction=hubtiploss_correction,
                        shed_starting=shed_starting,
                        shed_unsteady=shed_unsteady,
                        unsteady_shedcrit=unsteady_shedcrit,
                        extra_runtime_function=extra_runtime_function,
                        # ----- OUTPUT OPTIONS ------------------
                        save_path=save_path,
                        run_name=run_name)


        As_sparse = model.A
        Bs_sparse = model.B
        xdots_sparse = model.x_dots


        f = open(A_file, "a")
        for i in eachindex(As_sparse[:,1,1])
            for j in eachindex(As_sparse[1,:,1])
                for it in eachindex(As_sparse[1,1,:])
                    print(f, As_sparse[i,j,it], " ")
                end
                print(f, "\n")
            end
            print(f, "\n")
        end
        close(f)
        f = open(x_dot_file, "a")
        for i in eachindex(xdots_sparse[:,1])
            for it in eachindex(xdots_sparse[1,:])
                print(f, xdots_sparse[i,it], " ")
            end
            print(f, "\n")
        end
        close(f)
        

        As_d1 = beta*As_d1 + (1-beta)*As_sparse
        Bs_d1 = beta*Bs_d1 + (1-beta)*Bs_sparse
        xdots_d1 = beta*xdots_d1 + (1-beta)*xdots_sparse

        As_sparse = As_d1/(1-beta^(iter-1))
        Bs_sparse = Bs_d1/(1-beta^(iter-1))
        xdots_sparse = xdots_d1/(1-beta^(iter-1))



        # if !isnan(xs_step[1,end])
        #     j = 1
        #     k = 1
        #     for i in 2:nt+1
        #         xs[:,i] = x0
        #         us[:,i-1] = u0

        #         if i-1 == indices[k]
        #             xdots_sparse[:,k] =  f(xs_step[:,i-1], us_step[:,i-1])
        #             # xdots[:,i-1] = (xs[:,i] - xs[:,i-1])/dt
        #             As_sparse[j:j+nx-1, :] = A(xs_step[:,i-1], us_step[:,i-1]) #A(xs[:,i-1], us[:,i-1])
        #             Bs_sparse[j:j+nx-1, :] = B(xs_step[:,i-1], us_step[:,i-1]) #B(xs[:,i-1], us[:,i-1])
        #             k += 1
        #             j += nx
        #         end
        #         @show i
        #     end
        # end

        ts_fine = range(ts_sparse[1], ts_sparse[end], length(ts))

        As = interpolate_matrix(ts_sparse, As_sparse, ts_fine, As, 5)
        Bs = interpolate_matrix(ts_sparse, Bs_sparse, ts_fine, Bs, 5)
        xdots_hist[:,:,iter] = interpolate_xdots(ts_sparse, xdots_sparse, ts, xdots, 5)

        As[:,4:5] .= 0.0



        iter += 1
        step_d1 = deepcopy(step)
    end
    return xs_hist[:,:,1:iter-1], us_hist[:,:,1:iter-1], xdots_hist[:,:,1:iter-1]
end





function solve_QP(As, Bs, Q, R, Cf, x0, u0, xf, tf, x_0, u_0, xdot_0;
    lambda = 0.0,
    q = 0.0,
    u_lims = [-Inf*ones(length(u0)) Inf*ones(length(u0))],
    x_lims = [-Inf*ones(length(x0)) Inf*ones(length(x0))],
    Cx = [1.0 0.0 0.0 0.0 0.0 0.0;
          0.0 1.0 0.0 0.0 0.0 0.0],
    x_scaling = ones(length(x0)),
    u_scaling = ones(length(u0)),
    g_scaling = ones(length(x0)),
    o_scaling = 1.0)

    nx = length(As[1,:]); nu = length(Bs[1,:])
    nt = Int(length(As[:,1])/nx)
    nx_bar = nt*nx; nu_bar = nu*nt
    dt = tf/(nt-1)
    Q, H, q, l, u = generate_matrices(As, Bs, Q, R, Cf, x0, u0, xf, nx, nu, nx_bar, nu_bar, dt, nt, x_0, u_0, xdot_0, lambda, q, u_lims, x_lims, Cx)

    # sol = OP_osqp(Qs, Hs, qs, ls, us)
    sol = OP_osqp(Q, H, q, l, u)

    x_bar = sol[1:nx_bar]
    u_bar = sol[nx_bar+1:nx_bar+1+nu_bar-1]
    # sigma = sol[end-ni+1:end]
    return range(0.0, tf, length = nt), hcat(x0, reshape(x_bar, nx, nt)), reshape(u_bar, nu, nt)
end

function OP_osqp(Q, H, q, l, u)
    P = sparse(Q)
    A = sparse(H)
    prob = OSQP.Model()
    OSQP.setup!(prob; P=P, q=q, A=A, l=l, u=u, alpha=1, max_iter=20000, eps_abs = 1e-4, eps_rel = 1e-4, scaling = true)
    sol = OSQP.solve!(prob)
    return sol.x
end

function cost_matrix(Q, R, nx, nu, nx_bar, nu_bar, nt, dt, l, q)
    R_bar = zeros(eltype(dt), nu_bar, nu_bar)
    iu = 1
    for i in 1:nt
        if i == 1 || i == nt
            R_bar[iu:iu+nu-1, iu:iu+nu-1] = R + LA.diagm(l)
            if i == nt
                R_bar[iu:iu+nu-1, iu:iu+nu-1] *= 1
            end
        else
            R_bar[iu:iu+nu-1, iu:iu+nu-1] = R + 2*LA.diagm(l)
        end
        if i < nt
            R_bar[iu:iu+nu-1, iu+nu:iu+2*nu-1] = -LA.diagm(l)
            R_bar[iu+nu:iu+2*nu-1, iu:iu+nu-1] = -LA.diagm(l)
        end
        iu += nu
    end
    Q_bar = zeros(eltype(dt), nx_bar, nx_bar)
    ix = 1
    for i in 1:nt
        if i == 1 || i == nt
            Q_bar[ix:ix+nx-1, ix:ix+nx-1] = Q + LA.diagm(q)
            if i == nt
                Q_bar[ix:ix+nx-1, ix:ix+nx-1] *= 1
            end
        else
            Q_bar[ix:ix+nx-1, ix:ix+nx-1] = Q + 2*LA.diagm(q)
        end
        if i < nt
            Q_bar[ix:ix+nx-1, ix+nx:ix+2*nx-1] = -LA.diagm(q)
            Q_bar[ix+nx:ix+2*nx-1, ix:ix+nx-1] = -LA.diagm(q)
        end
        ix += nx
    end
    return [Q_bar zeros(eltype(dt), nx_bar, nu_bar)
            zeros(eltype(dt), nu_bar, nx_bar) R_bar]
end

function equality_constraint_matrix(As, Bs, Cf, nx, nu, nx_bar, nu_bar, dt)
        
    A_bar = zeros(nx_bar, nx_bar)
    B_bar = zeros(nx_bar, nu_bar)

    j = 1
    for i in 1:nx:nx_bar-nx
        A_bar[i+nx:i+2*nx-1,i:i+nx-1] = As[i+nx:i+2*nx-1,:]*dt + LA.I
        B_bar[i:i+nx-1,j:j+nu-1] = Bs[i:i+nx-1,:]*dt
        j += nu
    end
    B_bar[end-nx+1:end,end-nu+1:end] = Bs[end-nx+1:end,:]*dt

    Cf_bar = zeros(eltype(dt), length(Cf[:,1]), nx_bar); Cf_bar[:,end-length(Cf[1,:])+1:end] = Cf

    Ce = [1.0 0.0 0.0 0.0 0.0 0.0
          0.0 1.0 0.0 0.0 0.0 0.0
          0.0 0.0 1.0 0.0 0.0 0.0]
    Ce_bar = zeros(eltype(dt), length(Ce[:,1]), nx_bar)
    Ce_bar[1:length(Ce[:,1]), end-length(Ce[1,:])+1:end] = Ce
    Ce_bar[1:length(Ce[:,1]), end-2*length(Ce[1,:])+1:end-length(Ce[1,:])] = -Ce
    # Ce_bar[length(Ce[:,1])+1:end, end-2*length(Ce[1,:])+1:end-length(Ce[1,:])] = Ce
    # Ce_bar[length(Ce[:,1])+1:end, end-3*length(Ce[1,:])+1:end-2*length(Ce[1,:])] = -Ce


    H = [A_bar - LA.I(nx_bar) B_bar
        Cf_bar zeros(eltype(dt), length(Cf[:,1]), nu_bar)
        Ce_bar zeros(eltype(dt), length(Ce_bar[:,1]), nu_bar)]

    return H
end

function generate_matrices(As, Bs, Q, R, Cf, x0, u0, xf, nx, nu, nx_bar, nu_bar, dt, nt, x_0, u_0, xdot_0, lambda, q, u_lims, x_lims, Cx)
    Q = cost_matrix(Q, R, nx, nu, nx_bar, nu_bar, nt, dt, lambda, q)
    H = equality_constraint_matrix(As, Bs, Cf, nx, nu, nx_bar, nu_bar, dt)
    q = zeros(eltype(dt), nx_bar+nu_bar)
    b = zeros(eltype(dt), nx_bar)
    j = 1
    for i in 1:nx:nx_bar
        b[i:i+nx-1] = xdot_0[j,:]*dt - As[i:i+nx-1,:]*dt*x_0[j,:] - Bs[i:i+nx-1,:]*dt*u_0[j,:]
        j += 1
    end
    b[1:nx] += (As[1:nx,:]*dt + LA.I)*x_0[1,:]
    b = vcat(b , -Cf*xf, zeros(3))

    H_ineq = [zeros(nu_bar, nx_bar) LA.I(nu_bar)]
    H = [H; H_ineq]

    Cx_bar = zeros(eltype(dt), length(Cx[:,1])*nt, nx_bar)
    j = 1
    for i in 1:length(Cx[:,1]):length(Cx[:,1])*nt
        Cx_bar[i:i+length(Cx[:,1])-1,j:j+nx-1] = Cx
        j += nx
    end

    H = [H; Cx_bar zeros(eltype(dt), length(Cx_bar[:,1]), nu_bar)]


    limits = zeros(length(b) + nu_bar + length(Cx_bar[:,1]), 2)
    limits[1:length(b),1] = -b; limits[1:length(b),2] = -b

    # limits[nx_bar + 5,:] = [8e-5, 10.0]

    for i in length(b)+1:nu:length(b)+nu_bar
        limits[i:i+nu-1,:] = u_lims
    end
    for i in length(b)+nu_bar+1:length(Cx[:,1]):length(b)+nu_bar+length(Cx_bar[:,1])
        limits[i:i+length(Cx[:,1])-1,:] = x_lims
    end
    
    return Q, H, q, limits[:,1], limits[:,2]
end 





# num_curves = 6  # Number of line-scatter pairs
# colors = palette(:default, num_curves)  # Get distinct colors

# labels = ["vx", "vz", "θ_dot", "rx", "rz", "θ"]

# p = plot(title = "dAz/d_", legend = :outertopright, xlabel = "Time (s)")
# for i in 1:num_curves
#     y_line = As_r[2,i,:]  # Generate line data
#     y_scatter = A[2,i,:]
    
#     plot!(ts, y_line, label=labels[i], color=colors[i])
#     scatter!(ts, y_scatter, label=false, color=colors[i], marker=:circle)
# end

# display(p)
