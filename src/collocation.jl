import ForwardDiff as FD
# import FastDifferentiation as FsD
import FLOWMath as FM
import LinearAlgebra as LA
using Plots
using OSQP
using SparseArrays


function interpolate_A(ts, As, ts_fine, As_fine)

    nx = length(As[1,:])
    nt = length(ts)

    As_r = zeros(nx, nx, nt)
    A = zeros(nx,nx)

    j = 1
    for i in 1:nt
        As_r[:,:,i] = As[j:j+nx-1,:]
        j += nx
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


function interpolate_xdots(ts, xdots, ts_fine, xdots_fine)

    nx = length(xdots[:,1])

    xdot = zeros(nx)
    xdot_spline = (t)->begin
        for i in 1:nx
            xdot[i] = FM.akima(ts, xdots[i,:], t)
        end
        return xdot
    end

    for i in eachindex(ts_fine)
        xdots_fine[:,i] = xdot_spline(ts_fine[i])
    end

    return xdots_fine
end



function interpolate_B(ts, Bs, nx, ts_fine, Bs_fine)

    nu = length(Bs[1,:])
    nt = length(ts)

    Bs_r = zeros(nx, nu, nt)
    B = zeros(nx,nu)

    j = 1
    for i in 1:nt
        Bs_r[:,:,i] = Bs[j:j+nx-1,:]
        j += nx
    end

    B_spline = (t)->begin
        for i in eachindex(B[:,1])
            for j in eachindex(B[1,:])
                B[i,j] = FM.akima(ts, Bs_r[i,j,:],t)
            end
        end
        return B
    end

    j = 1
    for i in 1:nx:nx*length(ts_fine)
        Bs_fine[i:i+nx-1,:] = B_spline(ts_fine[j])
        j += 1
    end

    return Bs_fine
end


function collocation(model, sim, A0, B0, Cf, x0, u0, xf, ts;
    max_iter = 20,
    step_indices = [1],
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
    p_per_step=1,
    max_particles=Int(1e5),
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
    # ----- OUTPUT OPTIONS ------------------
    save_path="runs/",
    run_name="example_run")

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

    na = 10
    indices = Int.(round.(range(1, nt, na))) 
    ts_sparse = ts[indices]

    j = 1
    k = 1
    for i in 2:nt+1
        xs[:,i] = x0
        us[:,i-1] = u0

        if i-1 == indices[k]
            xdots_sparse[:,k] = zeros(length(x0))
            # xdots[:,i-1] = (xs[:,i] - xs[:,i-1])/dt
            As_sparse[j:j+nx-1, :] = A0 #A(xs[:,i-1], us[:,i-1])
            Bs_sparse[j:j+nx-1, :] = B0 #B(xs[:,i-1], us[:,i-1])
            k += 1
            j += nx
        end
    end

    As = interpolate_A(ts_sparse, As_sparse, ts, As)
    Bs = interpolate_B(ts_sparse, Bs_sparse, nx, ts, Bs)
    xdots = interpolate_xdots(ts_sparse, xdots_sparse, ts, xdots)

    xs_hist[:,:,1] = xs
    us_hist[:,:,1] = us
    xdots_hist[:,:,1] = xdots

    iter = 2
    step = 100
    step_d1 = step
    step_size = 0.1
    relax = 1.0
    while iter < max_iter && step >= 1.0 #&& step_size < 1.0

        l = ones(3)*10000
        q = [1.0, 1.0, 100.0, 0.0, 0.0, 100.0]*0

        Q = [1.0 0.0 0.0 0.0 0.0 0.0 
             0.0 1.0 0.0 0.0 0.0 0.0 
             0.0 0.0 1.0 0.0 0.0 0.0
             0.0 0.0 0.0 1.0 0.0 0.0 
             0.0 0.0 0.0 0.0 1.0 0.0 
             0.0 0.0 0.0 0.0 0.0 1.0]*0
        R = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 1.0]

        ts_c, xs_c, us_c = solve_QP(As, Bs, Q, R, Cf, x0, u0, xf, ts[end], xs_hist[:,:,iter-1]', us_hist[:,:,iter-1]', xdots_hist[:,:,iter-1]'; lambda = l, q=q, u_lims = u_lims, x_lims = x_lims, Cx = Cx, x_scaling=x_scaling, u_scaling=u_scaling, g_scaling=g_scaling, o_scaling=o_scaling)
        # ts_c, xs_c, us_c = solve_QP(As, Bs, R, Cf, x0, u0, xf, ts[end], xs_hist[:,:,iter-1]', us_hist[:,:,iter-1]', xdots_hist[:,:,iter-1]'; lambda = l, u_lims = u_lims, x_scaling=ones(7), u_scaling=ones(3), g_scaling=ones(7), o_scaling=1.0)

        px = xs_c .* x_scaling - xs_hist[:,:,iter-1]
        pu = us_c .* u_scaling - us_hist[:,:,iter-1]
        step, _ = findmax(abs.(px[step_indices,:]))


        if verbose @show iter, step, step_size end

        step_size = relax_step/findmax(abs.(px[relax_index,:]))[1]

        xs_step = xs_hist[:,:,iter-1] + 1*px
        us_step = us_hist[:,:,iter-1] + 1*pu

        p1 = plot(xs_step[4,:], xs_step[5,:], xlabel = "X (m)", ylabel = "Y (m)", label = false)
        p2 = plot(ts_c, us_step[1,:], xlabel = "Time (s)", ylabel = "Control Inputs", label = "Front", legend = :outertopright)
        plot!(ts_c, us_step[2,:], label = "Back")
        plot!(ts_c, us_step[3,:], label = "Pusher")

        p3 = plot(ts_c, xs_step[6,2:end], xlabel = "Time (s)", ylabel = "Pitch (deg)", label = false)

        p_step = plot(p1,p2,p3, layout = [1; 2])

        display(p_step)

        xs_hist[:,:,iter] = xs_step
        us_hist[:,:,iter] = us_step



        angle = ()
        RPM = (t->FM.akima(ts, us_step[1,:], t)*30/pi, 
               t->FM.akima(ts, us_step[2,:], t)*30/pi, 
               t->FM.akima(ts, us_step[3,:], t)*30/pi)

        v_vehicle = t->[FM.akima(ts, xs_step[1,:], t)
                        0.0
                        FM.akima(ts, xs_step[2,:], t)]

        angle_vehicle = t->[0.0
                            FM.akima(ts, xs_step[6,:], t)
                            0.0]

        sim.maneuver = uns.KinematicManeuver(angle, RPM, v_vehicle, angle_vehicle)

        uns.run_simulation(sim, n_steps_vpm;
                        Vinf=(x,t)->zeros(3),
                        rho=rho, mu=mu,
                        # ----- SOLVERS OPTIONS ----------------
                        p_per_step=p_per_step,
                        max_particles=max_particles,
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
                        extra_runtime_function=model,
                        # ----- OUTPUT OPTIONS ------------------
                        save_path=save_path,
                        run_name=run_name)

        As_sparse = model.A
        Bs_sparse = model.B
        xdots_sparse = model.x_dots

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

        As = interpolate_A(ts_sparse, As_sparse, ts, As)
        Bs = interpolate_B(ts_sparse, Bs_sparse, nx, ts, Bs)
        xdots_hist[:,:,iter] = interpolate_xdots(ts_sparse, xdots_sparse, ts, xdots)

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