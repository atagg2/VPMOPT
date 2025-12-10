# scp flowlab@10.32.117.231:/media/flowlab/ExternalStorage1/atagg/fall_2024/VPMOPT/A_3.txt .

# scp flowlab@10.32.117.231:/media/flowlab/ExternalStorage1/atagg/fall_2024/VPMOPT/x_dot_3.txt .

# scp flowlab@10.32.117.231:/media/flowlab/ExternalStorage1/atagg/fall_2024/VPMOPT/aero_3.txt .

# scp flowlab@10.32.117.231:/media/flowlab/ExternalStorage1/atagg/fall_2024/VPMOPT/thrust_3.txt .

# file = "_4_7_25_scaled_variables.txt"
# file2 = "_4_7_25_scaled_variables_cont.txt"
file = "_12_08_25.txt"
#file2 = "_7_28_25.txt"
#file3 = "_8_4_25.txt"

using DelimitedFiles
using Plots
import LinearAlgebra as LA
import FLOWMath as FM

A_2d = Float64.(readdlm("A"*file, ' ')[:,1:end-1])
B_2d = Float64.(readdlm("B"*file, ' ')[:,1:end-2])
x_dot_2d = Float64.(readdlm("x_dot"*file, ' ')[:,1:end-1])
# force = Float64.(readdlm("force"*file, ' ')[:,1:end-1])
x_hist_2d = Float64.(readdlm("x_hist"*file, ' ')[:,1:end-1])
u_hist_2d = Float64.(readdlm("u_hist"*file, ' ')[:,1:end-1])
# step2 = Float64.(readdlm("step"*file, ' ')[:,1:end-1])

#A_2d = [A_2d; Float64.(readdlm("A"*file2, ' ')[:,1:end-1])]
#x_dot_2d = [x_dot_2d; Float64.(readdlm("x_dot"*file2, ' ')[:,1:end-1])]
#x_hist_2d = [x_hist_2d; Float64.(readdlm("x_hist"*file2, ' ')[:,1:end-1])]
#u_hist_2d = [u_hist_2d; Float64.(readdlm("u_hist"*file2, ' ')[:,1:end-1])]


#A_2d = [A_2d; Float64.(readdlm("A"*file3, ' ')[:,1:end-1])]
#x_dot_2d = [x_dot_2d; Float64.(readdlm("x_dot"*file3, ' ')[:,1:end-1])]
#x_hist_2d = [x_hist_2d; Float64.(readdlm("x_hist"*file3, ' ')[:,1:end-1])]
#u_hist_2d = [u_hist_2d; Float64.(readdlm("u_hist"*file3, ' ')[:,1:end-1])]

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


tf = 10.0
nt = 400
ts = range(0.0, tf, nt)

dt = ts[2] - ts[1]

dt_vpm = 0.003 #30/(4*5400)
n_steps_vpm = Int(round(tf/dt_vpm))


nl = 15#Int(round(n_steps_vpm/300))
ns = 1

t_rng = collect(range(dt_vpm*100, tf-dt_vpm*ns, nl))

k = 1
t_sparse = zeros(nl*ns)
for i in 1:nl
    for j in 1:ns
        global k
        t_sparse[k] = t_rng[i] + dt_vpm*(j-1)
        k += 1
    end
end


ts_fine = range(t_sparse[1], t_sparse[end], length(ts))




nx = 6
nu = 7

ns = length(A_2d[1,:])

ni = Int(length(A_2d[:,1])/36)

A_3d = zeros(nx, nx, ns, ni)
B_3d = zeros(nx, nu, ns-1, ni)

k = 1
for i1 in 1:ni
    for i in 1:nx
        for j in 1:nx
            global k
            A_3d[i,j,:,i1] = A_2d[k,:]
            k += 1
        end
    end
end

k = 1
for i1 in 1:ni
    for i in 1:nx
        for j in 1:nu
           global k
            B_3d[i,j,:,i1] = B_2d[k,:]
            k += 1
        end
    end
end

nt = length(x_hist_2d[1,:])
n_iter = Int(length(x_hist_2d[:,1])/nx)

x_hist = zeros(nx, nt, n_iter)
u_hist = zeros(nu, nt-1, n_iter)
x_dot = zeros(nx, length(x_dot_2d[1,:]), n_iter-1)

dynamics = zeros(6, nt-1, n_iter-2)
xdots = zeros(nx, nt, n_iter)

k = 1
for it in 1:n_iter
    for i in 1:nx
        global k
        x_hist[i,:,it] = x_hist_2d[k,:]
        k += 1
    end
end

k = 1
for it in 1:n_iter
    for i in 1:nu
        global k
        u_hist[i,:,it] = u_hist_2d[k,:]
        k += 1
    end
end

k = 1
for it in 1:n_iter-2
    for i in 1:nx
        global k
        x_dot[i,:,it] = x_dot_2d[k,:]
        k += 1
    end
end

using HDF5

x_histf = h5read("/Users/atagg/Documents/research/trajectory-optimization-paper/optimal_trajectory/x_hist.h5", "x_hist")
u_histf = h5read("/Users/atagg/Documents/research/trajectory-optimization-paper/optimal_trajectory/u_hist.h5", "u_hist")

x_hist = cat(x_hist, x_histf, dims=3)
u_hist = cat(u_hist, u_histf, dims=3)

n_iter += length(x_histf[1,1,:])

steps = zeros(n_iter-1)
step_indices = [5,6]#[1,2,3,4,5,6]

objs = zeros(n_iter)

t = range(0.0, 10.0, nt-1)

x_hist_norm = deepcopy(x_hist)
u_hist_norm = deepcopy(u_hist)
# for i in 1:length(x_hist[:,1,1])
#     for j in 1:length(x_hist[1,:,1])
#         x_hist_norm[i,j,:] /= 

#     end
# end
# for i in 1:length(u_hist[:,1,1])
#     for j in 1:length(u_hist[1,:,1])
#         u_hist_norm[i,j,:] /= 
#     end
# end

indices = zeros(length(steps))

px_max = zeros(n_iter-1, 3)

if file == "_5_1_25.txt"
    l_relax = 0.01
    u_relax = 0.3
else
    l_relax = 0.03
    u_relax = 0.3
end

for i in 1:n_iter

    for j in 1:length(u_hist[1,:,i])-1
        objs[i] += abs(u_hist[1,j+1,i] - u_hist[1,j,i]) + abs(u_hist[2,j+1,i] - u_hist[2,j,i]) + abs(u_hist[3,j+1,i] - u_hist[3,j,i])
    end

    # @show objs[i]

    if i > 1

        #f_relax = (u_relax-l_relax)*exp(-0.2*i) + l_relax
        f_relax = 0.03


        px = (x_hist_norm[:,:,i] - x_hist_norm[:,:,i-1])/f_relax
        pu = (u_hist_norm[:,:,i] - u_hist_norm[:,:,i-1])/f_relax


        # px[1,:] /= 4.0
        # px[5,:] /= 15.0
        # px[6,:] /= 2.0
        # pu[1,:] /= 15.0
        # pu[2,:] /= 20.0
        # pu[3,:] /= 25.0

        px_max[i-1, 1], _ = findmax(abs.(px[4,1:end-1]))
        px_max[i-1, 2], _ = findmax(abs.(px[5,1:end-1]))
        px_max[i-1, 3], _ = findmax(abs.(px[6,1:end-1]))


        # px[1,:] /= 33.0
        # px[2,:] /= 4.0
        # px[3,:] /= 3.0
        # px[4,:] /= 250.0
        # px[5,:] /= 20.0
        # px[6,:] /= 10.0
        # pu[1,:] /= 90.0
        # pu[2,:] /= 80.0
        # pu[3,:] /= 180.0
        # steps[i-1], ind = findmax(abs.([px[:,1:end-1]; pu]))
        # indices[i-1] = ind[1]
        # @show steps[i-1] ind
        steps[i-1] = LA.norm([px[:,1:end-1]; pu])
        # steps[i-1] = LA.norm(px[[4,5,6],1:end-1])

        @show steps[i-1]
        # steps[i-1] = LA.norm(px[:,1:end-1])

        # if steps[i-1] < 0.3
        #     stop
        # end

    end


    # p1 = plot(x_hist[4,:,i], x_hist[5,:,i], xlabel = "X (m)", ylabel = "Y (m)", label = false)
    # p2 = plot(t, u_hist[1,:,i], xlabel = "Time (s)", ylabel = "Rotor Velocity (rad/s)", label = "Front", legend = :outertopright)
    # plot!(t, u_hist[2,:,i], label = "Back")
    # plot!(t, u_hist[3,:,i], label = "Pusher")

    # p3 = plot(t, x_hist[6,2:end,i], xlabel = "Time (s)", ylabel = "Pitch (deg)", label = false, ylims = [-0.1, maximum(x_hist[6,:,:]) + 0.1])

    # p_step = plot(p1,p2,p3, layout = [1; 2])

    # display(p_step)


    p1 = plot(x_hist[4,:,i], x_hist[5,:,i], xlabel = "X", ylabel = "Y", label = false)#, ylims = [-0.01, maximum(x_hist[4,:,:])], xlims = [-0.01, maximum(x_hist[5,:,:])])
    p2 = plot(t, x_hist[6,2:end,i], xlabel = "Time (s)", ylabel = "Pitch", label = false, ylims = [minimum(x_hist[6,:,:]), maximum(x_hist[6,:,:])])
    p3 = plot(t, u_hist[1,:,i], xlabel = "Time (s)", ylabel = "Rotor Velocity", label = "Front", legend = :outertopright, ylims = [minimum(u_hist[1:3,:,:]), maximum(u_hist[1:3,:,:])])
    plot!(t, u_hist[2,:,i], label = "Middle")
    plot!(t, u_hist[3,:,i], label = "Back")
    p4 = plot(t, u_hist[4,:,i], xlabel = "Time (s)", ylabel = "Tilt", label = "Front", legend = :outertopright, ylims = [minimum(u_hist[4:end,:,:])-0.1, maximum(u_hist[4:end,:,:])])
    plot!(t, u_hist[5,:,i], label = "Middle")
    plot!(t, u_hist[6,:,i], label = "Back")
    plot!(t, u_hist[7,:,i], label = "Elevator")

    p_step = plot(p1,p2,p3,p4)

    display(p_step)

    sleep(0.1)



    if i < length(x_dot[1,1,:])
        interpolate_xdots(t_sparse, x_dot[:,:,i], ts_fine, xdots[:,:,i], 2)
        dynamics[:,:,i] = (x_hist[:,2:end,i] - x_hist[:,1:end-1,i])/dt - xdots[:,1:end-1,i]
    end
end





xdots_sparse = x_dot[:,:,end-1]
xdots = zeros(nx, nt+1)

As_sparse = A_3d[:,:,:,end]
As = zeros(nx*nt, nx)



interpolate_xdots(t_sparse, xdots_sparse, ts_fine, xdots, 6)

x_dyn = zeros(nx, nt+1)

x_dyn[:,1] = x_hist[:,1,end]


for i in 1:nt
    x_dyn[:,i+1] = x_dyn[:,i] + xdots[:,i]*dt
end

x_scaling = [60.0, 10.0, 1.0, 300.0, 50.0, 5.0]
u_scaling = [400.0, 400.0, 400.0, pi/2, pi/2, pi/2, 25*pi/180]

#f = open("xs_0.txt", "a")
#for i in eachindex(x_hist[:,1,1])
#    for j in eachindex(x_hist[1,:,1])
#        print(f, x_hist[i,j,end]*x_scaling[i], " ")
#    end
#    print(f, "\n")
#end
#close(f)

#f = open("us_0.txt", "a")
#for i in eachindex(u_hist[:,1,1])
#    for j in eachindex(u_hist[1,:,1])
#        print(f, u_hist[i,j,end]*u_scaling[i], " ")
#    end
#    print(f, "\n")
#end
#close(f)

obj_diff = zeros(length(objs)-1)


for i in eachindex(obj_diff)
    obj_diff[i] = objs[i+1] - objs[i]
end

step_diff = zeros(length(steps)-1)

for i in eachindex(step_diff)
    step_diff[i] = steps[i+1] - steps[i]
end

# p = plot()
# for i in 2:length(x_hist[1,1,:])
#     f_relax = (0.3 - 0.1)*exp(-0.2*(i-1)) + 0.1
#     px = (x_hist[:,:,end-i+1] - x_hist[:,:,end-i])/0.1
#     plot!(px[1,:])
#     display(p)
#     sleep(0.1)
# end





# p1 = plot(xs_c[4,:], xs_c[5,:], xlabel = "X", ylabel = "Y", label = false, ylims = [-0.01, maximum(xs_c[4,:])], xlims = [-0.01, maximum(xs_c[5,:])])
# p2 = plot(t, xs_c[6,2:end], xlabel = "Time (s)", ylabel = "Pitch", label = false, ylims = [minimum(xs_c[6,:]), maximum(xs_c[6,:])])
# p3 = plot(t, us_c[1,:], xlabel = "Time (s)", ylabel = "Rotor Velocity", label = "Front", legend = :outertopright, ylims = [minimum(us_c[1:3,:]), maximum(us_c[1:3,:])])
# plot!(t, us_c[2,:], label = "Middle")
# plot!(t, us_c[3,:], label = "Back")
# p4 = plot(t, us_c[4,:], xlabel = "Time (s)", ylabel = "Tilt", label = "Front", legend = :outertopright, ylims = [minimum(us_c[4:end,:,:])-0.1, maximum(us_c[4:end,:])])
# plot!(t, us_c[5,:], label = "Middle")
# plot!(t, us_c[6,:], label = "Back")
# # plot!(t, u_hist[7,:,i], label = "Elevator")

# p_step = plot(p1,p2,p3,p4)


# ts = range(0.0, 10.0, ns)
# nt = 3324

# t = range(0.0, 10.0, nt)


# n_iter = Int(length(x_hist_2d[:,1])/nx)

# A11 = zeros(n_iter, ns)

# k = 1
# for i in 1:n_iter
#     A11[i,:] = A_2d[k,:]
#     k += 1
# end


# A_regressions = zeros(ns, n_iter)

# order = 5
# psi = zeros(ns, order+1)
# for i in 1:ns
#     for j in 1:order+1
#         psi[i,j] = ts[i]^(j-1)
#     end
# end


# p1 = plot(xlabel = "Time (s)", ylabel = "dA_x/dV_z", title = "Data", legend = false)

# for i in [1,8]
#     plot!(ts, A11[i,:], label = "iteration $i")
# end

# p2 = plot(xlabel = "Time (s)", ylabel = "dA_x/dV_z", title = "Regression")

# for i in [1,8]
#     w = (psi' * psi) \ psi' * A11[i,:]
#     A_regressions[:,i] = psi * w
#     plot!(ts, A_regressions[:,i], label = "iteration $i")
# end

# plot(p1, p2)




# samples = zeros(ns, n_iter)
# forces = zeros(nt, n_iter)
# nf = 0
# for i in 1:n_iter
#     forces[:,i] = force[nf+1:nf+nt]
#         for j in 1:ns
#            i_fine = findmin(abs.(ts[j] .- t))[2]
#            samples[j,i] = forces[i_fine,i]
#        end
#     nf += nt
# end

# regressions = zeros(ns, n_iter)

# order = 5
# psi = zeros(ns, order+1)
# for i in 1:ns
#     for j in 1:order+1
#         psi[i,j] = ts[i]^(j-1)
#     end
# end

# for i in 1:n_iter
#     w = (psi' * psi) \ psi' * samples[:,i]
#     regressions[:,i] = psi * w
# end

# p1 = plot(xlabel = "Time (s)", ylabel = "X Acceleration (m/s^2)", title = "Regression", legend=false)
# scatter!(ts, samples[:,1], color = :blue, label = "iteration 1")
# plot!(ts, regressions[:,1], color = :blue, label = false)
# scatter!(ts, samples[:,8], color = :red, label = "iteration 8")
# plot!(ts, regressions[:,8], color = :red, label = false)

# p2 = plot(t, forces[:,1], xlabel = "Time (s)", label = "iteration 1", color = :blue, title = "Data")
# plot!(t, forces[:,8], label = "iteration 8", color = :red)

# plot(p1, p2)

# A = zeros(nx, nx+nu, nt, n_iter)


# nx = 6
# nu = 3
# nt = Int(round(length(A_2d[:,1])/nx))

# A = zeros(nx, nx+nu, nt)

# j = 1
# for i in 1:nx:length(A_2d[:,1])
#     A[:,:,j] = A_2d[i:i+nx-1,:]
#     global j += 1
# end


# # aero = Float64.(readdlm("aero"*file, ' ')[:,1:end-1])

# # T = Float64.(readdlm("thrust"*file, ' ')[:,1:end-1])



# nx = 6
# nu = 3

# n_iter = Int(round(length(x_hist_2d[:,1])/nx))
# @show n_iter
# nt = length(x_hist_2d[1,:])

# x_hist = zeros(nx, nt, n_iter)

# j = 1
# for i in 1:nx:length(x_hist_2d[:,1])
#     x_hist[:,:,j] = x_hist_2d[i:i+nx-1,:]
#     global j += 1
# end

# n_iter = Int(round(length(u_hist_2d[:,1])/nu))
# nt = length(u_hist_2d[1,:])

# u_hist = zeros(nu, nt, n_iter)

# j = 1
# for i in 1:nu:length(u_hist_2d[:,1])
#     u_hist[:,:,j] = u_hist_2d[i:i+nu-1,:]
#     global j += 1
# end

# ts_c = range(0.0, 10.0, nt)

# for i in 1:n_iter

#     xs_step = x_hist[:,:,i]
#     us_step = u_hist[:,:,i]

#     p1 = plot(xs_step[4,:], xs_step[5,:], xlabel = "X (m)", ylabel = "Y (m)", label = false)
#     # p1 = plot(ts_c, xs_step[1,2:end], xlabel = "X (m)", ylabel = "Y (m)", label = false)

#     p2 = plot(ts_c, us_step[1,:], xlabel = "Time (s)", ylabel = "Control Inputs", label = "Front", legend = :outertopright, ylims = [0.0, 200.0])
#     plot!(ts_c, us_step[2,:], label = "Back")
#     plot!(ts_c, us_step[3,:], label = "Pusher")

#     p3 = plot(ts_c, xs_step[6,2:end], xlabel = "Time (s)", ylabel = "Pitch (deg)", label = false)

#     p_step = plot(p1,p2,p3, layout = [1; 2])

#     display(p_step)

#     sleep(0.5)

#     @show i

# end

# nt1 = length(x_dot[:,1])
# t1 = range(0.0, 10.0, nt1)
# dt1 = t1[2] - t1[1]
# xs = zeros(6, nt1)
# for i in 2:32

#     dx = [x_dot[i-1, 1:3]; xs[1:3,i-1]]
#     xs[:,i] = xs[:,i-1] + dx*dt1

# end 



# # time = 0.1

# # As = interpolate_matrix(ts_sparse, As_sparse, ts, As, 7; lambda = 0.0)
# # Bs = interpolate_matrix(ts_sparse, Bs_sparse, ts, Bs, 7; lambda = 0.0)
# # xdots_hist[:,:,iter] = interpolate_xdots(ts_sparse, xdots_sparse, ts, xdots, 7)

# # As_r = zeros(nx, nx, nt)
# # Bs_r = zeros(nx, nu, nt)

# # j = 1
# # for i in 1:nt
# #     Bs_r[:,:,i] = Bs[j:j+nx-1,:]
# #     As_r[:,:,i] = As[j:j+nx-1,:]
# #     j += nx
# # end 

# # pas = []
# # for i in 1:3
# #     p = plot()
# #     for j in 1:6
# #         plot!(ts, As_r[i,j,:], label = "A[$i,$j]", legend = :outertopright)
# #         scatter!(ts_sparse, As_sparse[i,j,:], label = false)
# #         display(p)
# #         sleep(time)
# #     end
# #     push!(pas, p)
# # end

# # pbs = []
# # for i in 1:3
# #     p = plot()
# #     for j in 1:3
# #         plot!(ts, Bs_r[i,j,:], label = "B[$i,$j]", legend = :outertopright)
# #         scatter!(ts_sparse, Bs_sparse[i,j,:], label = false)
# #         display(p)
# #         sleep(time)
# #     end
# #     push!(pbs, p)
# # end

# # pxd = plot()
# # for i in 1:6
# #     plot!(ts, xdots_hist[i,1:end-1,iter], label = "x[$i]", legend = :outertopright)
# #     scatter!(ts_sparse, xdots_sparse[i,:], label = false)
# #     display(pxd)
# #     sleep(time)
# # end
