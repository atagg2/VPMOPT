include("vehicle.jl")
include("vpm_model.jl")
include("collocation.jl")

import PreallocationTools as PT
using DelimitedFiles

x0 = zeros(6); x0[1] = 1e-5; x0[2] = 0.1
u0 = [200.0, 200.0, 200.0, pi/2, pi/2, pi/2, 0.0]

x_scaling = [60.0, 10.0, 1.0, 300.0, 50.0, 5.0]
u_scaling = [400.0, 400.0, 400.0, pi/2, pi/2, pi/2, 25*pi/180]

x0 ./= x_scaling
u0 ./= u_scaling

# define desired final state
xf = [60.0, 0.0, 0.0, 300.0, 50.0, 3.0] ./ x_scaling
Cf = [1.0 0.0 0.0 0.0 0.0 0.0
      0.0 1.0 0.0 0.0 0.0 0.0
      0.0 0.0 0.0 0.0 1.0 0.0
      0.0 0.0 0.0 0.0 0.0 1.0]


nx = length(x0)
nu = length(u0)
nt = 400
tf = 7.5
ts = range(0.0, tf, nt)
dt = ts[2] - ts[1]

R = zeros(nu, nu)

# create vpm simulation
Vref = 1.0
RPMref = 1.0
vehicle, r_rotors, dir_rotors, r_surfs = generate_joby_vehicle(Float64);
v_vehicle(t) = zeros(3)                      
angle_vehicle(t) = zeros(3)                  
angle = ((t)->0.0, (t)->0.0, (t)->0.0, (t)->0.0)                                  
RPM = ((t)->0.0, (t)->0.0, (t)->0.0) 
maneuver = uns.KinematicManeuver(angle, RPM, v_vehicle, angle_vehicle)
sim = uns.Simulation(vehicle, maneuver, Vref, RPMref, tf);


m = 1000.0  #confirm
I = 4000.0  #confirm
cog = [0.8075, 0.0, 0.0]
rho = 1.225
mu = 1.81e-5
g = 9.81# wing


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

A = zeros(nx, nx, nl*ns)
B = zeros(nx, nu, nl*ns)
x_dots = zeros(nx, nl*ns)

A0 = [-0.000465696   0.016747  -0.0048591  -0.0  -0.0  -0.228573
-0.471775     -0.258286   0.0196457   0.0   0.0   0.067471
 2.01242      -0.406216  -0.188692    0.0   0.0   0.699355
 0.2           0.0        0.0         0.0   0.0   0.0
 0.0           0.2        0.0         0.0   0.0   0.0
 0.0           0.0        1.0         0.0   0.0   0.0]

B0 = [0.0184839  0.0181258    0.0180737  -0.115951  -0.113708   -0.111225   -0.0
2.35147    2.30591      2.29928     0.032812   0.0321773   0.0314746   0.0
11.3098     5.84273    -13.8255      0.661582   0.276405    1.32398     0.0
0.0        0.0          0.0         0.0        0.0         0.0         0.0
0.0        0.0          0.0         0.0        0.0         0.0         0.0
0.0        0.0          0.0         0.0        0.0         0.0         0.0]

omega_ref = 100.0
RPM_ref = omega_ref*30/pi 

R = 1.524/4
b = 10.97/4

p_per_step      = 3                         # Sheds per time step
shed_starting   = false                     # Whether to shed starting vortex
shed_unsteady   = true                      # Whether to shed vorticity from unsteady loading
unsteady_shedcrit = 0.001                   # Shed unsteady loading whenever circulation
                                            #  fluctuates by more than this ratio

# Regularization of embedded vorticity
sigma_vlm_surf  = b/400                     # VLM-on-VPM smoothing radius
sigma_rotor_surf= R/20                      # Rotor-on-VPM smoothing radius
lambda_vpm      = 2.125                     # VPM core overlap
                                            # VPM smoothing radius
sigma_vpm_overwrite         = lambda_vpm * (2*pi*RPM_ref/60*R + xf[1])*dt_vpm / p_per_step
sigmafactor_vpmonvlm        = 1.0             # Shrink particles by this factor when
                                            #  calculating VPM-on-VLM/Rotor induced velocities
sigmafactor_vpm = 1.0

# Rotor solver
vlm_rlx                     = 0.2           # VLM relaxation <-- this also applied to rotors
hubtiploss_correction       = vlm.hubtiploss_correction_prandtl # Hub and tip correction
vlm_init = false

# Wing solver: actuator surface model (ASM)
vlm_vortexsheet             = false         # Whether to spread the wing surface vorticity as a vortex sheet (activates ASM)
vlm_vortexsheet_overlap     = 2.125         # Overlap of the particles that make the vortex sheet
vlm_vortexsheet_distribution= uns.g_pressure# Distribution of the vortex sheet
# vlm_vortexsheet_sigma_tbv = thickness*chord / 100  # Size of particles in trailing bound vortices
vlm_vortexsheet_sigma_tbv   = sigma_vpm_overwrite
                                            # How many particles to preallocate for the vortex sheet


# Wing solver: force calculation
KJforce_type                = "regular"     # KJ force evaluated at middle of bound vortices_vortexsheet also true)
include_trailingboundvortex = false         # Include trailing bound vortices in force calculations

include_unsteadyforce       = true          # Include unsteady force
add_unsteadyforce           = false         # Whether to add the unsteady force to Ftot or to simply output it

include_parasiticdrag       = true          # Include parasitic-drag force
add_skinfriction            = true          # If false, the parasitic drag is purely parasitic, meaning no skin friction
calc_cd_from_cl             = false         # Whether to calculate cd from cl or effective AOA
wing_polar_file             = "xf-n0012-il-500000-n5.csv"    # Airfoil polar for parasitic drag


# VPM solver
# vpm_integration = vpm.rungekutta3         # VPM temporal integration scheme
vpm_integration = vpm.euler

vpm_viscous     = vpm.Inviscid()            # VPM viscous diffusion scheme
# vpm_viscous   = vpm.CoreSpreading(-1, -1, vpm.zeta_fmm; beta=100.0, itmax=20, tol=1e-1)

vpm_SFS       = vpm.SFS_none              # VPM LES subfilter-scale model
# vpm_SFS         = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
#                                   alpha=0.999, maxC=1.0,
#                                   clippings=[vpm.clipping_backscatter],
#                                   controls=[vpm.control_directional, vpm.control_magnitude])

                                            

omit_shedding = []
sigmfactor_vpm = 1.0


data_cache = PT.GeneralLazyBufferCache((s)->begin
        T = eltype(s)
        vehicle, _, _, _ = generate_joby_vehicle(T)
        v_vehicle(t) = zeros(T, 3)                      
        angle_vehicle(t) = zeros(T, 3)                  
        angle = ((t)->zero(T), (t)->zero(T), (t)->zero(T), (t)->zero(T))                                   
        RPM = ((t)->zero(T), (t)->zero(T), (t)->zero(T)) 
        maneuver = uns.KinematicManeuver(angle, RPM, v_vehicle, angle_vehicle)
        sim = uns.Simulation(vehicle, maneuver, one(T), one(T), tf*one(T))
    
        max_particles = Int(1e5)
        vpm_solver = [
                            (:formulation, vpm.rVPM),
                            (:viscous, vpm.Inviscid()),
                            (:kernel, vpm.gaussianerf),
                            (:UJ, vpm.UJ_fmm),
                            (:SFS, vpm.SFS_none),
                            (:integration, vpm.rungekutta3),
                            (:transposed, true),
                            (:relaxation, vpm.pedrizzetti),
                            (:fmm, vpm.FMM(; p=4, ncrit=50, theta=0.4, nonzero_sigma=false)),
                        ]
    
        pfield = uns.vpm.ParticleField(max_particles, T;  Uinf=t->zeros(T, 3), vpm_solver...)
    
        max_staticp = 4*uns._get_m_static(vehicle)
        staticpfield = vpm.ParticleField(max_staticp, T; Uinf=t->[T(1e-5), 0.0, 0.0],
                                                                        vpm_solver...)
    
    
        static_particles_function = uns.generate_static_particle_fun(pfield, staticpfield,
                                                                    sim.vehicle, sigma_vlm_surf, sigma_rotor_surf;
                                                                    vlm_vortexsheet=vlm_vortexsheet,
                                                                    vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
                                                                    vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
                                                                    vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv,
                                                                    save_path=nothing,
                                                                    run_name=nothing, nsteps_save=0)
    
    
    
        return sim, pfield, staticpfield, static_particles_function
    end)
                                        

sim, pfield, staticpfield, static_particles_function = data_cache[1.0];


model = VPMModel(sim, pfield, data_cache, dt_vpm, m, I, g, rho, cog, 
                    collect(r_rotors), dir_rotors, r_surfs, A, B, x_dots, [1], t_sparse, 
                    vlm_rlx, sigma_vlm_surf, sigma_rotor_surf, hubtiploss_correction, vlm_init, sigmafactor_vpmonvlm, 
                    sigmafactor_vpm, sigma_vpm_overwrite, omit_shedding, x_scaling, u_scaling);


# xs_0 = zeros(nx, nt)
# us_0 = zeros(nu, nt)

# xs_0[1,:] = (xf[1] - x0[1])/tf * ts#range(x0[1], xf[1], nt)
# # xs_0[2,:] = (xf[5]*x_scaling[5] - x0[5]*x_scaling[5])*2/tf^2/x_scaling[2] * ts
# xs_0[4,:] = (xf[1]*x_scaling[1] - x0[1]*x_scaling[1])/2/tf * ts .^2 / x_scaling[4]#range(x0[4], xf[4], nt)
# xs_0[5,:] = (xf[5] - x0[5])/tf^2 * ts .^2 #range(x0[4], xf[5], nt)
# # xs_0[6,2:end] = xf[6] ./(1 .+ exp.(-2*(ts .- ts[end]/2)/2))
# # xs_0[6,2:end] .-= xs_0[6,2] 
# # xs_0[5,2:end] = xf[5] ./(1 .+ exp.(-2*(ts .- ts[end]/2)/2))
# # xs_0[5,2:end] .-= xs_0[5,2] 
# # for i in 2:nt+1
# #     xs_0[2,i] = (xs_0[5,i] - xs_0[5,i-1])/(ts[2]-ts[1])
# #     xs_0[3,i] = (xs_0[6,i] - xs_0[6,i-1])/(ts[2]-ts[1])
# # end
# # xs_0[2,2] = xs_0[2,3]/2
# us_0[1,:] = range(u0[1], u0[1], nt)
# us_0[2,:] = range(u0[2], u0[2], nt)
# us_0[3,:] = range(u0[3], u0[3], nt)

# us_0[4,:] = range(u0[4], 0.0, nt)
# us_0[5,:] = range(u0[5], 0.0, nt)
# us_0[6,:] = range(u0[6], 0.0, nt)
# us_0[7,:] = range(0.0, -1.0*pi/180 / u_scaling[7], nt)



xs_0 = Float64.(readdlm("xs_0_joby.txt", ' ')[1:6,2:end-1])
us_0 = Float64.(readdlm("us_0_joby.txt", ' ')[:,1:end-1])

xs_0 ./= x_scaling
us_0 ./= u_scaling


# x_scale = [33.0, 3.5, 2.7, 247.9, 20.0, 1.2]
# u_scale = [90.6, 76.7, 180.2]

# for i in eachindex(xs_0[1,:])
#     xs_0[:,i] .*= x_scale
# end

# for i in eachindex(us_0[1,:])
#     us_0[:,i] .*= u_scale
# end

u_lims = [0.0 600.0; 0.0 600.0; 0.0 600.0; 0.0 pi/2; 0.0 pi/2; 0.0 pi/2; -0.5 0.5] ./ u_scaling


run_name = "example_7_14_25"
save_path = nothing#"runs/"*run_name*"/"
                    
xs_hist, us_hist, xdots_hist = collocation(model, sim, A0, B0, Cf, x0, u0, xf, ts;
                                            xs_0 = xs_0,
                                            us_0 = us_0,
                                            na = nl,
                                            u_lims=u_lims,
                                            n_steps_vpm=n_steps_vpm,
                                            p_per_step=p_per_step,
                                            max_particles=Int(1e7),
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
                                            vpm_SFS=vpm_SFS,
                                            run_name=run_name,
                                            save_path=save_path)