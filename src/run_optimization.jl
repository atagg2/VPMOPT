include("vehicle.jl")
include("vpm_model.jl")
include("collocation.jl")

x0 = zeros(6); x0[1] = 0.000001
u0 = [122.75403765737931
102.38282567355749
100]

# define desired final state
xf = [70.0, 0.0, 0.0, 0.0, 100.0, 4.5218090718056425]
Cf = [1.0 0.0 0.0 0.0 0.0 0.0
      0.0 1.0 0.0 0.0 0.0 0.0
      0.0 0.0 1.0 0.0 0.0 0.0
      0.0 0.0 0.0 0.0 1.0 0.0]
    #   0.0 0.0 0.0 0.0 0.0 1.0 0.0]

R = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]


nx = length(x0)
nu = length(u0)
nt = 100
tf = 20.0
ts = range(0.0, tf, nt)
dt = ts[2] - ts[1]

# create vpm simulation
Vref = 1.0
RPMref = 1.0
vehicle, r_rotors, dir_rotors = generate_uli_vehicle(Float64)
v_vehicle(t) = zeros(3)                      
angle_vehicle(t) = zeros(3)                  
angle = ()                                  
RPM = ((t)->0.0, (t)->0.0, (t)->0.0) 
maneuver = uns.KinematicManeuver(angle, RPM, v_vehicle, angle_vehicle)
sim = uns.Simulation(vehicle, maneuver, Vref, RPMref, tf)


m = 3600.0    #confirm
I = 1e4          #confirm
cog = [3.631692, 0.0, 1.73257818]   #confirm
rho = 1.225
mu = 1.81e-5
g = 9.81# wing

nl = 10

n_steps_vpm = 5000
dt_vpm = tf/n_steps_vpm

t_sparse = collect(range(0.0, tf, n_steps_vpm))

A = zeros(nx, nx, nl)
B = zeros(nx, nu, nl)
x_dots = zeros(nx, nl)

A0 = [-3.28141e-11  -2.09263e-9   -1.15879e-12  -0.0  -0.0  -0.171217 
2.65775e-9   -0.00694101    7.68805e-10   0.0   0.0   0.0136669     
-8.46226e-11  -0.000628293  -1.31418e-10   0.0   0.0   1.09658e-11  
1.0           0.0           0.0           0.0   0.0   0.0          
0.0           1.0           0.0           0.0   0.0   0.0           
0.0           0.0           1.0           0.0   0.0   0.0]

B0 = [-0.0        -0.0         0.0156611
0.0961976   0.0762956   0.0
0.0102949  -0.0123433  -0.000380198
0.0         0.0         0.0
0.0         0.0         0.0
0.0         0.0         0.0]


omega_ref = 100.0
RPM_ref = omega_ref*30/pi 

R = 1.524
b = 7.7*2

p_per_step      = 1                         # Sheds per time step
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

                                        
omit_shedding = []
sigmfactor_vpm = 1.0


data_cache = PT.GeneralLazyBufferCache((s)->begin
        T = eltype(s)
        vehicle, _, _ = generate_uli_vehicle(T)
        v_vehicle(t) = zeros(T, 3)                      
        angle_vehicle(t) = zeros(T, 3)                  
        angle = ()                                  
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
        staticpfield = vpm.ParticleField(max_staticp, T; Uinf=t->Vinf(Xdummy, t),
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
                                        

sim, pfield, staticpfield, static_particles_function = data_cache[1.0]


model = VPMModel(sim, pfield, data_cache, dt, m, I, g, rho, cog, 
                    collect(r_rotors), dir_rotors, A, B, x_dots, 1, t_sparse, 
                    vlm_rlx, sigma_vlm_surf, sigma_rotor_surf, hubtiploss_correction, vlm_init, sigmafactor_vpmonvlm, 
                    sigmafactor_vpm, sigma_vpm_overwrite, omit_shedding)


xs_hist, us_hist, xdots_hist = collocation(model, sim, A0, B0, Cf, x0, u0, xf, ts;
                                            na = nl,
                                            p_per_step=1,
                                            max_particles=Int(1e5),
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
                                            unsteady_shedcrit=unsteady_shedcrit)