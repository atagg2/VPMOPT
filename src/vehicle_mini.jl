import FLOWUnsteady as uns

function generate_uli_vehicle(T)
    mpft = 0.3048/2# general parameters
    m = 200.0  #confirm
    I = 500.0  #confirm
    cog = [0.0, 0.0, 0.0]   #confirm
    rho = 1.225
    mu = 1.81e-5
    g = 9.81# wing
    xle_wing = [0.0, 0.0]*mpft
    yle_wing = [0.0, 25.265]*mpft
    zle_wing = yle_wing*-sind(2.0)
    c_wing = [5.44328, 2.84088]*mpft
    twist_wing = [0.0, 0.0] .+ 3.0*pi/180
    # airfoil = DF.readdlm("files/airfoils/LS(1)-0417.dat")        # check units
    camber_wing = fill((xc) -> 0, length(xle_wing))
    mirror_wing = true
    xle_tail = [0.0, 0.0] .+ 18.5*mpft
    yle_tail = [0.0, 6.52000]*mpft
    zle_tail = [0.0, 0.0]
    c_tail = [3.7875, 2.27250]*mpft           # double check this
    twist_tail = [0.0, 0.0]
    camber_tail = fill((xc) -> 0, length(xle_tail))
    mirror_tail = true
    Rtip = 5.0*mpft

    Ry = (theta)->[cosd(theta) 0.0 sind(theta)
                    0.0 1.0 0.0 
                    -sind(theta) 0.0 cosd(theta)]
    axis = Ry(90)*[-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]
    dir = [1.0, 0.0, 0.0]

    CW = [true, true, true]
        
    axis_p = [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]
    dir_p = [1.0, 0.0, 0.0]

    axes = zeros(3,3,3)
    for i in 1:2
        axes[:,:,i] = axis
    end 
    axes[:,:,end] = axis_p

    system = uns.vlm.WingSystem(; TF_trajectory=T, )


    b_wing = yle_wing[2]*2
    AR_wing = b_wing/c_wing[1]
    ns_wing = 12
    r_wing = [0.0, 1.0]
    clen_wing = c_wing/c_wing[1]
    sweep_wing = [atand(xle_wing[2]-xle_wing[1], yle_wing[2]-yle_wing[1])]
    dihed_wing =[atand(zle_wing[2]-zle_wing[1], yle_wing[2]-yle_wing[1])]
    wing = uns.vlm.complexWing(T(b_wing), T(AR_wing), ns_wing, T.(r_wing), T.(clen_wing), T.(twist_wing*180/pi), T.(sweep_wing), T.(dihed_wing), T)
    uns.vlm.setcoordsystem(wing, [xle_wing[1], yle_wing[1], zle_wing[1]], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    uns.vlm.addwing(system, "wing", wing)

    b_tail = yle_tail[2]*2
    AR_tail = b_tail/c_tail[1]
    ns_tail = 6
    r_tail = [0.0, 1.0]
    clen_tail = c_tail/c_tail[1]
    sweep_tail = [atand(xle_tail[2]-xle_tail[1], yle_tail[2]-yle_tail[1])]
    dihed_tail =[atand(zle_tail[2]-zle_tail[1], yle_tail[2]-yle_tail[1])]
    tail = uns.vlm.complexWing(T(b_tail), T(AR_tail), ns_tail, T.(r_tail), T.(clen_tail), T.(twist_tail*180/pi), T.(sweep_tail), T.(dihed_tail), T)
    uns.vlm.setcoordsystem(tail, [xle_tail[1], yle_tail[1], zle_tail[1]], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    uns.vlm.addwing(system, "tail", tail)

    vlm_system = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(vlm_system, "wing", wing)
    uns.vlm.addwing(vlm_system, "tail", tail)

    omega_ref = 100.0
    RPM_ref = omega_ref*30/pi
    n_ref = omega_ref/2/pi
    V_ref = 50.0
    D_ref = 2*Rtip
    J_ref = V_ref/n_ref/D_ref
    mu = 1.81e-5
    rho = 1.225

    r_surfs = [xle_wing[1] xle_tail[1]
               yle_wing[1] yle_tail[1]
               zle_wing[1] zle_tail[1]]


    pos = [-7.0 0.0 0.0
            11.0 0.0 0.0
            23.0 0.0 0.0]*mpft
    r_rotors = pos'
    dir_rotors = zeros(3,3)
    for i in 1:2
        dir_rotors[:,i] = [0.0, 0.0, 1.0]
    end
    dir_rotors[:,end] = [-1.0, 0.0, 0.0]

    rotor_files = fill("apc10x7_mini.csv", 3)
    data_path = "src/rotor/"
    rotors = Vector{uns.vlm.Rotor}(undef, length(rotor_files))
    for i in eachindex(rotors)
        @show i
        rotors[i] = uns.generate_rotor(rotor_files[i]; TF_trajectory=T, pitch = 0.0, altReD=[RPM_ref, J_ref, mu/rho],
                                        n=7, blade_r=1/20, CW=CW[i], 
                                        xfoil=false, verbose=true, data_path=data_path, plot_disc=false)
        uns.vlm.setcoordsystem(rotors[i], pos[i,:], axes[:,:,i])
        uns.vlm.addwing(system, "rotor$i", rotors[i])
    end

    rotor_systems = (rotors[[1]], rotors[[2]], rotors[[3]])
    

    wake_system = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(wake_system, "wing", wing)
    uns.vlm.addwing(wake_system, "tail", tail)
    for i in eachindex(rotors)
        uns.vlm.addwing(wake_system, "rotor$i", rotors[i])
    end

    # build vehicle
    vehicle = uns.VLMVehicle(system;
            V = zeros(T, 3),
            W = zeros(T, 3),
            rotor_systems=rotor_systems,
            vlm_system=vlm_system,
            wake_system=wake_system)

    return vehicle, r_rotors, dir_rotors, r_surfs
end