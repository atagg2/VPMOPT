import FLOWUnsteady as uns

function generate_uli_vehicle(T)
    mpft = 0.3048# general parameters
    m = 3600.0    #confirm
    I = 1e4          #confirm
    cog = [3.631692, 0.0, 1.73257818]   #confirm
    rho = 1.225
    mu = 1.81e-5
    g = 9.81# wing
    xle_wing = [0.0, 0.0] .+ 8.89*mpft
    yle_wing = [0.0, 25.265]*mpft
    zle_wing = yle_wing*-sind(2.0) .+ 8.5*mpft
    c_wing = [5.44328, 2.84088]*mpft
    twist_wing = [0.0, 0.0] .+ 5.0*pi/180
    # airfoil = DF.readdlm("files/airfoils/LS(1)-0417.dat")        # check units
    camber_wing = fill((xc) -> 0, length(xle_wing))
    mirror_wing = true
    xle_tail = [0.0, 0.0] .+ 27.428*mpft
    yle_tail = [0.0, 6.52000]*mpft
    zle_tail = [0.0, 0.0] .+ 8.008*mpft
    c_tail = [3.7875, 2.27250]*mpft           # double check this
    twist_tail = [0.0, 0.0]
    camber_tail = fill((xc) -> 0, length(xle_tail))
    mirror_tail = true
    Rtip = 5.0*mpft
    Rhub = 0.09*Rtip
    n = 2
    propgeom = [
            0.0909091  0.08      46.75
            0.173554   0.110336  47.0562
            0.256198   0.137292  46.9327
            0.338843   0.16024   46.3793
            0.421487   0.178556  45.3963
            0.504132   0.191612  43.9834
            0.586777   0.198783  42.1409
            0.669422   0.2       39.8685
            0.752066   0.2       37.1664
            0.834711   0.2       34.0345
            0.917355   0.2       30.4729
            1.0        0.13      26.4815]
    r = propgeom[:, 1] * Rtip
    c_rotor = propgeom[:, 2] * Rtip
    twist_rotor = propgeom[:, 3] * pi/180
    # af = CC.AlphaAF("examples/aurora/vlm_bem/naca4412.dat")
    af = (x)->false
    pos = [5.07 -18.75 6.73
            19.2 -18.75 9.01
            4.63 -8.13 7.04
            18.76 -8.45 9.3
            4.63 8.45 7.04
            18.76 8.45 9.3
            5.07 18.75 6.73
            19.2 18.75 9.01]*mpft


    Ry = (theta)->[cosd(theta) 0.0 sind(theta)
                    0.0 1.0 0.0 
                    -sind(theta) 0.0 cosd(theta)]
    axis = Ry(90)*[-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]
    dir = [1.0, 0.0, 0.0]

    CW = [true, true, true, true, false, false, false, false]
    pusher_propgeom = [
        0.0909091  0.08      46.75
        0.173554   0.110336  47.0562
        0.256198   0.137292  46.9327
        0.338843   0.16024   46.3793
        0.421487   0.178556  45.3963
        0.504132   0.191612  43.9834
        0.586777   0.198783  42.1409
        0.669422   0.2       39.8685
        0.752066   0.2       37.1664
        0.834711   0.2       34.0345
        0.917355   0.2       30.4729
        1.0        0.13      26.4815]
        
    Rtip_p = 4.5*mpft
    Rhub_p = 0.09*Rtip_p
    n_p = 6
    r_p = pusher_propgeom[:, 1] * Rtip_p
    c_p = pusher_propgeom[:, 2] * Rtip_p
    twist_p = pusher_propgeom[:, 3] * pi/180
    axis_p = [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]
    dir_p = [1.0, 0.0, 0.0]
    pos_p = [31.94, 0.0, 7.79]*mpft
    CW_p = true

    CW = vcat(CW, CW_p)
    pos = vcat(pos, pos_p')
    axes = zeros(3,3,9)
    for i in 1:8
        axes[:,:,i] = axis
    end 
    axes[:,:,end] = axis_p

    system = uns.vlm.WingSystem(; TF_trajectory=T, )


    b_wing = yle_wing[2]*2
    AR_wing = b_wing/c_wing[1]
    ns_wing = 24
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

    r_rotors = pos'
    dir_rotors = zeros(3,9)
    for i in 1:8
        dir_rotors[:,i] = [0.0, 0.0, 1.0]
    end
    dir_rotors[:,end] = [1.0, 0.0, 0.0]

    rotor_files = vcat(fill("apc10x7.csv", 8), "apc10x7_5blade.csv")
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

    rotor_systems = (rotors[[1,3,5,7]], rotors[[2,4,6,8]], rotors[[9]])
    

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