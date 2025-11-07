using FLOWUnsteady 
uns = FLOWUnsteady

function generate_joby_vehicle(T)
    mpft = 0.3048/2# general parameters
    m = 2400.0  #confirm
    I = 8694.0  #confirm
    cog = [0.8075, 0.0, 0.0]
    rho = 1.225
    mu = 1.81e-5
    g = 9.81
    # wing
    xle_wing = [0.0, -0.063, -0.2114]
    yle_wing = [0.0, 2.317, 5.486]
    zle_wing = [0.0, 0.289, 0.114]
    c_wing = [1.615, 1.204, 0.7925]
    twist_wing = [0.0, 0.0, 0.0]*pi/180 .+ 5*pi/180
    # airfoil = DF.readdlm("files/airfoils/LS(1)-0417.dat")        # check units
    camber_wing = fill((xc) -> 0, length(xle_wing))
    mirror_wing = true
    # xle_tail = [0.0, -0.74] .+ 3.17
    xle_tail = [0.0, 0.0] .+ 3.17
    yle_tail = [0.0, 2.618]
    zle_tail = [1.81, 1.81]/2
    c_tail = [1.006, 1.006]*0.75
    twist_tail = [0.0, 0.0]*pi/180
    camber_tail = fill((xc) -> 0, length(xle_tail))
    mirror_tail = true

    xle_cs = xle_tail + c_tail
    yle_cs = yle_tail
    zle_cs = zle_tail
    c_cs = c_tail*0.25
    twist_cs = [0.0, 0.0]*pi/180

    Rtip = 1.524

    Ry = (theta)->[cosd(theta) 0.0 sind(theta)
                    0.0 1.0 0.0 
                    -sind(theta) 0.0 cosd(theta)]
    axis = [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]

    CW = [true, false, true, false, true, false]
        
    axes = zeros(3,3,6)
    for i in 1:6
        axes[:,:,i] = axis
    end 

    system = uns.vlm.WingSystem(; TF_trajectory=T)


    b_wing = yle_wing[end]*2
    AR_wing = b_wing/c_wing[1]
    ns_wing = 12
    r_wing = [0.0, yle_wing[2]/yle_wing[end], 1.0]
    clen_wing = c_wing/c_wing[1]
    sweep_wing = [atand(xle_wing[2]-xle_wing[1], yle_wing[2]-yle_wing[1]), atand(xle_wing[3]-xle_wing[2], yle_wing[3]-yle_wing[2])]
    dihed_wing = [atand(zle_wing[2]-zle_wing[1], yle_wing[2]-yle_wing[1]), atand(zle_wing[3]-zle_wing[2], yle_wing[3]-yle_wing[2])]
    wing = uns.vlm.complexWing(T(b_wing), T(AR_wing), ns_wing, T.(r_wing), T.(clen_wing), T.(twist_wing*180/pi), T.(sweep_wing), T.(dihed_wing), T)
    uns.vlm.setcoordsystem(wing, [xle_wing[1], yle_wing[1], zle_wing[1]], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    uns.vlm.addwing(system, "wing", wing)

    b_tail = yle_tail[end]*2
    AR_tail = b_tail/c_tail[1]
    ns_tail = 6
    r_tail = [0.0, 1.0]
    clen_tail = c_tail/c_tail[1]
    sweep_tail = [atand(xle_tail[2]-xle_tail[1], yle_tail[2]-yle_tail[1])]
    dihed_tail =[atand(zle_tail[2]-zle_tail[1], yle_tail[2]-yle_tail[1])]
    tail = uns.vlm.complexWing(T(b_tail), T(AR_tail), ns_tail, T.(r_tail), T.(clen_tail), T.(twist_tail*180/pi), T.(sweep_tail), T.(dihed_tail), T)
    uns.vlm.setcoordsystem(tail, [xle_tail[1], yle_tail[1], zle_tail[1]], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    uns.vlm.addwing(system, "tail", tail)

    b_cs = yle_cs[end]*2
    AR_cs = b_cs/c_cs[1]
    ns_cs = 6
    r_cs = [0.0, 1.0]
    clen_cs = c_cs/c_cs[1]
    sweep_cs = [atand(xle_cs[2]-xle_cs[1], yle_cs[2]-yle_cs[1])]
    dihed_cs =[atand(zle_cs[2]-zle_cs[1], yle_cs[2]-yle_cs[1])]
    cs = uns.vlm.complexWing(T(b_cs), T(AR_cs), ns_cs, T.(r_cs), T.(clen_cs), T.(twist_cs*180/pi), T.(sweep_cs), T.(dihed_cs), T)
    uns.vlm.setcoordsystem(cs, [xle_cs[1], yle_cs[1], zle_cs[1]], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    uns.vlm.addwing(system, "cs", cs)

    vlm_system = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(vlm_system, "wing", wing)
    uns.vlm.addwing(vlm_system, "tail", tail)
    uns.vlm.addwing(vlm_system, "cs", cs)


    omega_ref = 100.0
    RPM_ref = omega_ref*30/pi
    n_ref = omega_ref/2/pi
    V_ref = 50.0
    D_ref = 2*Rtip
    J_ref = V_ref/n_ref/D_ref
    mu = 1.81e-5
    rho = 1.225

    r_surfs = [xle_wing[1] xle_tail[1] xle_cs[1]
               yle_wing[1] yle_tail[1] yle_cs[1]
               zle_wing[1] zle_tail[1] zle_cs[1]]


    r_rotors = [-1.13 yle_wing[2] zle_wing[2]
                -1.13 -yle_wing[2] zle_wing[2]
               xle_wing[3] yle_wing[3] zle_wing[3]
               xle_wing[3] -yle_wing[3] zle_wing[3]
               xle_tail[2] yle_tail[2] zle_tail[2]
               xle_tail[2] -yle_tail[2] zle_tail[2]]
    dir_rotors = zeros(3,6)
    for i in 1:6
        dir_rotors[:,i] = [0.0, 0.0, 1.0]
    end

    rotor_files = fill("apc10x7_mini.csv", 6)
    data_path = "src/rotor/"
    rotors = Vector{uns.vlm.Rotor}(undef, length(rotor_files))
    for i in eachindex(rotors)
        @show i
        rotors[i] = uns.generate_rotor(rotor_files[i]; TF_trajectory=T, pitch = 0.0, altReD=[RPM_ref, J_ref, mu/rho],
                                        n=7, blade_r=1/20, CW=CW[i], 
                                        xfoil=false, verbose=true, data_path=data_path, plot_disc=false)
        uns.vlm.setcoordsystem(rotors[i], r_rotors[i,:], axes[:,:,i])
        uns.vlm.addwing(system, "rotor$i", rotors[i])
    end

    rotor_systems = (rotors[1:2], rotors[3:4], rotors[5:6])
    

    wake_system = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(wake_system, "wing", wing)
    uns.vlm.addwing(wake_system, "tail", tail)
    uns.vlm.addwing(wake_system, "cs", cs)

    for i in eachindex(rotors)
        uns.vlm.addwing(wake_system, "rotor$i", rotors[i])
    end

    front_tilt = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(front_tilt, "front_right", rotors[1])
    uns.vlm.addwing(front_tilt, "front_left", rotors[2])

    middle_tilt = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(middle_tilt, "middle_right", rotors[3])
    uns.vlm.addwing(middle_tilt, "middle_left", rotors[4])

    back_tilt = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(back_tilt, "back_right", rotors[5])
    uns.vlm.addwing(back_tilt, "back_left", rotors[6])

    cs_tilt = uns.vlm.WingSystem(; TF_trajectory=T)
    uns.vlm.addwing(cs_tilt, "cs", cs)

    O1 = [front_tilt.wings[1]._wingsystem.O[1], 0.0, front_tilt.wings[1]._wingsystem.O[3]]
    for rotor in front_tilt.wings
        uns.vlm.setcoordsystem(rotor, [-0.48, rotor._wingsystem.O[2], 0.0], rotor._wingsystem.Oaxis)
    end
    uns.vlm.setcoordsystem(front_tilt, O1, front_tilt.Oaxis)

    O1 = [middle_tilt.wings[1]._wingsystem.O[1], 0.0, middle_tilt.wings[1]._wingsystem.O[3]]
    for rotor in middle_tilt.wings
        uns.vlm.setcoordsystem(rotor, [-0.48, rotor._wingsystem.O[2], 0.0], rotor._wingsystem.Oaxis)
    end
    uns.vlm.setcoordsystem(middle_tilt, O1, middle_tilt.Oaxis)

    O1 = [back_tilt.wings[1]._wingsystem.O[1], 0.0, back_tilt.wings[1]._wingsystem.O[3]]
    for rotor in back_tilt.wings
        uns.vlm.setcoordsystem(rotor, [-0.48, rotor._wingsystem.O[2], 0.0], rotor._wingsystem.Oaxis)
    end
    uns.vlm.setcoordsystem(back_tilt, O1, back_tilt.Oaxis)

    O1 = [cs_tilt.wings[1].O[1], 0.0, cs_tilt.wings[1].O[3]]
    uns.vlm.setcoordsystem(cs_tilt.wings[1], zeros(3), cs_tilt.wings[1].Oaxis)
    uns.vlm.setcoordsystem(cs_tilt, O1, cs_tilt.Oaxis)

    tilting_systems = (front_tilt, middle_tilt, back_tilt, cs_tilt)

    # build vehicle
    vehicle = uns.VLMVehicle(system;
            V = zeros(T, 3),
            W = zeros(T, 3),
            tilting_systems=tilting_systems,
            rotor_systems=rotor_systems,
            vlm_system=vlm_system,
            wake_system=wake_system)

    return vehicle, r_rotors, dir_rotors, r_surfs
end
