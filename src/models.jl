abstract type Model end
abstract type Parameters end

struct Conventional <: Parameters
    m::Float64
    I::Float64
    X::Float64
    L::Float64  #tail moment arm
    SWing::Float64
    bWing::Float64
    STail::Float64
    bTail::Float64
    rho::Float64
    mu::Float64
end

struct BiWingTailSitter <: Parameters
    m::Float64
    I::Float64
    X::Float64
    s::Float64  #distance between wings
    S::Float64
    b::Float64
    rho::Float64
    mu::Float64
end

struct LowFidel <: Model
    parameters::Parameters
    forces::Function  
end

struct HighFidel <: Model
    # #VLM BEM parameters
    # m::Float64
    # I::Float64
    # X::Float64
    # s::Float64
end

function polar_constructor(Cds,Cls,Cms,alphas,Res)
    function polar_function(alpha, Re)         #interpolation
        Cd = interp2d(linear, alphas, Res, Cds, alpha, Re)
        Cl = interp2d(linear, alphas, Res, Cls, alpha, Re)
        Cm = interp2d(linear, alphas, Res, Cms, alpha, Re)
        return [Cd;Cl;Cm]
    end
    return polar_function
end

function conventional_forces_constructor(wing_polar_function, tail_polar_function, parameters::Conventional)
    function forces_conventional(x, u)
        #current state and inputs
        vinf, gamma, thetadot, theta, posx, posy = x
        thrust, deflection = u
        #get physical and environmental parameters
        SWing = parameters.SWing
        bWing = parameters.bWing
        STail = parameters.STail
        bTail = parameters.bTail
        cWing = SWing/bWing
        cTail - STail/bWing
        L = parameters.L
        X = parameters.X
        m = parameters.m
        rho = parameters.rho
        mu = parameters.mu
        #current wing aerodynamics
        alpha = theta - gamma
        Re = rho*vinf*cWing/mu
        CdWing, ClWing, CmWing = wing_polar_function(alpha, Re)
        #current tail aerodynamics
        alpha += deflection
        Re *= cTail/cWing
        CdTail, ClTail, CmTail = tail_polar_function(alpha, Re)
        #forces on each surface
        f = [[CdWing, CdTail] [ClWing, ClTail] [CmWing, CmTail]]
        f[1:2,1] *= 1/2*rho*vinf^2*SWing
        f[1:2,2] *= 1/2*rho*vinf^2*STail
        f[3,1] *= 1/2*rho*vinf^2*SWing*cWing
        f[3,2] *= 1/2*rho*vinf^2*STail*cTail
        #add thrust
        f[1,:] .-= thrust*cosd(alpha)
        f[2,:] .+= thrust*sind(alpha)
        #transform forces into total force and moment
        f[1,:] .*= -1
        F = [sum(f[1,:]) - m*9.81*sind(gamma), sum(f[2,:]) - m*9.81*cosd(gamma)]
        M = sum(f[3,:]) + X*((f[1,1] - f[1,1])*cosd(alpha) + (f[2,2] - f[2,1])*sind(alpha)) + 
            X*((f[2,1] + f[2,2])*cosd(alpha) - (f[1,1] + f[1,2])*sind(alpha))
        return F, M
    end
    return forces_conventional
end

function biwing_tailsitter_forces_constructor(polar_function, parameters::BiWingTailSitter)
    function forces_CRC3(x, u)
        vinf, gamma, thetadot, theta, posx, posy = x
        #get physical and environmental parameters
        S = parameters.S  
        b = parameters.b 
        c = S/b
        s = parameters.s
        X = parameters.X
        m = parameters.m
        rho = parameters.rho
        mu = parameters.mu
        #current top state
        alpha = theta - gamma
        v = vinf - s/2*thetadot*cosd(alpha)
        Re = rho*v*c/mu
        #calculate aerodynamic forces
        #top aerodynamics at current time step 
        CdTop, ClTop, CmTop = polar_function(alpha,Re)
        #current bottom state
        v = vinf + s/2*thetadot*cosd(alpha)
        Re = rho*v*c/mu
        #bottom aerodynamics current time step 
        CdBottom, ClBottom, CmBottom = polar_function(alpha,Re)
        f = [[CdTop, CdBottom] [ClTop, ClBottom] [CmTop, CmBottom]]'
        #denormalization
        f[1:2,:] *= 1/2*rho*v^2*S               # change v here!!
        f[3,:] *= 1/2*rho*v^2*S*c
        #add thrust 
        f[1,:] .-= 2*u*cosd(alpha)
        f[2,:] .+= 2*u*sind(alpha)
        #forces in the x and y directions with respect to the inertial frame
        f[1,:] .*= -1
        F = [sum(f[1,:]) - m*9.81*sind(gamma),sum(f[2,:])-m*9.81*cosd(gamma)]
        #total moment due to forces and torques
        M = sum(f[3,:]) + s/2*((f[1,2] - f[1,1])*cosd(alpha) + (f[2,2] - f[2,1])*sind(alpha)) + 
            X*((f[2,1] + f[2,2])*cosd(alpha) - (f[1,1] + f[1,2])*sind(alpha))
        return F, M
    end
    return forces_CRC3
end