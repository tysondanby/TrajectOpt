abstract type Model end
abstract type Parameters end

struct Conventional{T1,T2} <: Parameters
    m::T1 #Float64
    I::T1
    X::T1
    L::T1  #tail moment arm
    SWing::T1
    bWing::T1
    STail::T1
    bTail::T1
    rho::T1
    mu::T1
    g::T1
    unames::T2#Matrix{String}

    function Conventional(m,I,X,L,Swing,bwing,Stail,bTail,rho,mu,g)
        unames= ["Thrust" "Elevator"]
        new{typeof(m),typeof(unames)}(m,I,X,L,Swing,bwing,Stail,bTail,rho,mu,g,unames)
    end
end

struct BiWingTailSitter{T1,T2} <: Parameters
    m::T1
    I::T1
    X::T1
    s::T1  #distance between wings
    S::T1
    b::T1#
    rho::T1 #TODO: maybe these should be defined somewhere else. they are properties of the atmosphere, not the aircraft.
    mu::T1
    g::T1
    unames::T2#Matrix{String}

    function BiWingTailSitter(m,I,X,s,S,b,rho,mu,g)
        unames= ["Top Thrust" "Bottom Thrust"]
        new{typeof(m),typeof(unames)}(m,I,X,s,S,b,rho,mu,g,unames)
    end

end

struct LowFidel{Parameters} <: Model
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
        cTail = STail/bWing
        L = parameters.L
        X = parameters.X
        m = parameters.m
        rho = parameters.rho
        mu = parameters.mu
        #current wing aerodynamics
        alphaWing = theta - gamma
        ReWing = rho*vinf*cWing/mu
        CdWing, ClWing, CmWing = wing_polar_function(alphaWing, ReWing)
        #current tail aerodynamics
        alphaTail = alphaWing + deflection
        ReTail = ReWing*cTail/cWing
        CdTail, ClTail, CmTail = tail_polar_function(alphaTail, ReTail)
        #forces on each surface
        f = [[CdWing, CdTail] [ClWing, ClTail] [CmWing, CmTail]]'
        f[1:2,1] *= 1/2*rho*vinf^2*SWing
        f[1:2,2] *= 1/2*rho*vinf^2*STail
        f[3,1] *= 1/2*rho*vinf^2*SWing*cWing
        f[3,2] *= 1/2*rho*vinf^2*STail*cTail
        #add thrust
        # f[1,:] .-= thrust*cosd(alpha)
        # f[2,:] .+= thrust*sind(alpha)
        #transform forces into total force and moment
        # f[1,:] .*= -1
        F = [u[1]*cosd(alphaWing) - sum(f[1,:]) - m*9.81*sind(gamma), u[1]*sind(alphaWing) + sum(f[2,:]) - m*9.81*cosd(gamma)]
        M = sum(f[3,:]) - X*(f[1,1]*sind(alphaWing) + f[2,1]*cosd(alphaWing)) - (X+L)*(f[1,2]*sind(alphaWing) + f[2,2]*cosd(alphaWing))
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
        # f[1,:] .-= 2*u*cosd(alpha)
        # f[2,:] .+= 2*u*sind(alpha)
        #forces in the TANGENTIAL and NORMAL directions with respect to the inertial frame
        # f[1,:] .*= -1
        F = [sum(u)*cosd(alpha) - sum(f[1,:]) - m*9.81*sind(gamma), sum(u)*sind(alpha) + sum(f[2,:]) - m*9.81*cosd(gamma)]
        #total moment due to forces and torques
        M = sum(f[3,:]) + s/2*(u[2] - u[1] + (f[1,1] - f[1,2])*cosd(alpha) + (f[2,2] - f[2,1])*sind(alpha)) -
            X*((f[2,1] + f[2,2])*cosd(alpha) + (f[1,1] + f[1,2])*sind(alpha))
        return F, M
    end
    return forces_CRC3
end
