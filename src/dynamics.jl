function dynamics_2D!(dx, x, p, t)  
    #get states
    vinf, gamma, thetadot, theta, posx, posy = x
    #get parameters, inputs, and forces
    uSpline, model = p
    u = [uSpline[1](t),uSpline[2](t)]
    F, M = model.forces(x,u)
    #get physical parameters from model
    m = model.parameters.m 
    I = model.parameters.I
    #calculate time derivatives from equtions of motion
    k = F[2]/(m*vinf^2)
    dx[1] = F[1]/m  
    dx[2] = vinf*k  
    dx[3] = M/I
    dx[4] = thetadot
    dx[5] = vinf*cosd(gamma)
    dx[6] = vinf*sind(gamma)
    return dx
end

function simulate(x0, uSpline, model, tSpan)
    prob = DE.ODEProblem(dynamics_2D!, x0, tSpan, (uSpline, model))
    sol = solve(prob, abstol = 1e-3, reltol = 1e-3)     
    return sol
end

function plot_simulation(path, uSpline)
    Vinf_points = zeros(length(path.u))
    gamma_points = zeros(length(path.u))
    theta_dot_points = zeros(length(path.u))
    theta_points = zeros(length(path.u))
    posx_points = zeros(length(path.u))
    posy_points = zeros(length(path.u))
    for i in 1:length(path.u)
        Vinf_points[i] = path.u[i][1]
        gamma_points[i] = path.u[i][2]
        theta_dot_points[i] = path.u[i][3]
        theta_points[i] = path.u[i][4]
        posx_points[i] = path.u[i][5]
        posy_points[i] = path.u[i][6]
    end

    Vinf_spline = Akima(path.t, Vinf_points)
    gamma_spline = Akima(path.t, gamma_points)
    theta_spline = Akima(path.t, theta_points)
    posy_spline = Akima(path.t, posy_points)
    posx_spline = Akima(path.t, posx_points)

    t = range(path.t[1], stop = path.t[end], length = 200)
    Vinf = Vinf_spline.(t)
    u_top = uSpline[1].(t)
    u_bottom = uSpline[2].(t)
    gamma = gamma_spline.(t)
    theta = theta_spline.(t)
    posy = posy_spline.(t)
    posx = posx_spline.(t)
    aoa = theta - gamma

    p1 = plot(t, Vinf, xlabel = "Time (s)", ylabel = "State", label = "Velocity", legend = :topleft)
    plot!(t, gamma, label = "Flightpath Angle")
    plot!(t, theta, label = "Pitch Angle")
    plot!(t, aoa, label = "Angle of Attack")
    plot!(t, posy, label = "Y Position")
    # plot!(t, posx, label = "X Position")
    # p2 = plot(t, u_top, xlabel = "Time (s)", ylabel = "Thrust (N)", label = "Top", legend = :topleft)
    p2 = plot(t, u_top, xlabel = "Time (s)", ylabel = "Input", label = "Thrust", legend = :topleft)
    plot!(t, u_bottom, label = "Elevator")
    p = plot(p1,p2,layout = 2)
    display(p)
end

# #=============Euler integration scheme=============#

# """
# dynamics(x,F,M,model)

# """

# function dynamicsEU(x,F,M,model)
#     vinf, gamma, thetadot, theta, posx, posy = x
#     #get physical parameters from model
#     m = model.m 
#     I = model.I
#     #calculate time derivatives from equtions of motion
#     k = F[2]/(m*vinf^2)
#     dx = zeros(6)   # inefficient
#     dx[1] = F[1]/m  
#     dx[2] = vinf*k*180/pi
#     dx[3] = M/I*180/pi
#     dx[4] = thetadot
#     dx[5] = vinf*cosd(gamma) 
#     dx[6] = vinf*sind(gamma)
#     #@show dx
#     return dx
# end

# # take the current state and inputs and predict the next step
# function step!(x,u,dt,model)
#     F,M = forces(x,u,model)
#     dx = dynamicsEU(x,F,M,model)
#     x += dx*dt          
#     return x, dx
# end

# function simulateEU(x,u,t,model)
#     x_hist = zeros(length(x),length(t))
#     x_hist[:,1] = x
#     for i in 2:length(t)
#         dt = t[i] - t[i-1]
#         x, dx = step!(x,u[:,i-1],dt,model)
#         x_hist[:,i] = x
#     end
#     Vinf = x_hist[1,:]
#     gamma = x_hist[2,:]
#     theta_dot = x_hist[3,:]
#     theta = x_hist[4,:]
#     posx = x_hist[5,:]
#     posy = x_hist[6,:]
#     alpha = theta - gamma
#     return x_hist #Vinf, gamma, theta_dot, theta, posx, posy, alpha
# end

# function stability_check(x,u_spline,us,t,model)
#     stability = true
#     i = 2
#     while (stability == true)&&(i < length(t))
#         p = (model,u_spline)
#         dt = t[i] - t[i-1]
#         dx = dynamics!(zeros(6),x,p,t[i])
#         x += dx*dt
#         if (abs(x[2]) >= 50)||(abs(x[4]) >= 50)
#             us[:,i] = us[:,i-1]
#             u_spline = [FM.Akima(t,us[1,:]),FM.Akima(t,us[2,:])]
#         end
#         i += 1
#     end
#     return u_spline
# end
