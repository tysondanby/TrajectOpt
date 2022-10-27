
function optimize_trim(x0, u0, model, final)
    trim_objective = trim_objective_constructor(final, model)
    designVars = zeros(length(u0) + length(x0))
    designVars[1:2] = u0
    designVars[3:end] = x0
    ng = 11
    lx = -Inf*ones(length(designVars))  # lower bounds on x
    ux = Inf*ones(length(designVars))  # upper bounds on x
    lg = zeros(ng)
    ug = [0,0,0,0,Inf,Inf,Inf,0,0,0,0]
    ip_options = Dict(
        "tol" => 1e-6,
        "max_iter" => 3000
        )
    solver = IPOPT(ip_options)
    options = Options(;solver)
    xopt, fopt, info = minimize(trim_objective, designVars, ng, lx, ux, lg, ug,options)
    return xopt, fopt
end

function trim_objective_constructor(final, model)#Use optimizer to find a trimmed state
    function trim_objective(g, designVars)#
        u = designVars[1:2]
        uSpline = [FM.Akima(0:2,u[1]*ones(3)), FM.Akima(0:2,u[2]*ones(3))]
        x = designVars[3:end]
        obj = u[1]#sum(u .^2)
        dx = dynamics_2D!(x, x, (uSpline, model), 1) #look up best practices
        g[1] = dx[1]
        g[2] = dx[2]
        g[3] = dx[3]
        g[4] = dx[4]
        g[5] = u[1] #TODO: why is this input constrained to zero? is it torque?
        # g[6] = u[2]
        g[7] = x[2]
        g[8] = x[6] - final[6]
        return obj
    end
    return trim_objective
end

function optimize_trajectory(initial, final, us, tFinal, model)
    designVars = vcat(us[1,:], us[2,:], tFinal)
    trajectory_objective = trajectory_objective_constructor(initial, final, model)
    ng = 3
    lu = zeros(length(designVars))
    uu = 5*ones(length(designVars))
    lu[Int((end-1)/2):end-1] .= -30
    uu[Int((end-1)/2):end-1] .= 30
    uu[end] = Inf
    lg = zeros(ng)
    ug = zeros(ng)
    ip_options = Dict(
        "tol" => 1e-1,
        "max_iter" => 1000
        )
    solver = SNOPT()
    options = Options(derivatives = ForwardAD(); solver)
    xopt, fopt, info = minimize(trajectory_objective, designVars, ng, lu, uu, lg, ug, options)
    return xopt, fopt, info
end

#Minimize Thrust squared with end trim state constrained
function trajectory_objective_constructor(initial, final, model)#Give an inital state, desired final state, and a model. constructs an objective function with constraints.
    function trajectory_objective(g,designVars)#g are constraints (residual functions, we want all g values to be zero) and designVars are the things we can change (thrust spline points and final time)
        #designVars[:1:end-1] = thrust spline points, designVars[end] = time
        # obj = sum(designVars[1:end-1].^2)
        obj = sum(designVars[1:Int((end-1)/2)].^2)   #Thrust squared
        #rearrange designVars
        us = transpose(reshape(designVars[1:end-1],Int((length(designVars)-1)/2),2))
        #get values from the dual numbers
        tSpan = [0,designVars[end]]
        #create splines
        t = range(0, stop = tSpan[2], length = length(us[1,:]))
        uSpline = [FM.Akima(t,us[1,:]), FM.Akima(t,us[2,:])]
        x = simulate(initial, uSpline, model, tSpan)
        # dx = dynamics!(x[:,end],x[:,end],(model, uSpline),t)
        #constrain final states to the desired final states
        g[1] = x[6,end] - final[6]
        g[2] = x[2,end] - final[2]
        g[3] = x[1,end] - final[1]

        # mappedXIndices = Int.(round.(range(1, stop = length(x[1,:]), length = length(t))))
        # g[2:length(t)+1] = x[1,mappedXIndices] .- initial[1]
        # g[2:length(t)+1] = x[2,mappedXIndices]
        # g[length(t)+1:end] = x[4,mappedXIndices]

        # Vinf = zeros(length(x.u))
        # gamma = zeros(length(x.u))
        # theta_dot = zeros(length(x.u))
        # theta = zeros(length(x.u))
        # posx = zeros(length(x.u))
        # posy = zeros(length(x.u))
        # for i in 1:length(x.u)
        #     Vinf[i] = x.u[i][1]
        #     gamma[i] = x.u[i][2]
        #     theta_dot[i] = x.u[i][3]
        #     theta[i] = x.u[i][4]
        #     posx[i] = x.u[i][5]
        #     posy[i] = x.u[i][6]
        # end
        # Vinf_spline = FM.Akima(x.t,Vinf)
        # gamma_spline = FM.Akima(x.t,gamma)
        # theta_spline = FM.Akima(x.t,theta)
        # posy_spline = FM.Akima(x.t,posy)
        # Vinf_map = Vinf_spline.(t)
        # gamma_map = gamma_spline.(t)
        # theta_map = theta_spline.(t)
        # posy_map = posy_spline.(t)
        # open("optimization_outputs.txt", "a") do io
        #     DF.writedlm(io, [t, us[1,:], us[2,:], Vinf_map, gamma_map, theta_map, posy_map])
        # end
        return obj
    end
end

function optimize_trajectory_by_segments(initial, final, us, tFinal, nSegs, model)
    #create segements
    segments = range(initial, stop = final, length = nSegs)
    # @show segments
    splineLength = Int((length(designVars)-1)/2)
    splineSegment = 1:splineLength

    #first segment optimization
    uopt, _, _ = optimize_trajectory(initial, segments[2], us, tFinal, model)
    u = zeros(2,nSegs*length(splineSegment))
    u[1,splineSegment] = uopt[1:Int((end-1)/2)]
    u[2,splineSegment] = uopt[Int((end-1)/2)+1:end-1]
    t = range(0,stop = uopt[end],length = length(u[1,splineSegment]))
    tSpan = zeros(2,length(segments))
    tSpan[:,1] = [0 t[end]]
    uSpline = [Akima(t,u[1,splineSegment]),Akima(t,u[2,splineSegment])]
    segmentPath = simulate(initial, uSpline, model, tSpan[:,1])

    for i in 2:nSegs-1
        #determine next segment
        state1 = segmentPath.u[end]
        state2 = segments[i+1]
        uNew = [u[1,splineSegment[end]],u[2,splineSegment[end]]]    #maybe trim here

        #set initial design variables
        designVars[1:Int((end-1)/2)] .= uNew[1]
        designVars[Int((end-1)/2)+1:end-1] .= uNew[2]
        designVars[end] = t[end]

        #optimize
        uopt = optimize_trajectory(state1, state2, uNew, t[end], model)
        splineSegment = splineLength*(i-1)+1:splineLength*i
        u[1,splineSegment] = uopt[1][1:Int((end-1)/2)]
        u[2,splineSegment] = uopt[1][Int((end-1)/2)+1:end-1]
        t = range(0,stop = uopt[1][end],length = length(u[1,splineSegment]))
        tSpan[:,i] = [0 t[end]]
        uSpline = [Akima(t,u[1,splineSegment]),Akima(t,u[2,splineSegment])]
        segmentPath = simulate(state1, uSpline, model, tSpan[:,i])
    end

    uopt, _, _ = optimize_trajectory(segmentPath.u[end], segments[end], uNew, t[end], model)
    u[1,splineSegment[end]+1:end] = uopt[1:Int((end-1)/2)]
    u[2,splineSegment[end]+1:end] = uopt[Int((end-1)/2)+1:end-1]
    t = range(0,stop = uopt[end],length = length(u[1,splineSegment[end]+1:end]))
    tSpan[:,end] = [0 t[end]]
    uSpline = [Akima(t,u[1,splineSegment[end]+1:end]),Akima(t,u[2,splineSegment[end]+1:end])]

    tFinal = sum(tSpan[2,:])
    t = range(0,stop = tFinal, length = length(u[1,:]))
    fopt = sum(u.^2)
    uoptSpline = [Akima(t,u[1,:]),Akima(t,u[2,:])]
    return t, uoptSpline, fopt, u
end

function optimization_tracker()

    Vinf = zeros(length(x.u))
    gamma = zeros(length(x.u))
    theta_dot = zeros(length(x.u))
    theta = zeros(length(x.u))
    posx = zeros(length(x.u))
    posy = zeros(length(x.u))
    for i in 1:length(x.u)
        Vinf[i] = x.u[i][1]
        gamma[i] = x.u[i][2]
        theta_dot[i] = x.u[i][3]
        theta[i] = x.u[i][4]
        posx[i] = x.u[i][5]
        posy[i] = x.u[i][6]
    end
    Vinf_spline = FM.Akima(x.t,Vinf)
    gamma_spline = FM.Akima(x.t,gamma)
    theta_spline = FM.Akima(x.t,theta)
    posy_spline = FM.Akima(x.t,posy)
    Vinf_map = Vinf_spline.(t)
    gamma_map = gamma_spline.(t)
    theta_map = theta_spline.(t)
    posy_map = posy_spline.(t)
    open("optimization_outputs.txt", "a") do io
        DF.writedlm(io, [t, us[1,:], us[2,:], Vinf_map, gamma_map, theta_map, posy_map])
    end
end



# function optimize_trajectory3(initial,final,designVars,model)
#     trajectory_objective = trajectory_objective_constructor3(initial,final,model)
#     ng = 1
#     lu = zeros(length(designVars))  # lower bounds on x
#     uu = 10*ones(length(designVars))  # upper bounds on x
#     uu[end] = Inf
#     lg = zeros(ng)
#     ug = zeros(ng)
#     ip_options = Dict(
#         "tol" => 1e-3,
#         "max_iter" => 1000
#         )
#     solver = IPOPT(ip_options)
#     options = Options(derivatives = ForwardAD();solver)
#     xopt, fopt, info = minimize(trajectory_objective, designVars, ng, lu, uu, lg, ug, options)
#     return xopt, fopt, info
# end

# function trajectory_objective_constructor3(initial,final,model)
#     function trajectory_objective(g,designVars)
#         obj = sum(designVars .^2)
#         s = model.s
#         uTop = designVars[1:end-1]
#         uBottom = zeros(length(uTop))
#         tspan = [0 designVars[end].value]
#         t = range(0,stop = tspan[2], length = length(uTop))
#         x = initial
#         for i in 1:length(t)
#             usCurrent = [uTop[i].value, 0]
#             if i < length(t)
#                 dt = t[i+1] - t[i]
#             else
#                 dt = t[i] - t[i-1]
#             end
#             F, M = forces(x,usCurrent,model)
#             if M < 0
#                 uBottom[i] = -M/s
#             end
#             usCurrent = [uTop[i].value, uBottom[i]]
#             F, M = forces(x,usCurrent,model)
#             x, dx = step!(x,usCurrent,dt,model)    #possibly use DE here too
#             @show i x dx usCurrent
#         end
#         us = zeros(2,length(uTop))
#         for i in 1:length(uTop)
#             us[1,i] = uTop[i].value
#         end
#         us[2,:] = uBottom
#         u = [FM.Akima(t,us[1,:]),FM.Akima(t,us[2,:])]
#         x = simulate(initial,model,u,tspan)
#         # @show x[6,:]
#         g[1] = x[6,end] - final[6]

#         Vinf = zeros(length(x.u))
#         gamma = zeros(length(x.u))
#         theta_dot = zeros(length(x.u))
#         theta = zeros(length(x.u))
#         posx = zeros(length(x.u))
#         posy = zeros(length(x.u))
#         for i in 1:length(x.u)
#             Vinf[i] = x.u[i][1]
#             gamma[i] = x.u[i][2]
#             theta_dot[i] = x.u[i][3]
#             theta[i] = x.u[i][4]
#             posx[i] = x.u[i][5]
#             posy[i] = x.u[i][6]
#         end
#         Vinf_spline = FM.Akima(x.t,Vinf)
#         gamma_spline = FM.Akima(x.t,gamma)
#         theta_spline = FM.Akima(x.t,theta)
#         posy_spline = FM.Akima(x.t,posy)
#         Vinf_map = Vinf_spline.(t)
#         gamma_map = gamma_spline.(t)
#         theta_map = theta_spline.(t)
#         posy_map = posy_spline.(t)
#         open("optimization_outputs.txt", "a") do io
#             DF.writedlm(io, [t, us[1,:], us[2,:], Vinf_map, gamma_map, theta_map, posy_map])
#         end
#         open("constraint_outputs.txt", "a") do io
#             DF.writedlm(io, [g[1].value])
#         end
#         return obj
#     end
# end

# function optimize_trajectory4(initial,final,designVars,nSegs,model)
#     for i in 1:nSegs

#     end
#     trajectory_objective = trajectory_objective_constructor2(initial,final,model)
#     ng = 3
#     lu = zeros(length(designVars))  # lower bounds on x
#     uu = 10*ones(length(designVars))  # upper bounds on x
#     uu[end] = Inf
#     lg = zeros(ng)
#     ug = zeros(ng)
#     ip_options = Dict(
#         "tol" => 1e-3,
#         "max_iter" => 1000
#         )
#     solver = IPOPT(ip_options)
#     options = Options(derivatives = ForwardAD();solver)
#     xopt, fopt, info = minimize(trajectory_objective, designVars, ng, lu, uu, lg, ug, options)
#     return xopt, fopt, info
# end





















# function optimize_bottom_thrust(x,u_top,model)
#     thrust_objective = thrust_objective_constructor(x,model)
#     ng = 1
#     lu = 0  # lower bounds on x
#     uu = Inf  # upper bounds on x
#     lg = zeros(ng)
#     ug = .001*ones(ng)
#     ip_options = Dict(
#         "tol" => 1e-3,
#         "max_iter" => 3
#         )
#     solver = IPOPT(ip_options)
#     options = Options(derivatives = CentralFD();solver)
#     xopt, fopt, info = minimize(thrust_objective, designVars, ng, lu, uu, lg, ug, options)
#     return xopt, fopt, info
# end

# function thrust_objective_constructor(x,u_top,model)
#     function thrust_objective(g,u_bottom)
#         obj = u_bottom^2
#         u_spline = [FM.Akima(0:2,u_top*ones(3)), FM.Akima(0:2,u_bottom*ones(3))]
#         dx = dynamics!(zeros(6),x,(model,u_spline),1)
#         g[1] = dx[1]
#         g[2] = dx[2]
#         g[3] = dx[4]
#         return obj
#     end
# end
