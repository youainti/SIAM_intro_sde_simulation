module ODE
export NumericalBase,FundamentalNumericalBase, ExtendedNumericalBase, euler_step, euler_method, euler_method_with_local_errors, heuns_step, heuns_method, heuns_method_with_local_errors,ab_method_with_local_errors,richardson_extrapolation_of_euler
#implement Euler
abstract type NumericalBase end

struct FundamentalNumericalBase <: NumericalBase
    a::Function
    Δt::Float64
end

struct ExtendedNumericalBase <: NumericalBase
    f::Function
    a::Function
    Δt::Float64
end
#Add some other construction methods.

function euler_step(fb::NumericalBase,x0,t)
    return x0 + fb.Δt*fb.a(x0,t)
end

function euler_method(nb::NumericalBase, x0::Float64, time_start::Float64, time_stop::Float64)
    
    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
    
    X = zeros(iterations) 
    X[1] = x0
    t = time_start
    
    for i in 2:iterations
        #println(i," ", iterations)
        X[i] = euler_step(nb,X[i-1],t)
        t+= nb.Δt
    end
    
    return X
end




function euler_method_with_local_errors(nb::ExtendedNumericalBase, x0::Float64, time_start::Float64, time_stop::Float64)
    
    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
    
    X = zeros(iterations) #Local Errors
    Y = zeros(iterations) #Actual Values
    Z = zeros(iterations) #Full estimation
    
    
    X[1] = x0
    Y[1] = x0
    Z[1] = x0
    
    t = time_start
    
    for i in 2:iterations
        X[i] = euler_step(nb,Y[i-1],t)
        Y[i] = nb.f(t)
        Z[i] = euler_step(nb,Z[i-1],t)
        t+= nb.Δt
    end
    
    return X,Y,Z
end




function heuns_step(nb::NumericalBase,y0,t)
    return y0 + 0.5 * ( nb.a(y0,t) + nb.a(euler_step(nb,y0,t),t)) * nb.Δt
end

function heuns_method_with_local_errors(nb::ExtendedNumericalBase, x0::Float64, time_start::Float64, time_stop::Float64)
    
    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
    
    X = zeros(iterations) #Local Errors
    Y = zeros(iterations) #Actual Values
    Z = zeros(iterations) #Full estimation
    
    
    X[1] = x0
    Y[1] = x0
    Z[1] = x0
    
    t = time_start
    
    for i in 2:iterations
        X[i] = heuns_step(nb,Y[i-1],t)
        Y[i] = nb.f(t)
        Z[i] = heuns_step(nb,Z[i-1],t)
        t+= nb.Δt
    end
    
    return X,Y,Z
end


function adam_brashford_step(
        nb::NumericalBase
        ,y2::Float64 #value at step 2
        ,y1::Float64 #value at step 1
        ,y0::Float64 #value at step 0
        ,t::Float64 #initial time
    )
    t1 = t - nb.Δt
    t2 = t1 - nb.Δt
    
    return y2 + ( 23*nb.a(y2,t) - 16*nb.a(y1,t1) + 5*nb.a(y0,t2))/12
end

function ab_method_with_local_errors(nb::ExtendedNumericalBase, x0::Float64, time_start::Float64, time_stop::Float64)
    
    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
    
    X = zeros(iterations) #Local Errors
    Y = zeros(iterations) #Actual Values
    Z = zeros(iterations) #Full estimation
    
    
    X[1] = x0
    Y[1] = x0
    Z[1] = x0
    
    t1 = time_start
    t2 = t1  + nb.Δt
    
    Z[2] = heuns_step(nb,Z[1],t1)
    Z[3] = heuns_step(nb,Z[2],t2)
    
    t = 0.0
    
    for i in 4:iterations
        X[i] = adam_brashford_step(nb,Y[i-1],Y[i-2],Y[i-3],time_start)
        Y[i] = nb.f(t)
        Z[i] = adam_brashford_step(nb,Z[i-1],Z[i-2],Z[i-3],time_start)
        t+= nb.Δt
    end
    
    return X,Y,Z
end

function richardson_extrapolation_of_euler(nb::NumericalBase, z0, start_time, end_time)
    nb2 = FundamentalNumericalBase(nb.a,nb.Δt/2)
    
    Yn = euler_method(nb,z0, start_time, end_time)
    Y2n = euler_method(nb2,z0, start_time, end_time)
    
    return 2*Y2n - Yn
end

end #end module