module ODE
export FundamentalNumericalBase, ExtendedNumericalBase, euler_step, euler_method, euler_method_with_local_errors, heuns_step, heuns_method
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
        println(i," ", iterations)
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

function make_method(step::Function)
    function ODE_method(
            nb::NumericalBase
            , x0::Float64
            , time_start::Float64
            , time_stop::Float64
        )

        iterations = Int(ceil((time_stop-time_start)/nb.Δt))

        X = zeros(iterations) 
        X[1] = x0
        t = time_start

        for i in 2:iterations
            X[i] = step(nb,X[i-1],t)
            t+= nb.Δt
        end

        return X
    end
end


function make_method_with_local_errors(step::Function)
    function ODE_method_with_errors(
        nb::ExtendedNumericalBase
        , x0::Float64
        , time_start::Float64
        , time_stop::Float64
    )

        iterations = Int(ceil((time_stop-time_start)/nb.Δt))

        X = zeros(iterations) #Incrementally corrected sequence
        Y = zeros(iterations) #Actual Values
        Z = zeros(iterations) #Full estimated sequence

        X[1] = x0
        Y[1] = x0
        Z[1] = x0

        t = time_start

        for i in 2:iterations
            X[i] = step(nb,Y[i-1],t)
            Y[i] = nb.f(t)
            Z[i] = step(nb,Z[i-1],t)
            t+= nb.Δt
        end

        return X,Y,Z
    end
    return ODE_method_with_errors
end


function heuns_step(fb::NumericalBase,y0,t)
    return y0 + 0.5 * ( nb.a(y0,t) + nb.a(euler_step(nb,y0,t),t)) * nb.Δt
end

heuns_method = make_method(heuns_step)
heuns_method_with_local_errors = make_method_with_local_errors(heuns_step)

end #end module