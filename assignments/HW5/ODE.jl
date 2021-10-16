module ODE
export NumericalBase, euler_step, euler_method,heuns_method
#implement Euler
struct NumericalBase
    a::Function
    Δt::Float64
end

function euler_step(fb::NumericalBase,x0,t)
    return x0 + fb.Δt*fb.a(x0,t)
end
function euler_method(nb::NumericalBase, x0::Float64, time_start::Float64, time_stop::Float64)
    
    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
    
    X = zeros(iterations) 
    X[1] = x0
    t = time_start
    
    for i in 2:iterations
        X[i] = euler_step(nb,X[i-1],t)
        t+= nb.Δt
    end
    
    return X
end

function heuns_step(fb::NumericalBase,y0,t)
    return y0 + 0.5 * ( nb.a(y0,t) + nb.a(euler_step(nb,y0,t),t)) * nb.Δt
end


end #end module