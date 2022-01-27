module ODE
export NumericalBase,FundamentalNumericalBase, ExtendedNumericalBase ,Significance, MachineLevel, SignificantDigits, SignificantFigures,euler_step, euler_method, euler_method_with_local_errors, heuns_step, heuns_method, heuns_method_with_local_errors,ab_method_with_local_errors,richardson_extrapolation_of_euler

#= Significant digits/figures type
This simplifies keeping track of rounding
=#
abstract type Significance end

struct MachineLevel end

struct SignificantDigits
    level::Int
end

struct SignificantFigures
    level::Int
end


function round(x::Float64, d::SignificantFigures)
    return round(x, sigdigits=d.level)
end
function round(x::Float64, d::SignificantDigits)
    return round(x, digits=d.level)
end
function round(x::Float64, d::MachineLevel)
    return x
end


#=
Way to keep track of numerical simulations to complete
=#
abstract type NumericalBase end
struct FundamentalNumericalBase <: NumericalBase
    a::Function
    Δt::Float64
    sig::Significance #used for estimation precision
end

struct ExtendedNumericalBase <: NumericalBase
    f::Function
    a::Function
    Δt::Float64
    sig::Significance #used for estimation precision
end
#Add some other construction methods.

#=
Euler Methods
=#
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
        X[i] = round(euler_step(nb,X[i-1],t), nb.sig)
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
        X[i] = round(euler_step(nb,Y[i-1],t), nb.sig)
        Y[i] = nb.f(t)
        Z[i] = round(euler_step(nb,Z[i-1],t), nb.sig)
        t+= nb.Δt
    end

    return X,Y,Z
end


#=
Heun's or the trapezoidal Method
=#

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
        X[i] = round(heuns_step(nb,Y[i-1],t), nb.sig)
        Y[i] = nb.f(t)
        Z[i] = round(heuns_step(nb,Z[i-1],t), nb.sig)
        t+= nb.Δt
    end

    return X,Y,Z
end


#=
3 step adam brashford
=#
function adam_brashford_step(
        nb::NumericalBase
        ,y2::Float64 #value at step 2
        ,y1::Float64 #value at step 1
        ,y0::Float64 #value at step 0
        ,t0::Float64 #initial time
    )
    t1 = t0 + nb.Δt
    t2 = t1 + nb.Δt

    return y2 + ( 23*nb.a(y2,t2) - 16*nb.a(y1,t1) + 5*nb.a(y0,t0))/12
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

    #Prep the full estimation
    Z[2] = round(heuns_step(nb,Z[1],t1), nb.sig)
    Z[3] = round(heuns_step(nb,Z[2],t2), nb.sig)

    #prep the actual values
    Y[2] = nb.f(t1)
    Y[3] = nb.f(t2)

    #Find the first few steps with local errors
    X[2] = round(heuns_step(nb,Y[1],t1), nb.sig)
    X[3] = round(heuns_step(nb,Y[2],t2), nb.sig)


    t = 0.0
    for i in 4:iterations
        X[i] = round(adam_brashford_step(nb,Y[i-1],Y[i-2],Y[i-3],time_start), nb.sig)
        Y[i] = nb.f(t)
        Z[i] = round(adam_brashford_step(nb,Z[i-1],Z[i-2],Z[i-3],time_start), nb.sig)
        t+= nb.Δt
    end

    return X,Y,Z
end

#=
Richardson Extrapolation of euler method
=#

function richardson_extrapolation_of_euler(nb::NumericalBase, z0, start_time, end_time)
    nb2 = FundamentalNumericalBase(nb.a,nb.Δt/2, nb.sig)

    Yn = euler_method(nb,z0, start_time, end_time)
    Y2n = euler_method(nb2,z0, start_time, end_time)

    return [2*Y2n[2*n] - Yn[n] for n=1:length(Yn)]
end

function richardson_extrapolation_of_euler_with_local_error(nb::ExtendedNumericalBase, z0, start_time, end_time)
    nb2 = ExtendedNumericalBase(nb.f, nb.a, nb.Δt/2, nb.sig)

    Xn,Yn,Zn = euler_method_with_local_errors(nb,z0, start_time, end_time)
    X2n,Y2n,Z2n = euler_method_with_local_error(nb2,z0, start_time, end_time)

    X_r =  [2*X2n[2*n] - Xn[n] for n=1:length(Yn)]
    Y_r =  [2*Y2n[2*n] - Yn[n] for n=1:length(Yn)]
    Z_r =  [2*Z2n[2*n] - Zn[n] for n=1:length(Yn)]

    #how do I pull local errors out of this?
end

end #end module
