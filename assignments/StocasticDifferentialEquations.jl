
#This is a module with useful functions for working with stochastic processes
module SDE
using LinearAlgebra, Statistics, Random

export WienerPath, StochasticPath, integrate, cumulative_integrate, refine_path, spread_path
#=
Information on importing data can be found here.
https://stackoverflow.com/questions/37200025/how-to-import-custom-module-in-julia

in effect just do
include("StocasticDifferentialEquations.jl")
using .StochasticDifferentialEquations
=#

#=
Create a struct to capture the type of exit, with a couple of helpful methods
=#
abstract type ExitCondition end
#note the implicit contract below
measured_time(exit_type::ExitCondition) = exit_type.time

struct Exit1dBox <: ExitCondition
    time::Int
    exited_top::Bool
end
exited_top_boundary(exit_type::Exit1dBox) = exit_type.exited_top
exited_bottom_boundary(exit_type::Exit1dBox) = !exit_type.exited_top



#= Wiener Paths
The following define a Wiener Process Structure and a set of useful functions 
for working with Weiner Paths
=#
struct WienerPath
    # the two key things to know are the sampled points from the path and 
    # how far apart they are in the coordinate system
    path::Vector{Float64}
    Δt::Float64
   
    #This overrides the default creation function.
    # It creates a wiener process using a given path, while 
    # enforcing the invariant that a WienerProcess starts at zero
    WienerPath(p::Vector{Float64},d::Float64) = p[1] != 0.0 ? new(vcat(0.0,p),d) : new(p,d)
end

# Generators
function WienerPath(size::Int, domain_length::Float64)
    #= 
    This function builds a wiener path out of randomly drawn steps and a domain length.
    It uses the domain length to calculate Δt.
    
    I am debating on turning this into a wiener_process function that will return M WienerPaths
    
    It returns a WienerPath
    =#
    return WienerPath(cumsum(randn(size)),domain_length/size)
end

# Other Methods
function get_steps(p::WienerPath)
    #=
    This pulls the steps out of the Wiener Path
    
    It depends on the invariates held up by the WienerPath struct
    
    It returns the steps as a Vector{Float64}
    =#
    len = length(p.path)
    steps = p.path[2:len] - p.path[1:len-1]
    
    return steps
end

function refine_path(p::WienerPath)
    #=
    This refines the temporal resolution of the given WienerPath.
    
    It does so according to the details on page 29 and 30 of the textbook mentioned in the README.
    
    It returns a new WienerPath with the new resolution.
    =#
    len = length(p.path)
    arr = zeros(2*len-1)
    
    for ind in 1:len-1
        arr[2*ind-1] = p.path[ind]
        arr[2*ind] =  (0.5 * (p.path[ind] + p.path[ind+1]) .+ (((√p.Δt)/2)*randn(1)))[1]
    end
    arr[2*len-1] = p.path[len]
    
    return WienerPath(arr, p.Δt/2)
end

function spread_path(p::WienerPath)
    #=
    This function spreads a given WienerPath to match a refinment without adding in the differences between them
    This allows you to easily compare the original path and the spread path when graphing, etc.
    
    It returns a new WienerPath at the new resolution, with intermetiate points made up of the average of surrounding points
    =#
    len = length(p.path)
    arr = zeros(2*len-1)
    
    for ind in 1:len-1
        arr[2*ind-1] = p.path[ind]
        arr[2*ind] =  (0.5 * (p.path[ind] + p.path[ind+1])) 
    end
    arr[2*len-1] = p.path[len]
    
    return WienerPath(arr, p.Δt/2)
end

function integrate(f::Function, wp::WienerPath )
    #=
        This function integrates across the time domain, by computing the riemann sum
        associated with the appropriate ito integral.

        Note that there is an alternate formulation using dot products, but I had already written this, 
        and my limited timing test suggested it was faster.

        it returns a Float64 value for the integral
    =#
        
    sum = 0.0
    time_tracker = 0.0
    
    steps = get_steps(wp)
    
    for i in 1:length(wp.path)-1
        #Calculate the running sum
        try
            sum += f(wp.path[i], time_tracker) * (steps[i]) 
        catch err
            if isa(err, MethodError)
                println("Method passed to integration has incorrect type signature")
                throw(err)
            else
                throw(err)
            end
        end
        #iterate time tracker
        time_tracker += wp.Δt
    end
    
    return sum
end


#=
Thoughts on writing a stochastic process module:
Start by working with just single differntial equation processes.
A StochProc struct should contain
 - Itself
    - a function and underlying wiener path + process.
 - What is needed to reconstruct it at higher resolutions.
    - which is just itself and a refine function that will recreate itself with a refined wiener process
Options are lazy and eagerly created:
 - lazy is there is a function to generate the stoch proc
    - recall, no member functions exist.
    - requires setting a new conventional function to call.
 - eager is that is that it's generation happens upon creation.
    - this is probably the way to go.
    - this will mean that in order to hang onto the stoch proc, you get the stoch proc and the underlying WP
After further thought I realized that these are, in a way, different concepts
A stochastic process contains the sets of functions etc that when composed and integrated, determine the final results.
A stochastic path is a single realization of the stochastic process.
So, a stochastic process should contain an ordered list of functions over which to integrate.
It should also have a method to generate stochastic paths.
A stochastic path should have the underlying functions and a wiener path that it is based on.
it should also 
=#
#=
stochastic process (lazy, it represents a generation procedure)
 It should contain
    bound and delta information
    the functions that make it up
Two important methods are 
   walk_until(exit_fn) -> ExitCondition(stochastic_path, exit_info): which creates the stochastic path by walking until an exit condition is met.
   walk_for(n) -> stochastic_path: which creates a stochastic path of length n
=#
struct StochasticProcess 
    Δt::Float64
    functions_list::Array{Function} #TODO: check later
end


function walk_until(sp::StochasticProcess, exit_fn::Function)
    return nothing #exit condition
end

function walk_for(sp::StochasticProcess,n::UInt)
    return nothing
end

#=
stochastic path
 it should contain
    wiener_path it is built on
    functions it is build out of
    the array of current values.
 important methods include
    integration *
    cumulative integration *
    refine_path *
    spread_path *
    time_period: the last period available 

There is also the question of whether or not to add a cumulative integration function.
it would need to add integration to it's functions list.

There is still a question of whether or not to build a stochastic path that doesn't have 
some of the fundamentals such as the base wiener path. 
This may be useful as a return type for some functions.
=#
struct StochasticPath 
    wiener_path::WienerPath
    functions_list::Vector{Function} #TODO: check later
    #These functions need to take in a path and then return a "path"
    # f(sp::StochasticPath) -> StochasticPath
    #See the cumulative_integration function for an example of what the new functions_list needs to look like.
    path::Vector{Float64}
    
    #change the creation method to verify the length of the WP and path are the same.
end

function path_to_process(sp::StochasticPath)
    return StochasticProcess(sp.wiener_path.Δt,sp.functions_list)
end

function cumulative_integrate(fn::Function, sp::StochasticPath)
    #=
        This function integrates across the time domain, by computing the riemann sum
        associated with the appropriate ito integral.
        
        the function must have the signature fn(X\_t, W\_t,\Delta t, t) -> Float
        Ideally fn = \mu dt + \sigma dW\_t
        
        it returns a Float64 value for the integral
    =#
        
    sum = 0.0
    time_tracker = 0.0
    cum_integral = zeros(length(sp.path)-1)
    
    steps = get_steps(sp.wiener_path)
    
    for i in 1:length(sp.wiener_path.path)-1
        #Calculate the running sum
        sum += fn(sp.path[i], steps[i], sp.wiener_path.Δt, time_tracker) 
        cum_integral[i] = sum
        #iterate time tracker
        time_tracker += sp.wiener_path.Δt
    end
    
    return StochasticPath(
        sp.wiener_path
        ,append!(vec(copy(sp.functions_list)), vec([x -> cumulative_integrate(fn,x)]))
        #Note^^: This adds a generic function to the array, i.e. unnamed when printed.
        ,cum_integral
    )
end
function integrate(fn::Function, sp::StochasticPath)
    #=
        This function integrates across the time domain, by computing the riemann sum
        associated with the appropriate ito integral.
        
        the function must have the signature fn(X\_t, W\_t,\Delta t, t) -> Float
        Ideally fn = \mu dt + \sigma dW\_t
        
        it returns a Float64 value for the integral
    =#
        
    sum = 0.0
    time_tracker = 0.0
    
    steps = get_steps(sp.wiener_path)
    
    for i in 1:length(sp.wiener_path.path)-1
        #Calculate the running sum
        sum += fn(sp.path[i], steps[i], sp.wiener_path.Δt, time_tracker) 
        #iterate time tracker
        time_tracker += sp.wiener_path.Δt
    end
    
    return sum
end

function refine_path(sp::StochasticPath)
    #=
    refines a stochastic path by refining the WP underling it and recalculating.
    
    #TODO: this needs some serious work. 
    Each function should take a stochastic path object and return a stochastic path
    I need to create some example functions
    =#
    wp = refine_path(sp.wiener_path)
    path = wp
    
    for fn in sp.functions_list
        path = fn(path)
    end
    
    return StochasticPath(
        wp
        ,sp.functions_list
        ,path
    )
end

function spread_path(sp::StochasticPath)
    #=
    =#
    wp = spread_path(sp.wiener_path)
    path = similar(wp)
    
    for fn in sp.functions_list
        path = fn(path)
    end
    
    return StochasticPath(
        wp
        ,sp.functions_list
        ,path
    )
end

#=
I also need to create a set of useful integrands.
=#

end #End of Module

module Examples
using LinearAlgebra, Statistics, Random, .SDE

export G


f(x,t) = x*t
g(x,t) = x
h(xt,Δwₜ,Δt,t) = xt*Δt + (xt*t)*Δwₜ

function G(sp::StochasticPath)
    return g(sp.path,0)
end

function H(sp::StochasticPath)
    return g(sp.path,0) #TODO: Finish this up somehow.
end

end

#module SDETesting
#=
This module is designed to contain tests of sdes.
=#
#using .SDE

# Testing constants
#const N1 = 1_500_000
#const N2 = 1_500
#const T = 1

#Example functions
#f1(x,t) = x

#Tests
#function t1() 
    #Check basic integral
#    w = WienerPath(N1,T)
#    @test integrate(f1,w) - 1/2*(w[N1]^2 - T) <= 1e2
#end

#end #end of module