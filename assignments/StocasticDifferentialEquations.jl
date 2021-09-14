
#This is a module with useful functions for working with stochastic processes
module StochasticDifferentialEquations

#=
Create a struct to capture the type of exit from a 1d box, with a couple of helpful methods
=#
struct ExitType1d
    time::Int
    exited_top::Bool
end

ExitType = Union[ExitType1d]

measured_time(exit_type::ExitType) = exit_type.time
exited_top_boundary(exit_type::ExitType) = exit_type.exited_top
exited_bottom_boundary(exit_type::ExitType) = !exit_type.exited_top



#=
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
    WienerPath(p,d) = p[1] != 0.0 ? new(vcat(0.0,p),d) : new(p,d)
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

function refine(p::WienerPath)
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
        arr[2*ind] =  (0.5 * (p.path[ind] + p.path[ind+1]) 
    end
    arr[2*len-1] = p.path[len]
    
    return WienerPath(arr, p.Δt/2)
end

function integrate(f::Function, wp::WienerProcess)
    #=
    This function integrates across the time domain, by computing the riemann sum
    associated with the appropriate ito integral.
        
    There is an alternate formulation using dot products, but I had already written this.
        
    it returns a Float64 value for the integral
    =#
        
    sum = 0
    steps = get_steps(wp)
    
    for i in 1:length(wp.path)-1
        sum += f(wp.path[i]) * (steps[i]) 
    end
    
    return sum
end
    

end #End of Module