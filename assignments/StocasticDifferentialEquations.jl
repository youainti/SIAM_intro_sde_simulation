
#This is a module with useful functions for working with stochastic processes
module StochasticDifferentialEquations

### FUNCTIONS ###

#This function refines a random walk
function refine_walk(walk)
    len = length(walk)
    arr = zeros(2*len-1)
    
    for ind in 1:len-1
        arr[2*ind-1] = walk[ind]
        arr[2*ind] =  (0.5 * (walk[ind] + walk[ind+1]) .+ (((√Δt)/2)*randn(1)))[1]
    end
    arr[2*len-1] = walk[len]
    
    return arr
end


### STRUCTS ###

#Create a struct to capture the type of of a 1d box, with a couple of helpful methods
struct ExitType
    time::Int
    exited_top::Bool
end

measured_time(exit_type::ExitType) = exit_type.time
exited_top_boundary(exit_type::ExitType) = exit_type.exited_top
exited_bottom_boundary(exit_type::ExitType) = !exit_type.exited_top



end #End of Module