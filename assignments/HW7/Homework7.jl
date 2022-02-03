### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 29924739-9883-4588-8262-2198a22f4cde
begin
	abstract type StochasticPaths end
	struct StochasticPath <: StochasticPaths
		path::Array{Float64}
		Δt::Float64
		t0::Float64
	end
	Base.getindex(wp::StochasticPath, i::Int) = wp.path[i]
	Base.length(wp::StochasticPath) = length(wp.path)
	Base.lastindex(p::StochasticPath) = lastindex(p.path)
end

# ╔═╡ cde68a0a-70e9-11ec-0612-79b52a798105
begin
abstract type StochasticProcesses end
#=
Stochastic Processes contain information and a generator for stochastic paths
=#
struct WienerProcess <: StochasticProcesses
	Δt::Float64
end
	
struct StochasticProcess <: StochasticProcesses
	Δt::Float64
	functions::Array{Function}
end



#=
Stochastic paths represent realizations of stochastic processes
=#

struct WienerPath <: StochasticPaths
	path::Array{Float64}
	Δt::Float64
end
Base.getindex(wp::WienerPath,i::Int) = wp.path[i]
Base.length(wp::WienerPath) = length(wp.path)
Base.lastindex(p::WienerPath) = lastindex(p.path)

function (wp::WienerProcess)(N::Int) 
	WienerPath(cumsum(sqrt(wp.Δt).*randn(N)), wp.Δt)
end

function steps(path::StochasticPaths)
	l = length(path)
	s = zeros(l-1)
	for i in 1:l-1
		s[i] = path[i+1]-path[i]
	end
	return s
end

end#begin

# ╔═╡ b3c670e2-2da1-4b6a-9e35-f3fca129206a
begin
	#Euler Maruyama method

	abstract type SDE end #represents the specification of the SDE
	struct SemiLinearMotion <: SDE
		μ::Function
		σ::Function
		σ′::Function
	end

	abstract type SDESemiLinearMethod end #represents the solution method
	#Each of the following structs contains the specifications required to run the associated solve function.
	struct EulerMaruyama <: SDESemiLinearMethod
		Δt::Float64
	end
	struct Milstein <: SDESemiLinearMethod
		Δt::Float64
	end
end

# ╔═╡ d8955067-6e28-4e1d-bab7-79fb41765e0d
begin
	function step(path::StochasticPath, em::EulerMaruyama, sde::SemiLinearMotion, step_number::Int)
		#the path gets updated
		current_step = path[step_number]
		current_time = step_number*em.Δt .+ path.t0
		wiener_step = sqrt(em.Δt) * randn()

		next_step = current_step + sde.μ(current_step,current_time) * em.Δt + sde.σ(current_step,current_time) * wiener_step

		return next_step
	end

	
	function step(path::StochasticPath, em::Milstein, sde::SemiLinearMotion, step_number::Int)
		current_step = path[step_number]
		current_time = step_number*em.Δt .+ path.t0
		wiener_step = sqrt(em.Δt) * randn()

		next_step = current_step +
			sde.μ(current_step,current_time) * em.Δt + 
			sde.σ(current_step,current_time) * wiener_step +
			sde.σ(current_step,current_time) * sde.σ′(current_step,current_time)*
				(wiener_step^2 - em.Δt)/2

		return next_step
	end
end

# ╔═╡ 5a711ec0-0e18-44ab-830d-1f6dd7e3b74e
function iterate(
	sde::SemiLinearMotion
	,method::SDESemiLinearMethod
	,start_time::Float64
	,n_steps::Int
	,starting_value::Float64
)
	#create a stochastic path with the required starting values
	p = zeros(n_steps+1)
	p[1] = starting_value
	path = StochasticPath(p,method.Δt,start_time)
	
	#take the required number of steps
	for i in 1:n_steps
		path.path[i+1] = step(path, method, sde, i)
	end

	return path
end

# ╔═╡ 9b2ab708-9186-469f-8850-1834679cd7ea
struct IncorrectBounds <: Exception
	#just because I'm likely to make mistakes
	message::String
	upper::Float64
	lower::Float64
end

# ╔═╡ 376a54f7-6d59-4794-87d7-8c8e211e69d2
function solve(
	sde::SemiLinearMotion
	,method::SDESemiLinearMethod
	,start_time::Float64
	,stop_time::Float64
	,starting_value::Float64
	;number_iterations::Int = 1
)
	#check bounds
	if start_time >= stop_time
		throw(IncorrectBounds("Stop time less than or equal to start time", stop_time, start_time))
	end

	#get number of steps
	n_steps = Int(ceil((stop_time - start_time)/method.Δt))

	if number_iterations == 1
		return iterate(sde,method,start_time,n_steps,starting_value)		
	else
		#Parallelize....

		#preallocate results array
		results = Array{StochasticPath}(undef,number_iterations)

		#Assign tasks to threads and save to results.
		Threads.@threads for i in 1:number_iterations
			results[i] = iterate(sde,method,start_time,n_steps,starting_value)
		end

		return results
	end
	
end

# ╔═╡ c34a054a-be38-4b50-a7cd-78e318d863a2
begin #setup the problem
	GeometricBrownianMotion(μ,σ) = SemiLinearMotion(
		(x,t) -> μ*x
		,(x,t) -> σ*x
		,(x,t) -> σ
	)
	gmb = GeometricBrownianMotion(0.2,0.6)
	t0 = 0.0
	x0 = 1.0
	n = 10
end

# ╔═╡ d4233b44-02e0-472c-b20f-fa6acbf88067
a = solve(gmb, EulerMaruyama(0.01), t0,1.0,x0)

# ╔═╡ fa822ee1-c345-4d01-aaf9-1f3ea0df88a3
b = solve(gmb, Milstein(0.02), t0,1.0,x0, number_iterations=1_000_000)

# ╔═╡ ee3b3365-78b8-49e4-a04f-4071c942afce
last.(b)

# ╔═╡ 45b55800-84d8-4c32-be99-c1d6bef67f8f
function SampleMeanSquareConvergence(X::Array{StochasticPath})
	return sum(last.(b).^2)
end

# ╔═╡ d019935a-ee8d-4d5e-979a-6c9908191e4b
SampleMeanSquareConvergence(b)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╠═cde68a0a-70e9-11ec-0612-79b52a798105
# ╠═29924739-9883-4588-8262-2198a22f4cde
# ╠═b3c670e2-2da1-4b6a-9e35-f3fca129206a
# ╠═d8955067-6e28-4e1d-bab7-79fb41765e0d
# ╠═5a711ec0-0e18-44ab-830d-1f6dd7e3b74e
# ╠═9b2ab708-9186-469f-8850-1834679cd7ea
# ╠═376a54f7-6d59-4794-87d7-8c8e211e69d2
# ╠═c34a054a-be38-4b50-a7cd-78e318d863a2
# ╠═d4233b44-02e0-472c-b20f-fa6acbf88067
# ╠═fa822ee1-c345-4d01-aaf9-1f3ea0df88a3
# ╠═ee3b3365-78b8-49e4-a04f-4071c942afce
# ╠═45b55800-84d8-4c32-be99-c1d6bef67f8f
# ╠═d019935a-ee8d-4d5e-979a-6c9908191e4b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
