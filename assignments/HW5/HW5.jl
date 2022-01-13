### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 93723e67-fc2f-4394-80f0-d8918c5823ad
using LinearAlgebra, Plots,Statistics, ForwardDiff

# ╔═╡ 05713008-3c40-11ec-2e26-9b6306549123
begin
#=
Functions for managing significant digits
=#
	abstract type Significance end
	
	struct MachineLevel <: Significance end
	
	struct SignificantDigits <: Significance
	    level::Int
	end
	
	struct SignificantFigures <: Significance
	    level::Int
	end
	
	
	function round(x::Float64, d::SignificantFigures)
	    return Base.round(x, sigdigits=d.level)
	end
	
	function round(x::Float64, d::SignificantDigits)
	    return Base.round(x, digits=d.level)
	end
	
	function round(x::Float64, d::MachineLevel)
	    return x
	end
end

# ╔═╡ 2b12fc31-b479-42b3-a1de-8c90b8113294
begin
	#the rounding level used in most problems
	sig4 = SignificantDigits(4)
	sig5 = SignificantDigits(5)
	mcl = MachineLevel()
end

# ╔═╡ 49b41f5e-89e0-4e22-88fe-358ceedcd4d1
begin
	#=
	Way to keep track of numerical parameters
	=#
	abstract type NumericalBase end
	struct FundamentalNumericalBase <: NumericalBase
	    a::Function
	    Δt::Float64
	    sig::Significance #used for estimation precision
	end
	
	struct ExtendedNumericalBase <: NumericalBase
	    f::Function #known function
	    a::Function #from differential equation
	    Δt::Float64 #timestep
	    sig::Significance #used for estimation precision
	end
	#Add some other construction methods.
	
	
end

# ╔═╡ 6e0de80f-deee-4052-9969-1bf49960fc03
md"""
## PC-Exercise 8.1.1 : Euler Method

Using Machine level precision at two different levels of $\Delta t$.
"""

# ╔═╡ b166ac96-d934-4a95-809b-b9236a5ab356
begin
	#=
	Euler Methods
	=#
	function euler_step(fb::NumericalBase,z0,t)
	    return z0 + fb.Δt*fb.a(z0,t)
	end
	
	function euler_method(nb::NumericalBase, z0::Float64, time_start::Float64, time_stop::Float64)
	    
	    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
	    
	    Z = zeros(iterations) 
	    Z[1] = z0
	    t = time_start
	    
	    for i in 2:iterations
	        #println(i," ", iterations)
	        Z[i] = round(euler_step(nb,Z[i-1],t), nb.sig)
	        t+= nb.Δt
	    end
	    
	    return Z
	end
	
function euler_method_with_local_errors(
		nb::ExtendedNumericalBase
		, z0::Float64
		, time_start::Float64
		, time_stop::Float64
	)
	    
	    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
	    
	    X = zeros(iterations) #Local Errors
	    Y = zeros(iterations) #Actual Values
	    Z = zeros(iterations) #Full estimation
	    
	    
	    X[1] = z0
	    Y[1] = z0
	    Z[1] = z0
	    
	    t = time_start
	    
	    for i in 2:iterations
	        X[i] = round(euler_step(nb,Y[i-1],t), nb.sig)
	        Y[i] = nb.f(t)
	        Z[i] = round(euler_step(nb,Z[i-1],t), nb.sig)
	        t+= nb.Δt
	    end
	    
	    return X,Y,Z
	end
	
end

# ╔═╡ 760dcc6c-4656-45ab-a9b2-ae9a73a6e702
begin
	#model
	a8_1_1(x,t) = -5x;
	#analytic solution
	x8_1_1(t) = exp(-5t);
	
	#setup estimation parameters
	nb8_1_1_a = FundamentalNumericalBase(a8_1_1,2^-3,mcl);
	nb8_1_1_b = FundamentalNumericalBase(a8_1_1,2^-5,mcl);
end

# ╔═╡ d7619a16-e89b-4c73-bf62-153b0854b805
begin
	p = plot(0.0:nb8_1_1_b.Δt:0.99
		,[
			[x8_1_1(t) for t=0.0:nb8_1_1_b.Δt:0.99]
			,euler_method(nb8_1_1_b, 1.0,0.0,1.0)
		]
		,label = ["actual" "high resolution"]
	)
	plot!(p,euler_method(nb8_1_1_a, 1.0,0.0,1.0),0.0:nb8_1_1_a.Δt:0.99
		,label = "low resolution"
	)
end

# ╔═╡ 2d6437c2-9f34-4335-9f9b-59978ed94a1f
md"""
## PC-Exercise 8.1.2 : Global Error

For the IVP in PC-Exercise 8.1.1, calculate the global discretization error at time $t=1$ for the euler method with time steps of equal length $\Delta = 2^{-n}, n\in \{0,1,\dots,13\}$, rounding off to 5 digits. 
Plot the $\log_2$ of the errors against $\log_2 \Delta$ and determine the slope of the resulting curve.
"""

# ╔═╡ 4e1e8635-2c50-47be-a598-0e53b7b83245
begin
	#8_1_2: make a list of numerical bases
	deltas = [2.0^-i for i = 0:13]
	bases = [ExtendedNumericalBase(x8_1_1,a8_1_1,delta,sig5) for  delta = deltas];
	
	global_errors_812 = zeros(length(bases));
end

# ╔═╡ 6f7af5f7-5ca1-4221-a0e2-ee91bba52ce3
for (i,base) in enumerate(bases)

    #Calculate method
    x,y,z = euler_method_with_local_errors(base, 1.0, 0.0, 1.0)
    
    #issue
    global_errors_812[i] = last(y)-last(z)
    
end

# ╔═╡ 66ada403-5e02-46db-9a26-9557056ea726
global_errors_812a = log2.(abs.(global_errors_812))

# ╔═╡ 983c2710-236e-4eae-b8c6-cf2957fc4f3e


# ╔═╡ 7c568c35-4129-44d7-92f0-b38865e384b3
plot(log2.(deltas)[2:14]
    ,[global_errors_812a[2:14],log2.(deltas)[2:14]]
    ,labels = ["Global Errors" "Delta Line"]
	,legend=:topleft
)

# ╔═╡ d21ce1b4-c9ce-41ff-ba54-f0f3f1037250
md"""
## PC-Exercise 8.1.3
Repeat PC-Exercise 8.1.2 using Heun's method. Compare results with eulers method.
Use machine level precision.
"""

# ╔═╡ b6d99bcb-9dc5-4156-ab2f-1c6854555b3c
begin
	bases_813 = [ExtendedNumericalBase(x8_1_1,a8_1_1,delta,mcl) for  delta = deltas];
	
	global_errors_813 = zeros(length(bases));
end

# ╔═╡ 3b34e92a-6ffa-4ab8-8821-dd233eb2c692
begin
	
	function heuns_step(nb::NumericalBase,y0,t)
	    return y0 + 0.5 * ( nb.a(y0,t) + nb.a(euler_step(nb,y0,t),t)) * nb.Δt
	end
	
	function heuns_method_with_local_errors(
			nb::ExtendedNumericalBase
			, x0::Float64
			, time_start::Float64
			, time_stop::Float64
		)
	
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
	
	
end

# ╔═╡ 1da23c46-394c-4e15-ad93-8ea5b0fa66da
begin
	#find the list of estimates.
	
	global_errors2 = zeros(14)
	cumulative_local_errors_813 = zeros(14)
	
	for (i,base) in enumerate(bases_813)
	
	    #Calculate method
	    x,y,z = heuns_method_with_local_errors(base, 1.0, 0.0, 1.0)
	    
    	global_errors_813[i] = last(y)-last(z)
		cumulative_local_errors_813[i] = sum(x)
	end
end

# ╔═╡ 3779a360-528c-47f2-b475-59352ec9cb5b
global_errors_813a = log2.(abs.(global_errors_813))

# ╔═╡ 3d69ca06-bfc4-45ce-b3ba-e5e4f308bc92
-log2.(cumulative_local_errors_813)

# ╔═╡ 206959df-f13d-4901-bf37-78dbbcb2a422
plot(log2.(deltas)[2:14]
    ,[
		global_errors_813a[2:14]
		,log2.(deltas)[2:14]
		,-log2.(cumulative_local_errors_813)[2:14]
	]
    ,labels = ["Global Errors" "Delta Line" "cumulative local errors"]
	,legend=:topleft
)

# ╔═╡ 260527a1-0979-41cb-9a10-6f21a9920703
md"""
I decided to throw in the cumulative local errors just to see how it differers
"""

# ╔═╡ 4ae88583-3fcf-4737-8505-cc0e444aa948
md"""
## PC-Exercise 8.1.5 :Adams Brashford
Repeat PC-Exercise 8.1.3 using the 3-step Adams-Bashford Method with the Heun method as its starting routine.
This includes using machine level precision.
"""

# ╔═╡ 4eb6000e-8f62-45ac-a1ab-e4ad16e3bc72
begin
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
	
	    return y2 + nb.Δt * ( 23*nb.a(y2,t2) - 16*nb.a(y1,t1) + 5*nb.a(y0,t0))/12
	end
	
	function ab_method_with_local_errors(
			nb::ExtendedNumericalBase
			, z0::Float64
			, time_start::Float64
			, time_stop::Float64
		)
	
	    iterations = Int(ceil((time_stop-time_start)/nb.Δt))
	
	    X = zeros(iterations) #Local Errors
	    Y = zeros(iterations) #Actual Values
	    Z = zeros(iterations) #Full estimation
	
	
	    X[1] = z0
	    Y[1] = z0
	    Z[1] = z0
	
	    t1 = time_start
	    t2 = t1  + nb.Δt
	
	    #Prep the full estimation
	    Z[2] = round(euler_step(nb,Z[1],t1), nb.sig)
	    Z[3] = round(heuns_step(nb,Z[2],t2), nb.sig)
	
	    #prep the actual values
	    Y[2] = nb.f(t1)
	    Y[3] = nb.f(t2)
	
	    #Find the first few steps with local errors
	    X[2] = round(euler_step(nb,Y[1],t1), nb.sig)
	    X[3] = round(heuns_step(nb,Y[2],t2), nb.sig)
	
		#get starting time
	    t = t2+nb.Δt
		#begin iterations
	    for i in 4:iterations
	        X[i] = round(adam_brashford_step(nb,Y[i-1],Y[i-2],Y[i-3],t), nb.sig)
	        Y[i] = round(nb.f(t), nb.sig)
	        Z[i] = round(adam_brashford_step(nb,Z[i-1],Z[i-2],Z[i-3],t), nb.sig)
	        t+= nb.Δt
	    end
	
	    return X,Y,Z
	end
	
	
end

# ╔═╡ e0fef77e-8a9c-4f13-ae63-331aa952a29d
begin
	#find the list of estimates.
	global_errors_815 = zeros(14)
	
	for (i,base) in enumerate(bases_813)
		if i<4
			continue
		end
	    #Calculate method
	    x,y,z = ab_method_with_local_errors(base, 1.0, 0.0, 1.0)
	    
	    #Get final results
    	global_errors_815[i] = last(y)-last(z)
	end
end

# ╔═╡ 8436ab86-bc96-40bd-8939-5270b964b5a8
plot(log2.(deltas[4:14])
    ,[log2.(global_errors_815[4:14]) log2.(deltas[4:14])]
    ,labels = ["Global Errors" "Delta Line"]
	,legend=:topleft
)

# ╔═╡ f8433ea7-1a67-4dae-819e-ae839df1798e
md"""
## PC-Exercise 8.1.7
Compare the error of the Euler and Richardson/Romberg extrapolation approximations of X(1) for the solution of the initial value problem
$\frac{\partial x}{\partial t} = -x, x(0)=1$
for equal time steps $\Delta = 2^{-3}, \dots, 2^{-10}$.
Plot $\log_2$ of the errors against $\log_2 \Delta$.
"""

# ╔═╡ cefd8f6b-188e-440b-ac4b-6ba7a987f6fb
begin
	#=
	Richardson Extrapolation of euler method
	=#
	
	function richardson_extrapolation_of_euler(nb::NumericalBase, z0, start_time, end_time)
	    nb2 = FundamentalNumericalBase(nb.a,nb.Δt/2, nb.sig)

		#solve the low resolution euler method
		Y1 = euler_method(nb,z0, start_time, end_time)
		#solve the higher resolution euler method
		Y2 = euler_method(nb2,z0, start_time, end_time)

		#get the needed terms
		Z = [2*Y2[2*itt] - yn for (itt,yn) in enumerate(Y1)]

		return (Richardson=Z, LowEuler=Y1, HighEuler=Y2)
	end


end

# ╔═╡ 58593717-bdbd-47a0-83c1-ee0f79802e7e
begin
	#Values of the exact solution
	f8_1_7(t) = exp(-t);
	a8_1_7(x,t) = -x;
	
	exact_solutions_817 = [f8_1_7.(0:(2.0^(-n)):1) for n=3:10];
	
	time_steps = [(2.0^(-n)) for n=3:10];
	
	bases_817 = [ExtendedNumericalBase(f8_1_7,a8_1_7,2.0^(-n),mcl) for n=3:10];
end

# ╔═╡ d0271b52-fcd7-4fa2-a7bd-7228fa169857
begin
	global_errors_817 = (
		Richardson = zeros(length(bases_817))
		, LowEuler = zeros(length(bases_817))
		, HighEuler = zeros(length(bases_817))
	)
	
	for (i,base) in enumerate(bases_817)
	
	    #Calculate method
	    zyy = richardson_extrapolation_of_euler(base, 1.0, 0.0, 1.0)
	    
	    Z = last(zyy.Richardson)
		Xl = last(zyy.LowEuler)
		Xh = last(zyy.HighEuler)
		Y = last(exact_solutions_817[i])
	
		global_errors_817.Richardson[i] = Y-Z
		global_errors_817.LowEuler[i] = Y-Xl
		global_errors_817.HighEuler[i] = Y-Xh
	end
end

# ╔═╡ 9bbdabec-0444-4e0f-8f51-f1669a9e661b
#plot the errors
plot(log2.(time_steps)
    ,[log2.(time_steps) log2.(abs.(global_errors_817.Richardson)) log2.(abs.(global_errors_817.LowEuler)) log2.(abs.(global_errors_817.HighEuler)) ]
    ,labels = ["Δt" "Richardson" "Low resolution Euler" "high resolution euler"]
	,legend=:topleft
)

# ╔═╡ 752864fe-2822-415d-a9ee-a3e808f28124
md"""
Note that I took the $\log_2$ of the absolute value of the errors as they tended to be positive.
"""

# ╔═╡ 5b43aba2-cfb3-4a39-9f12-a9ef58b47c3c
md"""
## PC-Exercise 8.2.1
Use the 2nd Order truncated taylor method with equal length time steps $\Delta = 2^{-3}, \dots, 2^{-10}$
to calculate approximations to the solution:

$x(t)=\frac{2}{1+e^{-t^2}}$

of the initial value problem

$\frac{\partial x}{\partial t} = tx(2-x) ~~~~~~~~~ x(0)=1$ 

over the interval $0\leq t \leq 0.5$
Repeat the calculations using the 3rd order truncated taylor method.
Plot $\log_2$ of the global discretiation errors against $\log_2 \Delta$.
"""

# ╔═╡ 0fbd8b0a-7101-4484-b8bd-0ea46befaeb6
struct DifferentiatedNumericalBase <: NumericalBase
	f::Function #known function
	a::Dict{String, Function} #This contains functions representing the partial derivatives wrt time and space. Rows represent time and columns represent space. Thus A[2,2] is ȧ′
	Δt::Float64 #timestep
	sig::Significance #used for estimation precision
end

# ╔═╡ 38015486-ffb6-4d66-8144-266666294958
begin
	function truncated_taylor_2nd_order_step(dnb::DifferentiatedNumericalBase, z0, t)
		z1 = z0 + 
			dnb.a["a"](z0,t) * dnb.Δt + 
			(
				dnb.a["ȧ"](z0,t) + dnb.a["a′"](z0,t) * dnb.a["a"](z0,t)
			)*(dnb.Δt^2)/2
		return z1
	end

	
	function truncated_taylor_3rd_order_step(dnb::DifferentiatedNumericalBase, z0, t)
		z1 = truncated_taylor_2nd_order_step(dnb, z0, t) +
				(
					dnb.a["ȧ̇"](z0,t) + 
					2*dnb.a["ȧ′"](z0,t) * dnb.a["a"](z0,t) +
					dnb.a["a′′"](z0,t) * dnb.a["a"](z0,t)^2 +
					dnb.a["ȧ"](z0,t) * dnb.a["a′"](z0,t) +
					dnb.a["a′"](z0,t)^2 * dnb.a["a"](z0,t)
				)*(dnb.Δt^3)/6
		return z1
	end
end

# ╔═╡ 59210a08-2bfb-4a3a-b268-8d4929ff5d7b
function solve(
	dnb::DifferentiatedNumericalBase
	,z0::Float64
	,time_start::Float64
	,time_stop::Float64
	,step_algorithm::Function
)
	#=
	This is generic across taylor algorithms (assuming dnb.a has enough dimensions/derivatives)
	=#
	
	iterations = Int(ceil((time_stop-time_start)/dnb.Δt))

	Z = zeros(iterations)

	Z[1] = z0
	t = time_start
	
	for i in 2:iterations
		Z[i] = step_algorithm(dnb,Z[i-1],t)
		t=+ dnb.Δt
	end

	return Z
end

# ╔═╡ 13b340e5-2eec-474e-9d5d-1493bbaee1a3
begin
	#get the list of derivatives that are required as functions
	f(x,t)=2/(1+exp(-t^2))
	f(t)=f(nothing,t)
	a(x,t) = t*x*(2-x)

	#time partials
	ȧ(x,t) = x*(2-x)
	ȧ̇(x,t) = 0
	
	#space partials
	a′(x,t) = t*(2-2x)
	a′′(x,t) = -2*t

	#mixed partials
	ȧ′(x,t) = 2 - 2x

	A = Dict(
		"a" => a
		,"ȧ" => ȧ
		,"ȧ̇" => ȧ̇
		,"a′" => a′
		,"a′′" => a′′
		,"ȧ′" => ȧ′
	)

end

# ╔═╡ 96889dbd-2306-40ff-bf40-34025c140f38
begin
	#setup deltas
	deltas_821 = [2.0^-x for x in 3:10]
	#setup numerical bases

	bases_821 = [DifferentiatedNumericalBase(f,A,x,mcl) for x in deltas_821]

	#setup error vectors
	global_errors_821a = zero(deltas_821)
	global_errors_821b = zero(deltas_821)

	#get the actual result
	actual_result_821 = f(nothing,0.5)
	
	for (i,base) in enumerate(bases_821)
		z = solve(base,1.0,0.0,0.5,truncated_taylor_2nd_order_step)
		global_errors_821a[i] = actual_result_821 - last(z)

		zz = solve(base,1.0,0.0,0.5,truncated_taylor_3rd_order_step)
		global_errors_821b[i] = actual_result_821 - last(zz)
	end
end

# ╔═╡ eae6db4d-8cab-4876-94aa-0fa4f806a4b9
plot(log2.(deltas_821)
    ,[log2.(deltas_821) log2.(abs.(global_errors_821a)) log2.(abs.(global_errors_821b)) ]
    ,labels = ["Δt" "2nd order taylor" "3rd order taylor"]
	,legend=:topleft
)

# ╔═╡ 98e693de-d244-4f56-8c43-bc5169948dd9
md"""
Honestly, I'm not sure what mistake I made here. It *should* work, but something is off. I suspect it is part of the 2nd order system, but I honestly don't know.
"""

# ╔═╡ 7ba7f265-92ec-44b7-bba4-3a3930cda6d4
md"""
## PC-Exercise 8.2.2
Repeat PC-Exercise 8.2.1 using the 4th order Runge-Kutta method with equal length time steps $\Delta = 2^{-3}, \dots, 2^{-7}$
"""

# ╔═╡ 45d55719-da29-4f95-9817-e2e531aa5b7b
function RK4_step(nb::NumericalBase,z,t)
	k₁ = nb.a(z,t)
	k₂ = nb.a(z + nb.Δt * k₁/2, t + nb.Δt/2)
	k₃ = nb.a(z + nb.Δt * k₂/2, t + nb.Δt/2)
	k₄ = nb.a(z + nb.Δt * k₃, t + nb.Δt)
	
	
	zₙ₊₁ = z + nb.Δt * (k₁ + 2*k₂ + 2*k₃ + k₄)/6

	return zₙ₊₁
end

# ╔═╡ 775166f6-f4ee-489b-9ef8-fbae8ec867df
function RK4(
	nb::NumericalBase
	,z0
	,start_time
	,stop_time
)
	iterations = Int(ceil((stop_time-start_time)/nb.Δt))

	Z = zeros(iterations)
	Z[1] = z0
	t = start_time
	
	for i in 2:iterations
		Z[i] = RK4_step(nb,Z[i-1],t)
		t=+nb.Δt
	end

	return Z
end

# ╔═╡ 85a36dcf-cf8e-4bcb-bc61-d6b424df0a15
#just comparing methods
RK4(bases_813[13],1.0,0.0,1.0) .- euler_method(bases_813[13],1.0,0.0,1.0)

# ╔═╡ a391166d-ab73-4519-935d-dcaf4a1d2ee4
bases_822 = [ExtendedNumericalBase(f,a,x,mcl) for x in deltas_821]

# ╔═╡ 73f15bfd-9397-4558-82ac-591bc58ffbf8
begin

	#setup error vectors
	global_errors_822 = zero(deltas_821)
	euler_error_822 = zero(deltas_821)

	#get the actual result
	actual_result_822 = f(nothing,0.5)
	
	for (i,base) in enumerate(bases_822)
		z = RK4(base,1.0,0.0,0.5)
		global_errors_822[i] = actual_result_822 - last(z)

		euler_error_822[i] = actual_result_822 - last(euler_method(base,1.0,0.0,0.5))

	end
end

# ╔═╡ b5e44e8d-5733-409b-b1cf-b4e10a174665
RK4(bases_822[8],1.0,0.0,0.5) .- euler_method(bases_822[8],1.0,0.0,0.5)

# ╔═╡ 37cfee15-12e9-41c1-a9a7-b813d7d0d6a8
log2.(euler_error_822)

# ╔═╡ 1a6db083-dc19-4ee8-83b8-25bbda2080bd
plot(log2.(deltas_821)
    ,[log2.(deltas_821) log2.(abs.(global_errors_822)) log2.(euler_error_822) ]
    ,labels = ["Δt" "RK4 global error" "euler error (for comparison)"]
	,legend=:topleft
)

# ╔═╡ 1e9c7f61-c7fc-428e-9c6b-cbc6a4819dc5
md"""
As you can see, there is some issue with my RK4 and taylor implementation, as the euler method performs significantly better. 
I have not been able to track down what is going on.
"""

# ╔═╡ 5d091b28-e5d3-4d38-a625-f3cb92b65fe8
plot([RK4(bases_822[8],1.0,0.0,0.5) euler_method(bases_822[8],1.0,0.0,0.5) ])

# ╔═╡ 2d37049d-0174-4a69-9ecb-1d4d233dfecb
plot([f(t) for t in 0.0:deltas_821[8]:0.5])

# ╔═╡ 9fb405a3-61dd-4342-bc71-766faab0c089
md"""
## PC-Exercise 8.2.3
"""

# ╔═╡ f2e1b3f7-a6c5-4426-bd55-188bb31e2a66
begin
	function midpoint_step(nb::NumericalBase,z0,z1,t)
		return z0 + 2*nb.a(z1,t) * nb.Δt
	end
	
	function midpoint_method(nb::NumericalBase,z0,start_time,stop_time)
		iterations = Int(ceil((stop_time-start_time)/nb.Δt))
	
		Z = zeros(iterations)
		Y = zeros(iterations)
		X = zeros(iterations)
		
		Z[1] = z0
		Y[1] = z0
		X[1] = z0
		t = start_time
		
		#starting step
		Z[2] = euler_step(nb,Z[1],t)
		X[2] = Z[2]
	
		
		for i in 3:iterations
			Z[i] = midpoint_step(nb,Z[i-2],Z[i-1],t)
			Y[i] = nb.f(t)
			X[i] = midpoint_step(nb,Y[i-2],Y[i-1],t)
			
			t=+nb.Δt
		end
	
		return Z,Y,X
	end
end

# ╔═╡ 5e4e3b67-4354-4454-ac31-6b231802457b
#problem setup
nb_823 = ExtendedNumericalBase(
	t -> 2/3 * exp(-3*t) + 1/3
	,(x,t) -> -3*x + 1
	,0.1
	,mcl)

# ╔═╡ aa37da4d-565c-47d7-8d3a-e4fc482d063d
mid_results, actual1, mid_error = midpoint_method(nb_823,1.0,0.0,1.0)

# ╔═╡ 7060af6b-ca10-4d58-9b77-0f82459fcefd
eul_results, actual2, eul_error = euler_method_with_local_errors(nb_823,1.0,0.0,1.0)

# ╔═╡ 6ca93e0d-a862-4e7e-a4aa-a829776404a5
plot([mid_results eul_results])

# ╔═╡ e2e1cbc9-d0ab-4893-8e0a-c9c7eb095e01
md"""
There is quite a difference between the two, particularly in the nature of convergence.
"""

# ╔═╡ e341c000-6243-45e6-9592-51646c300bff
md"""
## PC-Exercise 8.4.1
Calculate 300 iterates of

with initial value $y_0 = 0.1$ using the prescribed arithematic of the PC, 
at each step rounding the value of $Y_{n+1}$ to the first 4 significant figures.
Plot the relative frequencies of the roundoff errors in a histogram on 
using 40 equal bins.
"""

# ╔═╡ 1e7feaa4-18f5-438d-92e4-5a6760370bfc
function iter8_4_2(yn::Float64)
    return π/3 * yn
end

# ╔═╡ a0efd6fd-4dc6-4728-9275-9f421a2a3ca2
begin
	step_y = 0.1
	N8_4_1 = 300
	rounding_errors = zeros(N8_4_1)
	
	for i in 2:N8_4_1
	    yi = iter8_4_2(step_y)
	    rounding_errors[i] = round(yi,MachineLevel()) - round(yi, sig4)
	    step_y = yi
	end
end

# ╔═╡ a886c490-7710-4bbd-afbf-68c770266413
md"""
Mean: $(mean(rounding_errors))

Stdev: $(std(rounding_errors))
"""

# ╔═╡ 58178a7d-4f2b-4239-a8e2-cea744c9c415
histogram(rounding_errors,bins=40)

# ╔═╡ eb845364-f78b-448f-9d0e-6db8d6bee4ba
md"""
## PC-Exercise 8.4.2
Use the Euler method with equal length time steps $\Delta = 2^{-2}$ for the differential equation 

$\partial x = x \partial t$

over the interval $0 \leq t \leq 1$ for 1,000 different initial values $x(0) \in [0.4,0.6]$.

Use both the prescribed arithematic and round to 4 decimal places to determine the final 
accumulative roundoff error in each case,
plotting the roundoff error in a histogram with 40 equal subintervals.
In addition, calculate the sample mean and sample variance.
"""

# ╔═╡ d06c65e6-5921-4966-8ab6-d3bf54413218
begin
	#record the known solution and get the differential version.
	f8_4_2(x) = exp(x);
	f8_4_2′(x,t) = x;
	
	N8_4_2 = 1000
	
	nb8_4_2 = ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-2,MachineLevel());
	x0 = (2*rand(N8_4_2).+4.0)/10;
	R = zeros(N8_4_2);
end

# ╔═╡ 53f8b28a-4c72-4215-b2d9-76551fd64dfe
function record8_4(nb::NumericalBase,x0,N)
    for i in 1:N
		#FIX: something here needs to change
        yi = last(euler_method(nb,x0[i],0.0,1.0))
        R[i] = yi - round(yi, sig4)
    end
    
    return R
end

# ╔═╡ 894ebec6-5640-4b61-b3c0-f28e1a9da1a2
err_x = record8_4(nb8_4_2,x0,N8_4_2)

# ╔═╡ f7e1b876-e329-42da-afab-ce8cb8d2a130
md"""
Mean:  $(mean(err_x))

stdev: $(std(err_x))
"""

# ╔═╡ 9b29cefd-0176-4e7f-b921-56398f53a7d0
histogram(R,bins=40)

# ╔═╡ f2aabd6c-b8f7-459b-8daf-d01ad60cfbbd
md"""
## PC-Exercise 8.4.3
Repeat PC-Exercise 8.4.2 with N = 200 and with equal length time steps $\Delta = 2^{-3}, \dots, 2^{-5}$
, determine the roundoff error in each case. plot the 90% confidence intervales for the mean value of the error against $\Delta$.
"""

# ╔═╡ 65a9af69-fb9b-4a9b-b1ab-885f5339c43e
begin
	nb8_4_3 = [
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-2,mcl),
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-3,mcl),
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-4,mcl),
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-5,mcl)]

	rnb8_4_3 = [
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-2,sig4),
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-3,sig4),
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-4,sig4),
		ExtendedNumericalBase(f8_4_2,f8_4_2′,2^-5,sig4)
	]
end;

# ╔═╡ c3f8aaab-22ed-445d-9f8d-07e7e44bb3d4
X = 0.4 .+ 0.2.*rand(200)

# ╔═╡ 70d14706-4de7-4f29-883b-4709fa6ae6de
x=0.5

# ╔═╡ 0bd4b413-4fb4-40d6-8227-d31ef346ebcf
euler_method_with_local_errors(nb8_4_3[1],X[3],0.0,1.0)

# ╔═╡ 629d4be7-5f7b-4e19-9ebb-3379dc4f05c0
begin
	d = Dict()
	for (n,r) in zip(nb8_4_3,rnb8_4_3)
		roundoff_error = []
		for x in X
			n1,n2,n3 = euler_method_with_local_errors(n,x,0.0,1.0)
			r1,r2,r3 = euler_method_with_local_errors(r,x,0.0,1.0)
	
			append!(roundoff_error, n3-r3) #todo: Calculate roundoff error
		end
	
		#process roundoff error
		d[n.Δt] = roundoff_error
	end

end

# ╔═╡ 412c196c-f3a7-49f5-af49-bf3eff226184
begin
	top = []
	bottom = []
	delts = [2.0^(-2) 2.0^(-3) 2.0^(-4) 2.0^(-5)]
	for s in delts
		bott,tops = quantile(d[s],[0.05,0.95])
		append!(top,tops)
		append!(bottom,bott)
	end
end

# ╔═╡ 35f239b0-8581-4917-9f7a-5895ac055cb4
bottom,top

# ╔═╡ bc1455a5-e844-4e10-b71d-125375323f76
plot([top bottom],label=["95% quantile" "5% quantile"], type=:line)

# ╔═╡ 5bba4c9f-210d-4897-840a-281df5e87614
md"""
each integer point above represents a single step along $delts.
It turns out the effect of roundoff error grows as you take more steps, i.e. round more often.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
ForwardDiff = "~0.10.24"
Plots = "~1.23.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "3533f5a691e60601fe60c90d8bc47a27aa2907ec"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2b72a5624e289ee18256111657663721d59c143e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.24"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "d189c6d2004f63fd3c91748c458b09f26de0efaa"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.61.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d911b6a12ba974dabe2291c6d450094a7226b372"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.1"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "ca7d534a27b1c279f05cd094196cb70c35e3d892"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.23.2"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═93723e67-fc2f-4394-80f0-d8918c5823ad
# ╠═05713008-3c40-11ec-2e26-9b6306549123
# ╠═2b12fc31-b479-42b3-a1de-8c90b8113294
# ╠═49b41f5e-89e0-4e22-88fe-358ceedcd4d1
# ╟─6e0de80f-deee-4052-9969-1bf49960fc03
# ╠═b166ac96-d934-4a95-809b-b9236a5ab356
# ╠═760dcc6c-4656-45ab-a9b2-ae9a73a6e702
# ╠═d7619a16-e89b-4c73-bf62-153b0854b805
# ╟─2d6437c2-9f34-4335-9f9b-59978ed94a1f
# ╠═4e1e8635-2c50-47be-a598-0e53b7b83245
# ╠═6f7af5f7-5ca1-4221-a0e2-ee91bba52ce3
# ╠═66ada403-5e02-46db-9a26-9557056ea726
# ╠═983c2710-236e-4eae-b8c6-cf2957fc4f3e
# ╠═7c568c35-4129-44d7-92f0-b38865e384b3
# ╟─d21ce1b4-c9ce-41ff-ba54-f0f3f1037250
# ╠═b6d99bcb-9dc5-4156-ab2f-1c6854555b3c
# ╠═3b34e92a-6ffa-4ab8-8821-dd233eb2c692
# ╠═1da23c46-394c-4e15-ad93-8ea5b0fa66da
# ╠═3779a360-528c-47f2-b475-59352ec9cb5b
# ╠═3d69ca06-bfc4-45ce-b3ba-e5e4f308bc92
# ╠═206959df-f13d-4901-bf37-78dbbcb2a422
# ╠═260527a1-0979-41cb-9a10-6f21a9920703
# ╠═4ae88583-3fcf-4737-8505-cc0e444aa948
# ╠═4eb6000e-8f62-45ac-a1ab-e4ad16e3bc72
# ╠═e0fef77e-8a9c-4f13-ae63-331aa952a29d
# ╠═85a36dcf-cf8e-4bcb-bc61-d6b424df0a15
# ╠═8436ab86-bc96-40bd-8939-5270b964b5a8
# ╟─f8433ea7-1a67-4dae-819e-ae839df1798e
# ╠═cefd8f6b-188e-440b-ac4b-6ba7a987f6fb
# ╠═58593717-bdbd-47a0-83c1-ee0f79802e7e
# ╠═d0271b52-fcd7-4fa2-a7bd-7228fa169857
# ╠═9bbdabec-0444-4e0f-8f51-f1669a9e661b
# ╟─752864fe-2822-415d-a9ee-a3e808f28124
# ╟─5b43aba2-cfb3-4a39-9f12-a9ef58b47c3c
# ╠═0fbd8b0a-7101-4484-b8bd-0ea46befaeb6
# ╠═38015486-ffb6-4d66-8144-266666294958
# ╠═59210a08-2bfb-4a3a-b268-8d4929ff5d7b
# ╠═13b340e5-2eec-474e-9d5d-1493bbaee1a3
# ╠═96889dbd-2306-40ff-bf40-34025c140f38
# ╠═eae6db4d-8cab-4876-94aa-0fa4f806a4b9
# ╠═98e693de-d244-4f56-8c43-bc5169948dd9
# ╠═7ba7f265-92ec-44b7-bba4-3a3930cda6d4
# ╠═45d55719-da29-4f95-9817-e2e531aa5b7b
# ╠═775166f6-f4ee-489b-9ef8-fbae8ec867df
# ╠═a391166d-ab73-4519-935d-dcaf4a1d2ee4
# ╠═73f15bfd-9397-4558-82ac-591bc58ffbf8
# ╠═b5e44e8d-5733-409b-b1cf-b4e10a174665
# ╠═37cfee15-12e9-41c1-a9a7-b813d7d0d6a8
# ╠═1a6db083-dc19-4ee8-83b8-25bbda2080bd
# ╟─1e9c7f61-c7fc-428e-9c6b-cbc6a4819dc5
# ╠═5d091b28-e5d3-4d38-a625-f3cb92b65fe8
# ╠═2d37049d-0174-4a69-9ecb-1d4d233dfecb
# ╠═9fb405a3-61dd-4342-bc71-766faab0c089
# ╠═f2e1b3f7-a6c5-4426-bd55-188bb31e2a66
# ╠═5e4e3b67-4354-4454-ac31-6b231802457b
# ╠═aa37da4d-565c-47d7-8d3a-e4fc482d063d
# ╠═7060af6b-ca10-4d58-9b77-0f82459fcefd
# ╠═6ca93e0d-a862-4e7e-a4aa-a829776404a5
# ╟─e2e1cbc9-d0ab-4893-8e0a-c9c7eb095e01
# ╠═e341c000-6243-45e6-9592-51646c300bff
# ╠═1e7feaa4-18f5-438d-92e4-5a6760370bfc
# ╠═a0efd6fd-4dc6-4728-9275-9f421a2a3ca2
# ╟─a886c490-7710-4bbd-afbf-68c770266413
# ╠═58178a7d-4f2b-4239-a8e2-cea744c9c415
# ╟─eb845364-f78b-448f-9d0e-6db8d6bee4ba
# ╠═d06c65e6-5921-4966-8ab6-d3bf54413218
# ╠═53f8b28a-4c72-4215-b2d9-76551fd64dfe
# ╠═894ebec6-5640-4b61-b3c0-f28e1a9da1a2
# ╠═f7e1b876-e329-42da-afab-ce8cb8d2a130
# ╠═9b29cefd-0176-4e7f-b921-56398f53a7d0
# ╟─f2aabd6c-b8f7-459b-8daf-d01ad60cfbbd
# ╠═65a9af69-fb9b-4a9b-b1ab-885f5339c43e
# ╠═c3f8aaab-22ed-445d-9f8d-07e7e44bb3d4
# ╠═70d14706-4de7-4f29-883b-4709fa6ae6de
# ╠═0bd4b413-4fb4-40d6-8227-d31ef346ebcf
# ╠═629d4be7-5f7b-4e19-9ebb-3379dc4f05c0
# ╠═412c196c-f3a7-49f5-af49-bf3eff226184
# ╠═35f239b0-8581-4917-9f7a-5895ac055cb4
# ╠═bc1455a5-e844-4e10-b71d-125375323f76
# ╠═5bba4c9f-210d-4897-840a-281df5e87614
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
