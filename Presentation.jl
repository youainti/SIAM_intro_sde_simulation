### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 31c87e00-4cdc-11ec-1658-d1171bbae966
md"""
# Problem Description
In BLANK, XYZ proposes a deterministic model for the net present value of a 
satellite, conditional on the rate at which the satellite operator launches satellites, as:

$v = \int_0^T  \lambda e^{-\lambda t} \int_0^t pd^{-\rho \tau} d\tau  dt$


I am extending this to a stochastic price problem, where the instantaneous revenue (price of service) to the operator is a stochastic process.
Writing the problem in differential form, we get:

 - Prices $p$ are stochastic: $dp = \theta(\mu-p) dt + \sigma dW_t$
 - Discounted returns: $dr_\tau = p e^{-\rho \tau} d\tau$
 - Net present value depends on survival: $dv = r_t \cdot \lambda e^{-\lambda t} dt$



## Major questions
 1. What is the Expected Net Present Value $E(v)$?
 1. What is the distribution of net present values $v$?
 1. What happens when the price process is not a mean-reverting OU process?
    - Other mean reverting processes
        - Mean reverting geometric brownian motion
    - Non-mean reverting processes
        - Geometric brownian motion
        - Ornstein-Uhlenbeck motion


I am going to answer the first two questions analytically, and all three questions
numerically. 
In context, a standard goal in economics is to answer: **"What is the best policy when $\lambda$ is controllable, and how does it differ between a social planner and individual decision makers"**


## Exact solutions

### Exact Stochastic integral solution
Using standard Itô calculus, we get

### Expected Net Present Value
Using the Feynman-Kac formula, we get

note: are the discounting and survival interchangeable? 
If not, the best we can do with the feynman-kac formula
directly is put an upper bound on the expected net present value.
I believe they are interchangeable i.e. discounted expected value vs expected discounted value.
If so, then the suvival function is G and the feynman-kac formula solves directly.

### Distribution of Net Present Value
Using the (FKE or BKE?), we get

## Numerical simulations
Overview of Numerical simulation

"""

# ╔═╡ Cell order:
# ╠═31c87e00-4cdc-11ec-1658-d1171bbae966
