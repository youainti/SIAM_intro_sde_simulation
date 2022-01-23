# Problem 5.2

$dX_t = (a_1 X_t + a_2)dt + g(X_t,t) dW_t$

in integral form this is

$X(t) - X(0) = \int^t_0(a_1 X_s + a_2)ds + \int^t_0g(X_s,s) dW_s$

Thus

$E(X(t) - X(0)) = E \int^t_0(a_1 X_s + a_2)ds + E(\int^t_0g(X_s,s) dW_s)$

By the martingale property and linearity of expectation and integration

$E(X(t) )- E( X(0)) =E(X(t) - X(0)) = \int^t_0(a_1 E(X_s) + a_2)ds + 0$

By substituting $m_t = E(X_t)$ we get

$m_t - m_0 = \int^t_0 a_1 m_s + a_2 ds$

Which in differential form is the initial value problem

$\frac{d m_s}{ds} = a_1 m_s + a_2$        $m_0 = E(X_0)$

I am terrible at solving ODE's but using guess and check (based on the general form given in the book),
the solution takes the form 

$m(s) = b + e^{ds}$, thus $m^\prime(s) = d e^{ds}$

If that is true, then




# Problem 6.4



# Problem 6.5

Consider the SDE

$dS_t = dW_t$, i.e. $S_t = W_t$

Then for $Y_t(x) = \phi(x,t) = \frac{1}{1-x}$ by the ito formula, we get that

$dY_t = 0 dt + \phi_x(S_t,t) dS_t + \frac{1}{2} 1^2 \phi_{xx}(S_t,t) dt$

Thus as $\phi_x = (1-x)^{-2}$ and $\phi_{xx}= 2(1-x)^{-3}$ we get 


$dY_t = (1-S_t)^{-2}dS_t + \frac{2}{2}(1-S_t)^{-3} dt$

# Problem 7.1
## Expectation
Prove
$E\left[  \sum q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i}) ^2  -0.5  \sum q(W_{t_i})   \delta t\right] = 0$

Note that due to independence of $q(W_{t_i})$ and $(W_{t_{i+0.5}} - W_{t_i})$

$\sum E\left[ q(W_{t_i})\right]  E\left[ (W_{t_{i+0.5}} - W_{t_i}) ^2\right]  - 0.5  \delta t E\left[ \sum q(W_{t_i})  \right]$

And note that $(W_{t_{i+0.5}} - W_{t_i}) \sim N(0,\frac{1}{\sqrt 2} \delta t)$, thus

$\sum E\left[ q(W_{t_i})\right] 0.5 \delta t  - 0.5  \delta t E\left[ \sum q(W_{t_i})  \right] = 0$

## Variance
Prove 
$V\left[  \sum q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i}) ^2  -0.5  \sum q(W_{t_i})   \delta t\right]$ is of order $O(\delta t)$

Note the three major terms from factoring
$E\left[ ( \sum q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i}) ^2  -0.5  \sum q(W_{t_i})   \delta t)^2 \right]$ 

$$
A = \sum_i q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i}) ^2 \times \sum_j q(W_{t_j})(W_{t_{j+0.5}} - W_{t_j})^2 \\
= 2 \sum_{i < j} \sum_j q(W_{t_j})q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i})^2(W_{t_{j+0.5}} - W_{t_j})^2
+ \sum_i q(W_{t_i})^2 (W_{t_{i+0.5}} - W_{t_i})^4
$$ 


$B = -2 \sum_i q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i}) ^2   \times 0.5  \sum_j q(W_{t_j})   \delta t$


$C = 0.25  \sum_i q(W_{t_i})  \sum_j q(W_{t_j}) \delta t^2$




# Problem 7.7 
# Problem PC 7.1
# Problem PC 7.2

