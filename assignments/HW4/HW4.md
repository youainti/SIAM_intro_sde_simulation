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

$\frac{d m_s}{ds} = a_1 m_s + a_2$    
with initial condition
$m_0 = E(X_0)$

I am terrible at solving ODE's, but Referencing Lessons 11A and 11B on integrating factors from Tenenbaum and Pollard's ODE's 
I get:

$$
\frac{\partial m_t}{\partial t} -a_1 m_t =a_2
$$

With integrating factor $e^{-a_1 t}$ we get 

$$
m_t = \frac{-a_2}{a_1}\left(1 - e^{a_1 t} \right) + m_0 e^{a_1 t}
$$

an equivalent formulation of the solution.



# Problem 6.4
Given the SDE $dX_t = \lambda (\mu -X_t) dt + \sigma dW_t$  we get coefficients $a_1=-\lambda, a_2 = \lambda \mu$ matching exercise 5.2.

Let $Y_t = \phi(X_t, t) = X_t - m_t$. 
Thus as $\phi_t = d m_t, \phi_x = 1, \phi_{xx} = 0$ we get $dY_t = dX_t - dm_t$.

Thus as $dm_t = \lambda(\mu - m_t)dt$  from 5.2, it must be that

$$
dY_t = \lambda (\mu -X_t) dt + \sigma dW_t -  \lambda(\mu - m_t)dt \\
dY_t = -\lambda (X_t-m_t) dt + \sigma dW_t \\
dY_t = -\lambda (Y_t) dt + \sigma dW_t \\
$$
This is known to have the solution

$$
Y_t = \sigma e^{-\lambda t} \int_0^t e^{\lambda s} dW_s
$$

By transformation then,
$$
X_t =\sigma e^{-\lambda t} \int_0^t e^{\lambda s} dW_s + m_t \\
X_t =\sigma e^{-\lambda t} \int_0^t e^{\lambda s} dW_s + \mu(1-e^{-\lambda t}) + m_0e^{-\lambda t}\\
$$
As $m_0 = x_0$ this result is equivalent to the result in the book.


# Problem 6.5

Consider the SDE

$dS_t = dW_t$, i.e. $S_t = W_t$

Then for $Y_t = \phi(x,t) = \frac{1}{1-x}$ by the ito formula, we get that

$dY_t = 0 dt + \phi_x(S_t,t) dS_t + \frac{1}{2} 1^2 \phi_{xx}(S_t,t) dt$

Thus as $\phi_x = (1-x)^{-2}$ and $\phi_{xx}= 2(1-x)^{-3}$ we get 


$dY_t = (1-S_t)^{-2}dS_t + \frac{2}{2}(1-S_t)^{-3} dt = Y_t^2 dW_t + Y_t^3 dt$

Thus $Y_t = \frac{1}{1-W_t}$ solves the differential equation above. 



# Problem 7.1
## Expectation
Prove
$E\left[  \sum q(W_{t_i})(W_{t_{i+0.5}} - W_{t_i}) ^2  -0.5  \sum q(W_{t_i})   \delta t\right] = 0$

Note that due to independence of $q(W_{t_i})$ and $(W_{t_{i+0.5}} - W_{t_i})$

$\sum E\left[ q(W_{t_i})\right]  E\left[ (W_{t_{i+0.5}} - W_{t_i}) ^2\right]  - 0.5  \delta t E\left[ \sum q(W_{t_i})  \right]$

And note that 

$$
(W_{t_{i+0.5}} - W_{t_i}) \sim N(0,\frac{1}{\sqrt 2} \delta t)
$$

thus

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


I realized I stopped this problem halfway through.

# Problem 7.7 

Given an Ito SDE 
$$
dX_t = rX_t dt + \gamma X_t dW_t
$$
and corresponding Stratonovich SDE
$$
d\hat X_t = r \hat X_t dt + \gamma \hat X_t \circ dW_t
$$
we note that by equation 7.9, there is an ito SDE equivalent to the stratonovich SDE
written as
$$
d\bar X_t = \left( r - \frac{\gamma^2}{2}\right) \bar X_t dt + \gamma \bar X_t dW_t
$$
Giving the key relationship that for any function $F(\cdot)$, it mus t be that $F(\hat x) = F(\bar x)$

From the discussion surrounding 5.6, we know that 
$$
E(X_t) = E(X_0) e^{rt}
$$
Thus from the discussion above
$$
E(\hat X_t) = E(\bar X_t) = E(\hat X_0) e^{t(r-\frac{\gamma^2}{2})}
$$

Similarly, as Exercise 6.3 gives us that
$$
E(\log X_t) = \log X_0 + t(r-\frac{\gamma^2}{2})
$$
Thus by subsitution
$$
E(\log \hat X_t) = E(\log \bar X_t) = \log \hat X_0 
    + {t(r-\frac{\gamma^2}{2})} - t\frac{\gamma^2}{2} \\
E(\log \hat X_t) = E(\log \bar X_t) = \log \hat X_0 + tr
$$

# Problem PC 7.1 and Problem PC 7.2

See Submitted julia notebooks, html results, and pdf. 
They all say the same thing, but are different formats.
I suggest reading the HTML.

