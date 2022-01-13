# Homework 7
Background is 
## Question 1: Weak Convergence
My understanding of the question is that I should "prove" that EM converges in the weak sense. 
I don't recall if these were analytical or computational questions.
### OU

### GBM

### Brownian Bridge

## Question 2 

### P10.1

iter
###  P10.2

## Question 3

###  11.1

###  11.2

## Question 4

Find condition under which the Milstein Method is Mean Square Stable for Geometric Brownian Motion.
[Background on stability](20220107214350.md)

As the milstein method is summarized by the recursive definition
$$
Y_{n+1} 
    = Y_n + \mu(Y_n,t_n)\Delta t + \sigma(Y_n,t_n)\Delta W_t
    + \frac{1}{2} \sigma(Y_n,t_n) \sigma_{Y}(Y_n,t_n) (\Delta W_n^2 - \Delta t) \\
    = Y_n + \mu Y_n \Delta t + \sigma Y_n \Delta W_t
    + \frac{1}{2} \sigma Y_n \sigma (\Delta W_n^2 - \Delta t) \\
$$
Thus the condition for convergence is that $\lim E((Y_{n+1})^2) = 0$. 
Thus
$$
\left(
    Y_n + \mu Y_n \Delta t + \sigma Y_n \Delta W_t
    + \frac{1}{2} \sigma Y_n \sigma (\Delta W_n^2 - \Delta t)
\right)
\left(
   Y_n + \mu Y_n \Delta t + \sigma Y_n \Delta W_t
    + \frac{1}{2} \sigma Y_n \sigma (\Delta W_n^2 - \Delta t)
\right)
$$git s
$$
Y_n^2 \left(
    1 + \mu  \Delta t + \sigma  \Delta W_t
    + \frac{1}{2}  \sigma^2 (\Delta W_n^2 - \Delta t)
\right)^2 \\
Y_n^2 \left(
    1 + \mu^2  \Delta t^2 + \sigma^2  \Delta W_t^2
    + \frac{1}{4}  \sigma^4 (\Delta W_n^2 - \Delta t)^2
    + 2\left[
        1 
        + \mu  \Delta t 
        + \sigma  \Delta W_t
        + \frac{1}{2}  \sigma^2 (\Delta W_n^2 - \Delta t)
        + \mu  \Delta t \sigma  \Delta W_t
        + \frac{1}{2} \mu  \Delta t \sigma^2 (\Delta W_n^2 - \Delta t)
        +\sigma  \Delta W_t \frac{1}{2}  \sigma^2 (\Delta W_n^2 - \Delta t)
    \right]
\right)
\\
E(Y_n^2) \left(
    1 + \mu^2  \Delta t^2 + \sigma^2  E(\Delta W_t^2)
    + \frac{1}{4}  \sigma^4 E((\Delta W_n^2 - \Delta t)^2)
    + 2\left[
        1 
        + \mu  \Delta t 
        + \sigma  E(\Delta W_t)
        + \frac{1}{2}  \sigma^2 (E(\Delta W_n^2) - \Delta t)
        + \mu  \Delta t \sigma  E(\Delta W_t)
        + \frac{1}{2} \mu  \Delta t \sigma^2 (E(\Delta W_n^2) - \Delta t)
        +\sigma   \frac{1}{2}  \sigma^2 (E(\Delta W_n^3) - E(\Delta W_t)\Delta t)
    \right]
\right)
\\
E(Y_n^2) \left(
    1 + \mu^2  \Delta t^2 + \sigma^2  \Delta t
    + \frac{1}{4}  \sigma^4 E((\Delta W_n^2 - \Delta t)^2)
    + 2\left[
        1 
        + \mu  \Delta t 
        + 0
        + \frac{1}{2}  \sigma^2 (\Delta t - \Delta t)
        + \mu  \Delta t \sigma 0
        + \frac{1}{2} \mu  \Delta t \sigma^2 (\Delta t - \Delta t)
        +\sigma   \frac{1}{2}  \sigma^2 (0 - 0\Delta t)
    \right]
\right)
\\
$$
Recall that $E(\Delta W_n) = 0$, $E(\Delta W_n^2) = \Delta t$, $E(\Delta W_n^3) = 0$, $E(\Delta W_n^4) = 0$ . Simplifying gives:

$$
E(Y_{n+1}^2) = E(Y_n^2) \left(
    1 + \mu^2  \Delta t^2 + \sigma^2  \Delta t
    + \frac{1}{4}  \sigma^4 E((\Delta W_n^2 - \Delta t)^2)
    + 2 \left[ 1 + \mu  \Delta t \right]
\right)
\\
$$
Note that
$$
E((\Delta W_n^2 - \Delta t)^2) = E(\Delta W_n^4 -2 \Delta W_n^2 \Delta t + \Delta t^2) = \Delta t^2
$$
Thus
$$
E(Y_{n+1}^2) = E(Y_n^2) \left(
    1 + \mu^2  \Delta t^2 + \sigma^2  \Delta t
    + \frac{1}{4}  \sigma^4 \Delta t^2
    + 2 \left[ 1 + \mu  \Delta t \right]
\right)
\\
$$
Thus the limit condition for mean square stability of the milstein method for GMB is:
$$
1 + \mu^2  \Delta t^2 + \sigma^2  \Delta t
    + \frac{1}{4}  \sigma^4 \Delta t^2
    + 2 \left[ 1 + \mu  \Delta t \right]
< 0
\\
1 + \mu^2  \Delta t^2 + \sigma^2  \Delta t
    + \frac{1}{4}  \sigma^4 \Delta t^2
    + 2 \left[ 1 + \mu  \Delta t \right]
< 0
$$

## Question 5


### P11.1

### P11.2