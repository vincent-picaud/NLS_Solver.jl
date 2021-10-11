# Kunisch-Rendl solver

## Introduction

This file implements the [Kunisch-Rendl method](https://epubs.siam.org/doi/10.1137/S1052623400376135). This methods is a primal-dual active set method for quadratic problems with bound constraints.

**Attention:** this method is not always convergent. But:
- in our context this is not a problem as we can enforce strong
  diagonal dominance of the matrix and the method is convergent.
- if required there is a modification of the algorithm,
  [Hungerlaender-Rendl
  paper](http://www.optimization-online.org/DB_HTML/2016/09/5644.html),
  that insures convergence.


## Algorithm

The goal is to minimize a bound constrained quadratic form where ``Q``
is a **symmetric definite positive** matrix.

```math
\min\limits_{a\le x \le b} \frac{1}{2} x^t.Q.x + q^t.x
```

For well conditionned problems the method is quite effective. This
method is also intuitive and easy to implement.

The idea of the method is to guess active constraints, computes the
associated solution with its **a posteriori multipliers**. Then we
check the validity of our initial guess and we make the necessary
corrections and we restart with this new guess. It can be shown, in
exact arithmetic, that this process ends in a finite number of
iterations.

To do this we need to write the Lagragian associated to problem

```math
\mathcal{L} = \frac{1}{2} x^t.Q.x + q^t.x + \lambda^t (a-x) + \mu^t (x-b)
```

and its gradient

```math
\nabla_x \mathcal{L} = Q.x + q -\lambda + \mu = Q.x+q+\tau = 0
```
with ``\tau =  -\lambda + \mu`` storing multipliers associated to the constraints ``x\in[a,b]``

The KKT conditions are
```math
\nabla_x \mathcal{L} = Q.x+q-\lambda+\mu =0
```
```math
(\lambda \ge 0)\wedge(\lambda\odot(a-x))=0 
```
```math
(\mu \ge 0)\wedge(\mu\odot(x-b))=0
```

We introduce a ``Z\in\{-1,0,1\}^n`` vector to encode our guess on
which contraints are active or not:

```math
Z_i = \left\{
\begin{array}{c|c|cc}
-1 & x_i=a_i & \lambda_i = -\tau_i &  \mu_i = 0 \\
0  & a_i \le x_i \le b_i & \lambda_i = 0 & \mu_i = 0 \\
+1 &  x_i=b_i & \lambda_i = 0 &  \mu_i = \tau_i \\
\end{array}
\right.
``` 

### Inputs
 
``Z^{(0)}\in\{-1,01\}^n`` initial guess 
``x^{(0)}`` initial guess

### Result 

``x^{(k)}`` and ``\tau^{(k)}`` solution of the problem. ``\tau``
stores Lagrange multipliers in a compact way. ``\tau[i]>0`` means that
the upper bound constraint is active, the associated multiplier is
``\tau[i]``. ``\tau[i]<0`` means that the lower bound constraint is
active, the associated multiplier is ``-\tau[i]``. ``\tau[i]=0`` means
that bounds constraints are inactive.
	
### Loop

Repeat until ``Z^{(k)}\neq Z^{(k+1)}``

#### Modifies ``Q`` and ``q`` such that ``\tilde{x}^{(k)}`` fulfills active constraints

```math
\tilde{x}^{(k)}= \arg\min\limits_{x}\frac{1}{2} x^t.Q.x + q^t.x
```

```math
\begin{array}{lll}
x_i=a_i & \text{if} & Z_i=-1 \\ 
x_i=b_i & \text{if} & Z_i=+1
\end{array}
```

#### Compute Multiplicateurs a posteriori

```math
\tau^{(k)} = -( Q.\tilde{x}^{(k)}+q )
```

#### Update ``Z`` (corrections -> our new guess)
```math
\begin{array}{lrc}
  
  \text{If\ }Z_i^{(k)}=-1 & \text{then} & Z_i^{(k+1)}=
  \left\{
  \begin{array}{rl}
    -1 & \text{if\ }\tau_i\le 0 \\
    0 & \text{otherwise}
  \end{array}
  \right.  \\

  \text{If\ }Z_i^{(k)}=0 & \text{then} & Z_i^{(k+1)}=
  \left\{
  \begin{array}{rl}
    -1 & \text{if\ }x_i<a_i \\
    0 & \text{if\ }a_i\le x_i  \le b_i \\
    +1 & \text{if\ }b_i<x_i
  \end{array}
  \right.  \\

  \text{If\ }Z_i^{(k)}=+1 & \text{then} & Z_i^{(k+1)}=
  \left\{
  \begin{array}{rl}
    1 & \text{if\ }\tau_i\ge 0 \\
    0 & \text{otherwise}
  \end{array}
  \right.  \\
  
\end{array}
```

## Documentation

```@autodocs
Modules = [NLS_Solver]
Pages = ["QuadSolvers/Kunisch-Rendl.jl"]
```
