```@meta
CurrentModule = NLS_Solver
```

```@setup session
using NLS_Solver
```

# Nonlinear regressions

We have seen that this package solves problems of this form:
```math
\min\limits_{\theta} \frac{1}{2}\|r(\theta)\|^2
```

Classical fitting problems fall under this category. If you have data
constituted by two vectors ``X, Y`` of ``\mathbb{R}^{n_S}``, then you
can define a model

```math
\begin{align*}
m:& \mathbb{R}\times\mathbb{R}^{n_θ}\to\mathbb{R} \\
  & (x,θ) \mapsto y = m(x,θ)
\end{align*}
```
and a residue vector ``r`` of ``\mathbb{R}^{n_S}`` of components ``r_i(θ)=Y_i - m(X_i,θ)``.

## To be continued 

dsfdqsfqdsf
