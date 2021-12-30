var documenterSearchIndex = {"docs":
[{"location":"todo/#TODO-List","page":"TODO List","title":"TODO List","text":"","category":"section"},{"location":"todo/#Factorize-ResultTypes","page":"TODO List","title":"Factorize ResultTypes","text":"","category":"section"},{"location":"todo/","page":"TODO List","title":"TODO List","text":"Our Abstract_Solver_Result, Abstract_BC_Solver_Result, AbstractQuadResult result structures all have some common methods like converged(), iteration_count() or solution()","category":"page"},{"location":"todo/","page":"TODO List","title":"TODO List","text":"-> one must define and use a Base type with these common functions.","category":"page"},{"location":"api/#API-index","page":"API index","title":"API index","text":"","category":"section"},{"location":"api/","page":"API index","title":"API index","text":"The global API index is as follows:","category":"page"},{"location":"api/","page":"API index","title":"API index","text":"Modules = [NLS_Solver]","category":"page"},{"location":"api/#Public","page":"API index","title":"Public","text":"","category":"section"},{"location":"api/","page":"API index","title":"API index","text":"Modules = [NLS_Solver]\nPrivate = false\t","category":"page"},{"location":"api/#NLS_Solver.AbstractNLS","page":"API index","title":"NLS_Solver.AbstractNLS","text":"abstract type AbstractNLS end \n\nDefines an abstract non-linear least squares problem (NLS). In our context such problem is essentially a differentiable function mathbfr:\n\nmathbfr thetainmathbbR^n_θmapsto mathbfr(mathbftheta)inmathbbR^n_S\n\nwhere:\n\nmathbfr(mathbfθ)mathbbR^n_S is the residue vector,\nmathbfθmathbbR^n_θ is the parameter vector to be optimized\n\nThe objective function to minimize is:\n\nf(θ)=frac12 mathbfr(θ) ^2\n\nThe classical approach uses a linear approximation of mathbfr:\n\nmathbfr(mathbfθ+δmathbfθ)approx mathbfr(mathbfθ) + mathbfJ(mathbfθ)cdot δmathbfθ\n\nwhere mathbfJ is the Jacobian:\n\nmathbfJ_ij=partial_j r^i(mathbfθ) iin1n_S jin1n_θ\n\nThis leads to\n\nf(mathbfθ+δmathbfθ)approx f(mathbfθ) + langle nabla f δmathbfθ rangle + frac12  langle nabla^2 f cdot δmathbfθ  δmathbfθ rangle\n\nWhere the gradient nabla f is mathbfJ^t mathbfr and the (approximate) Hessian nabla^2 f is mathbfJ^t mathbfJ.\n\nTo implement such model, you must define the following functions:\n\nparameter_size : returns n_θ\nresidue_size : returns n_S\neval_r : compute mathbfr\neval_r_J : compute (mathbfr mathbfJ)\n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.Abstract_BC_Solver_Conf","page":"API index","title":"NLS_Solver.Abstract_BC_Solver_Conf","text":"abstract type Abstract_BC_Solver_Conf end\n\nAbstract solver configuration. These are the solvers to be used to solve bound constrained nonlinear least squares:\n\nminlimits_theta_l le theta le theta_u  frac12r(theta)^2\n\nImplementations:\n\nLevenbergMarquardt_BC_Conf \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.Abstract_Solver_Conf","page":"API index","title":"NLS_Solver.Abstract_Solver_Conf","text":"abstract type Abstract_Solver_Conf end\n\nAbstract solver configuration. These are the solvers to be used to solve unconstrained nonlinear least squares:\n\nminlimits_theta frac12r(theta)^2\n\nImplementations:\n\nLevenbergMarquardt_Conf \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.BoundConstraints","page":"API index","title":"NLS_Solver.BoundConstraints","text":"Store bound constraints lu\n\nPresence of NaN component and the lle u condition is checked at construction time. Note however that some components can be infinite.\n\nThe following constructors are available:\n\nConstruct 0010^n \n\nBoundConstraints(n)\n\nConstruct T(0)T(1)^n where components are of type T\n\nBoundConstraints(T,n)\n\nConstruct lu where l and u are lower and upper bound vectors\n\nBoundConstraints(l,u)\n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.LevenbergMarquardt_BC_Conf","page":"API index","title":"NLS_Solver.LevenbergMarquardt_BC_Conf","text":"mutable struct LevenbergMarquardt_BC_Conf <: Abstract_Solver_Conf\n    ...\nend\n\nUse this constructor\n\nLevenbergMarquardt_BC_Conf()\n\nto initialize the bound constrained Levenberg-Marquardt solver default configuration parameters.\n\nTo solve a problem with this method, you must then call  solve(nls::AbstractNLS, θ_init::AbstractVector, bc::BoundConstraints, conf::Abstract_BC_Solver_Conf) \n\nSee: \n\nset_max_iteration!(conf::LevenbergMarquardt_BC_Conf,max_iter::Int) \nset_ε_grad_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,ε_grad_Inf_norm::Float64) \nset_ε_step_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,ε_step_Inf_norm::Float64) \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.LevenbergMarquardt_Conf","page":"API index","title":"NLS_Solver.LevenbergMarquardt_Conf","text":"mutable struct LevenbergMarquardt_Conf <: Abstract_Solver_Conf\n    ...\nend\n\nUse this constructor\n\nLevenbergMarquardt_Conf()\n\nto initialize the Levenberg-Marquardt solver default configuration parameters.\n\nTo solve a problem with this method, you must then call  solve(nls::AbstractNLS, θ_init::AbstractVector, conf::Abstract_Solver_Conf) \n\nSee: \n\nset_max_iteration!(conf::LevenbergMarquardt_Conf,max_iter::Int) \nset_ε_grad_Inf_norm!(conf::LevenbergMarquardt_Conf,ε_grad_Inf_norm::Float64) \nset_ε_step_Inf_norm!(conf::LevenbergMarquardt_Conf,ε_step_Inf_norm::Float64) \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.NLS_ForwardDiff","page":"API index","title":"NLS_Solver.NLS_ForwardDiff","text":"struct NLS_ForwardDiff <: AbstractNLS\n    ...\nend\n\n\nA specialization that uses the ForwardDiff package to compute the Jacobian.\n\nBy comparison with AbstractNLS you only have to define these functions:\n\nparameter_size : returns n_θ\nresidue_size : returns n_S\neval_r : computation of mathbfr\n\nSee: create_NLS_problem_using_ForwardDiff \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.converged-Tuple{NLS_Solver.Abstract_BC_QuadSolver_Result}","page":"API index","title":"NLS_Solver.converged","text":"converged(::Abstract_BC_QuadSolver_Result)\n\nReturn true if the solver converged\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.converged-Tuple{NLS_Solver.Abstract_Solver_Result}","page":"API index","title":"NLS_Solver.converged","text":"converged(::Abstract_Solver_Result)\n\nReturn true if the solver converged\n\nSee: Abstract_Solver_Result \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.create_NLS_problem_using_ForwardDiff-Tuple{Function, Pair{Int64, Int64}}","page":"API index","title":"NLS_Solver.create_NLS_problem_using_ForwardDiff","text":"create_NLS_problem_using_ForwardDiff(r::Function;domain_image_dim::Pair{Int,Int})\n\nCreates an AbstractNLS specialized instance where the eval_r_J function is automatically defined using automatic differentiation.\n\nr is a function that maps a parameter vector θ to its residue. The Jacobian matrix is computed using the ForwardDiff package.\ndomain_image_dim is a pair of the form θ length => r length that defines domain and codomain dimensions.\n\nUsage example\n\nAn example defining the Rosenbrock function\n\nfrac12r(theta)^2text where r = sqrt2 left( beginarrayc  1-theta_1  10(theta_2-theta_1^2) endarray right)\n\nnls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ\n     sqrt(2) sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]\nend\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.eval_nls_fobj-Tuple{AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.eval_nls_fobj","text":"eval_nls_fobj(r::AbstractVector{T}) -> f(θ)\n\nCompute f(θ)=frac12 mathbfr(mathbfθ) ^2\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.eval_nls_∇fobj-Tuple{AbstractVector{T} where T, AbstractMatrix{T} where T}","page":"API index","title":"NLS_Solver.eval_nls_∇fobj","text":"eval_nls_∇fobj(r,J) -> ∇fobj\n\nCompute the gradient: nabla f(mathbfθ) = mathbfJ^tmathbfr\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.eval_nls_∇∇fobj-Tuple{AbstractMatrix{T} where T}","page":"API index","title":"NLS_Solver.eval_nls_∇∇fobj","text":"eval_nls_∇∇fobj(J) -> ∇∇fobj\n\nCompute the (approximate) Hessian: nabla^2 f(mathbfθ) = mathbfJ^tmathbfJ\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.eval_r-Tuple{AbstractNLS, AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.eval_r","text":"eval_r(nls::AbstractNLS,\n        θ::AbstractVector) -> r\n\nCompte the residual vector mathbfr\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.eval_r_J-Tuple{AbstractNLS, AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.eval_r_J","text":"eval_r_J(nls::AbstractNLS,θ::AbstractVector) -> (r,J)\n\nCompute the residual the vector mathbfr and its Jacobian mathbfJ \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.iteration_count-Tuple{NLS_Solver.Abstract_BC_QuadSolver_Result}","page":"API index","title":"NLS_Solver.iteration_count","text":"iteration_count(::Abstract_BC_QuadSolver_Result)\n\nReturn the number of consumed iteration\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.iteration_count-Tuple{NLS_Solver.Abstract_Solver_Result}","page":"API index","title":"NLS_Solver.iteration_count","text":"iteration_count(::Abstract_Solver_Result)\n\nReturn the number of consumed iteration\n\nSee: Abstract_Solver_Result \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.lower_bound-Tuple{BoundConstraints}","page":"API index","title":"NLS_Solver.lower_bound","text":"lower_bound(bc::BoundConstraints)\n\nReturn lower bound l\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.objective_value-Tuple{NLS_Solver.Abstract_BC_QuadSolver_Result}","page":"API index","title":"NLS_Solver.objective_value","text":"objective_value(::Abstract_BC_QuadSolver_Result)\n\nReturns objective value at the point solution.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.objective_value-Tuple{NLS_Solver.Abstract_Solver_Result}","page":"API index","title":"NLS_Solver.objective_value","text":"objective_value(::Abstract_Solver_Result)\n\nReturns objective value at the point solution.\n\nSee: Abstract_Solver_Result \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.parameter_size-Tuple{AbstractNLS}","page":"API index","title":"NLS_Solver.parameter_size","text":"parameter_size(nls::AbstractNLS)\n\nReturn the dimension n_θ of the parameter vector θ.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.project!-Union{Tuple{N}, Tuple{AbstractArray{var\"#s1\", N} where var\"#s1\"<:Real, BoundConstraints{var\"#s9\", N, LBT, UBT} where {var\"#s9\"<:Real, LBT<:AbstractArray{var\"#s9\", N}, UBT<:AbstractArray{var\"#s9\", N}}}} where N","page":"API index","title":"NLS_Solver.project!","text":"project!(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N})\n\nProject x such that x in lu is fullfiled.\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.residue_size-Tuple{AbstractNLS}","page":"API index","title":"NLS_Solver.residue_size","text":"sample_size(nls::AbstractNLS)\n\nReturn the dimension n_S of the residue vector r.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.set_max_iteration!-Tuple{LevenbergMarquardt_BC_Conf, Int64}","page":"API index","title":"NLS_Solver.set_max_iteration!","text":"set_max_iteration!(conf::LevenbergMarquardt_BC_Conf,\n                   max_iter::Int)\n\nModify the maximum number of iterations\n\nSee: LevenbergMarquardt_BC_Conf \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.set_max_iteration!-Tuple{LevenbergMarquardt_Conf, Int64}","page":"API index","title":"NLS_Solver.set_max_iteration!","text":"set_max_iteration!(conf::LevenbergMarquardt_Conf,\n                   max_iter::Int)\n\nModify the maximum number of iterations\n\nSee: LevenbergMarquardt_Conf \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.set_ε_grad_Inf_norm!-Tuple{LevenbergMarquardt_BC_Conf, Float64}","page":"API index","title":"NLS_Solver.set_ε_grad_Inf_norm!","text":"set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,\n                     ε_grad_Inf_norm::Float64)\n\nModify the stopping criterion nabla f_inftyleepsilon\n\nSee: LevenbergMarquardt_BC_Conf \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.set_ε_grad_Inf_norm!-Tuple{LevenbergMarquardt_Conf, Float64}","page":"API index","title":"NLS_Solver.set_ε_grad_Inf_norm!","text":"set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_Conf,\n                     ε_grad_Inf_norm::Float64)\n\nModify the stopping criterion nabla f_inftyleepsilon\n\nSee: LevenbergMarquardt_Conf \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.set_ε_step_Inf_norm!-Tuple{LevenbergMarquardt_BC_Conf, Float64}","page":"API index","title":"NLS_Solver.set_ε_step_Inf_norm!","text":"set_ε_step_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,\n                     ε_step_Inf_norm::Float64)\n\nModify the stopping criterion delta x_inftyleepsilon\n\nSee: LevenbergMarquardt_BC_Conf \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.set_ε_step_Inf_norm!-Tuple{LevenbergMarquardt_Conf, Float64}","page":"API index","title":"NLS_Solver.set_ε_step_Inf_norm!","text":"set_ε_step_Inf_norm!(conf::LevenbergMarquardt_Conf,\n                     ε_step_Inf_norm::Float64)\n\nModify the stopping criterion delta x_inftyleepsilon\n\nSee: LevenbergMarquardt_Conf \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.solution-Tuple{NLS_Solver.Abstract_BC_QuadSolver_Result}","page":"API index","title":"NLS_Solver.solution","text":"solution(::Abstract_BC_QuadSolver_Result)\n\nReturns the founded solution \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.solution-Tuple{NLS_Solver.Abstract_Solver_Result}","page":"API index","title":"NLS_Solver.solution","text":"solution(::Abstract_Solver_Result)\n\nReturns the founded solution \n\nSee: Abstract_Solver_Result \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.solve-Tuple{AbstractNLS, AbstractVector{T} where T, Abstract_Solver_Conf}","page":"API index","title":"NLS_Solver.solve","text":"solve(nls::AbstractNLS,\n      θ_init::AbstractVector,\n      conf::Abstract_Solver_Conf)::Abstract_Solver_Result\n\nGeneric interface to solve an AbstractNLS problem.\n\nThe used algorithm is defined through Abstract_Solver_Conf specializations.\n\nThe method returns a Abstract_Solver_Result specialization.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.solve-Tuple{AbstractNLS, AbstractVector{T} where T, BoundConstraints, Abstract_BC_Solver_Conf}","page":"API index","title":"NLS_Solver.solve","text":"solve(nls::AbstractNLS,\n      θ_init::AbstractVector,\n      bc::BoundConstraints,\n      conf::Abstract_BC_Solver_Conf) -> Abstract_Solver_Result\n\nGeneric interface to solve a AbstractNLS problem with bound constraints.\n\nThe used algorithm is defined through Abstract_BC_Solver_Conf specializations.\n\nThe method returns a Abstract_BC_Solver_Result specialization.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.upper_bound-Tuple{BoundConstraints}","page":"API index","title":"NLS_Solver.upper_bound","text":"upper_bound(bc::BoundConstraints)\n\nReturn upper bound u\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Private","page":"API index","title":"Private","text":"","category":"section"},{"location":"api/","page":"API index","title":"API index","text":"Modules = [NLS_Solver]\nPublic = false\t","category":"page"},{"location":"api/#NLS_Solver.Abstract_BC_QuadSolver_Result","page":"API index","title":"NLS_Solver.Abstract_BC_QuadSolver_Result","text":"abstract type Abstract_BC_QuadSolver_Result end\n\nQuadratic solver result abstraction\n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.Abstract_BC_Solver_Result","page":"API index","title":"NLS_Solver.Abstract_BC_Solver_Result","text":"The structure returned by solve when using the LevenbergMarquardt_BC_Conf method.\n\nSee Abstract_Solver_Result \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.Abstract_Solver_Result","page":"API index","title":"NLS_Solver.Abstract_Solver_Result","text":"abstract type Abstract_Solver_Result end\n\nThis is the base type returned by the solve method. It contains the information related to the founded solution.\n\nInterface\n\nconverged \niteration_count \nobjective_value \nsolution \n\nImplementations\n\nLevenbergMarquardt_Result \nLevenbergMarquardt_BC_Result \n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.BoundConstraint_Enum","page":"API index","title":"NLS_Solver.BoundConstraint_Enum","text":"@enum(BoundConstraint_Enum,\n      ACTIVE_LB = -1,\n      INACTIVE_BC = 0,\n      ACTIVE_UB = +1)\n\nAn enum to store bound constraint state.\n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.LevenbergMarquardt_BC_Result","page":"API index","title":"NLS_Solver.LevenbergMarquardt_BC_Result","text":"struct LevenbergMarquardt_BC_Result{T<:Real} <:  Abstract_BC_Solver_Result\n    ...\nend \n\nThis structure subtypes Abstract_BC_Solver_Result\n\n\n\n\n\n","category":"type"},{"location":"api/#NLS_Solver.LevenbergMarquardt_Result","page":"API index","title":"NLS_Solver.LevenbergMarquardt_Result","text":"struct LevenbergMarquardt_Result{T<:Real} <: Abstract_Solver_Result\n   ...\nend\n\nThis structure subtypes Abstract_Solver_Result\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.:+-Tuple{BoundConstraints, AbstractArray}","page":"API index","title":"Base.:+","text":"Translate bound constraints\n\na+τb+τ = ab+τ\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Base.:--Tuple{BoundConstraints, AbstractArray}","page":"API index","title":"Base.:-","text":"Translate bound constraints\n\na-τb-τ = ab-τ\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Base.axes-Tuple{BoundConstraints}","page":"API index","title":"Base.axes","text":"axes(bc::BoundConstraints)\n\nReturn bound axes\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Base.eltype-Union{Tuple{BoundConstraints{ELT, N, LBT, UBT} where {N, LBT<:AbstractArray{ELT, N}, UBT<:AbstractArray{ELT, N}}}, Tuple{ELT}} where ELT","page":"API index","title":"Base.eltype","text":"eltype(bc::BoundConstraints)\n\nReturn bound element type\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Base.in-Union{Tuple{N}, Tuple{AbstractArray{var\"#s6\", N} where var\"#s6\"<:Real, BoundConstraints{var\"#s7\", N, LBT, UBT} where {var\"#s7\"<:Real, LBT<:AbstractArray{var\"#s7\", N}, UBT<:AbstractArray{var\"#s7\", N}}}} where N","page":"API index","title":"Base.in","text":"in(bc::BoundConstraints)\n\nCheck if xin lu\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Base.length-Tuple{BoundConstraints}","page":"API index","title":"Base.length","text":"length(bc::BoundConstraints)\n\nReturn bound length\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#Base.size-Tuple{BoundConstraints}","page":"API index","title":"Base.size","text":"size(bc::BoundConstraints)\n\nReturn bound size\n\nSee: BoundConstraints \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.check_first_order-Tuple{AbstractVector{var\"#s5\"} where var\"#s5\"<:Real, AbstractVector{var\"#s3\"} where var\"#s3\"<:Real, BoundConstraints{var\"#s2\", 1, LBT, UBT} where {var\"#s2\"<:Real, LBT<:AbstractVector{var\"#s2\"}, UBT<:AbstractVector{var\"#s2\"}}}","page":"API index","title":"NLS_Solver.check_first_order","text":"check_first_order(∇f::AbstractVector{<:Real},\n                  xstar::AbstractVector{<:Real},\n                  bc::BoundConstraints{<:Real,1})\n\ncheck_first_order(Q::Symmetric{<:Real},\n                  q::AbstractVector{<:Real},\n                  xstar::AbstractVector{<:Real},\n                  bc::BoundConstraints{<:Real,1})\n\nCheck First-Order Conditions  (see Bound Constrained Optimization slides)\n\nIf x^star=argmin f(x) xinlu then:\n\npartial_i f(x^star) = leftbeginarrayll\nge 0  textif  x^stari = li \n= 0  textif  li le x^stari le ui \nle 0  textif  x^stari = ui \nendarray\nright\n\nThis is equivalent to:\n\nx^star = P_lu(x^star-nabla f(x^star))\n\nAccording to the previous result, this function returns:\n\nmax mid x^star - P_lu(x^star-(Qx^star+q)) mid\n\nFor a local stationary point this quantity must be null \n\nThe second function is a wrapper that computes f=Qx^star+q\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.clean_τ!-Tuple{AbstractArray{var\"#s9\", N} where {var\"#s9\"<:Real, N}, AbstractArray{NLS_Solver.BoundConstraint_Enum, N} where N}","page":"API index","title":"NLS_Solver.clean_τ!","text":"clean_τ!(τ::AbstractArray{<:Real},              \n         Z::AbstractArray{BoundConstraint_Enum})\n\nBy definition τ=-Qx-q. If the algorithm converged, then one must have τ[i]=0 when the constraint is inactive.\n\nThis function updates τ by overwriting τ[i]=0 when Z[i]=inactive. \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.compute_δL_constrained-Tuple{AbstractVector{T} where T, Real, AbstractVector{T} where T, AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.compute_δL_constrained","text":"Same idea than compute_δL_unconstrained, however when bound constraints are present h is such that:\n\n(nabla^2 f + mu I)h + nabla f + tau = 0\n\nit follows that:\n\nδL = L(0)-L(h) = frac12 langle h mu h + tau - nabla f rangle\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.compute_δL_unconstrained-Tuple{AbstractVector{T} where T, Real, AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.compute_δL_unconstrained","text":"Compute δL = L(0)-L(h) where L is the quadratic model\n\nL(h)=f(theta) + langle nabla f h rangle + frac12langle nabla^2 f h h rangle\n\nwith f(theta)=frac12 r(θ) _2^2, nabla f = J^t r and nabla^2 f = J^t J\n\nA direct computation gives:\n\nδL = L(0)-L(h) = -left(  langle J^tr h rangle + frac12langle nabla^2 f h h rangle right)\n\nHowever one can avoid the computation of nabla^2 f h if one uses the fact that h is solution of:\n\n(nabla^2 f + mu I)h + nabla f = 0\n\nWith this hypothesis, one gets:\n\nδL = L(0)-L(h) = frac12 langle h mu h - nabla f rangle\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.compute_δf-Tuple{AbstractVector{T} where T, AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.compute_δf","text":"Compute true variation of the real model: δf = frac12(r^t(θ)r(θ)-r^t(θ+h)r(θ+h))\n\nContrary to δL things are simpler. However a trick is to use an equivalent formulation:\n\nδf = frac12(r^t(θ)r(θ)-r^t(θ+h)r(θ+h)) = frac12(r(θ)-r(θ+h))^t(r(θ)+r(θ+h))\n\nthat has a better numerical behavior. \n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.initialize_x_Z-Tuple{AbstractArray, BoundConstraints}","page":"API index","title":"NLS_Solver.initialize_x_Z","text":"initialize_x_Z(x_init::AbstractArray,\n               bc::BoundConstraints)\n\nCreate (x,Z) from initial guess x_init and bound constraints bc\n\nZ is created by recording how x_init fulfils the bound constraints bc.\n\nx is the projection of x_init on the bounded domain [l,b].\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.multiplier_τ-Tuple{NLS_Solver.Abstract_BC_QuadSolver_Result}","page":"API index","title":"NLS_Solver.multiplier_τ","text":"multiplier_τ(::Abstract_BC_QuadSolver_Result)\n\nReturns the multipliers stored in a compact form (see τ definition, TODO)\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.quadratic_subproblem-Tuple{LinearAlgebra.Symmetric{var\"#s6\", S} where {var\"#s6\"<:Real, S<:(AbstractMatrix{var\"#s814\"} where var\"#s814\"<:var\"#s6\")}, AbstractVector{var\"#s5\"} where var\"#s5\"<:Real, AbstractVector{var\"#s3\"} where var\"#s3\"<:Real, BoundConstraints{var\"#s2\", 1, LBT, UBT} where {var\"#s2\"<:Real, LBT<:AbstractVector{var\"#s2\"}, UBT<:AbstractVector{var\"#s2\"}}, NLS_Solver.Abstract_BC_QuadSolver_Conf, NLS_Solver.LM_Damping, Int64}","page":"API index","title":"NLS_Solver.quadratic_subproblem","text":"Solve\n\nminlimits_hθ^l-θθ^u-θfrac12h^t(H+μI)h + f^th\n\nWe use the quadratic model of f, the bound contraints are such that the step h makes the update x+h falls in the θ^lθ^u bound.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.restrict_to_inactive!-Tuple{LinearAlgebra.Symmetric, AbstractVector{T} where T, AbstractVector{NLS_Solver.BoundConstraint_Enum}, AbstractVector{var\"#s3\"} where var\"#s3\"<:Real, AbstractVector{var\"#s2\"} where var\"#s2\"<:Real}","page":"API index","title":"NLS_Solver.restrict_to_inactive!","text":"restrict_to_inactive!(Q::Symmetric,\n                      q::AbstractVector,\n                      Z::AbstractVector{BoundConstraint_Enum},\n                      lb::AbstractVector{<:Real},\n                      ub::AbstractVector{<:Real})\n\nfunction restrict_to_inactive!(Q::Symmetric,\n                               q::AbstractVector,\n                               Z::AbstractVector{BoundConstraint_Enum},\n                               bc::BoundConstraints{<:Real,1})\n\nIn-place modification of (Qq) that produces (tildeQtildeq) such that the initial optimization problem:\n\ntildex = argmin frac12 x^t Q x + q^tx\n\nunder these constraints:\n\nxi = leftbeginarrayll \nli  textif  Zi = -1 \nui  textif  Zi = +1 \nendarrayright\n\nis transformed into this linear system:\n\ntildeQtildex+tildeq=0\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.update_Z!-Tuple{AbstractVector{T} where T, AbstractVector{T} where T, AbstractVector{NLS_Solver.BoundConstraint_Enum}, AbstractVector{T} where T, AbstractVector{T} where T}","page":"API index","title":"NLS_Solver.update_Z!","text":"update_Z!(x::AbstractVector,\n          τ::AbstractVector,\n          Z::AbstractVector{BoundConstraint_Enum},\n          lb::AbstractVector,\n          ub::AbstractVector)\n\nupdate_Z!(x::AbstractVector,\n          τ::AbstractVector,\n          Z::AbstractVector{BoundConstraint_Enum},\n          bc::BoundConstraints{<:Real,1})\n\nThis function updates Z according to x, τ and bounds lb, ub values.\n\nIt also count how many changes have be done during this update.\n\nNo change means that the algorithm has converged.\n\nNote: this function only modifies Z and return the number of bad hypothesis.\n\n\n\n\n\n","category":"method"},{"location":"api/#NLS_Solver.update_x!-Tuple{AbstractArray, AbstractArray{NLS_Solver.BoundConstraint_Enum, N} where N, BoundConstraints}","page":"API index","title":"NLS_Solver.update_x!","text":"update_x!(x::AbstractArray,\n          Z::AbstractArray{BoundConstraint_Enum},\n          bc::BoundConstraints)\n\nUpdate x value such that:\n\nxi = leftbeginarrayll \nli  textif  Zi = -1 \nui  textif  Zi = +1 \nendarrayright\n\nWhen Zi=0 the xi value is unaffected.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = NLS_Solver","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using NLS_Solver","category":"page"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This package goal is to solve nonlinear least squares problems. It currently supports two kind of problems:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"classical nonlinear least squares:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"minlimits_theta frac12r(theta)^2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"bound constrained nonlinear least squares","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"minlimits_theta_llethetaletheta_u frac12r(theta)^2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This package is easy to use. It follows a basic template where we have a solve() function of the form:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"solve(problem, algorithm_conf)::algorithm_result","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"problem is problem dependant data\nalgorithm_conf defines the chosen algorithm\nalgorithm_result is a structure that contains the result.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For the moment there are only two implemented methods. The classical Levenberg-Marquardt method for unconstrained problems. To use this method call LevenbergMarquardt_Conf. The other implemented method is a modification of the Levenberg-Marquardt where the inner quadratic problem is solved by the Kunisch-Rendl method to handle bound constraints. To use this method call LevenbergMarquardt_BC_Conf.  Please read this package README.org for extra details.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The two following examples illustrate how to solve a classical nonlinear least squares problem and a bound constrained one.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"In both cases, the objective function is the Rosenbrock function:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(1-θ_1)^2 + 100 (θ_2 - θ_1^2)^2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The classical problem is solved as follows:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"# define the objective methods\nnls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ\n  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]\nend\n\n# choose the method\nconf = LevenbergMarquardt_Conf()\n\n# initial value for θ\nθ_init = zeros(2)\n\n# solve it\nresult = solve(nls, θ_init, conf)\n\n# use the returned result\n@assert converged(result)\nθ_solution = solution(result)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Now to solve the same problem, but with bound constraints, proceed as follows:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"# choose the method\nconf = LevenbergMarquardt_BC_Conf()\n\n# define bound constraints\nθl = Float64[2,2]\nθu = Float64[4,4]\n\nbc = BoundConstraints(θl,θu)\n\n# solve it\nresult = solve(nls, θ_init, bc, conf)\n\n# use the returned result\n@assert converged(result)\nθ_solution = solution(result)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For furthers details see:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"solve(nls::AbstractNLS, θ_init::AbstractVector, conf::Abstract_Solver_Conf).\nsolve(nls::AbstractNLS, θ_init::AbstractVector,bc::BoundConstraints, conf::Abstract_Solver_Conf).","category":"page"}]
}