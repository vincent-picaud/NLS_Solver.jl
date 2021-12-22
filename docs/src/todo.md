# TODO List

## Factorize ResultTypes

Our `Abstract_Solver_Result`, `Abstract_BC_Solver_Result`, `AbstractQuadResult` result
structures all have some common methods like `converged()`,
`iteration_count()` or `solution()`

-> one must define and use a **Base** type with these **common**
functions.


