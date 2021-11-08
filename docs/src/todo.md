# TODO List

## DONE deploy key

Added (laptop: id_solver_jl_doc_deploykey.pub)

ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIGgAnccrpKqTRP4ykDQfkW90wB1yKe+A2Fr04KzCgmnT picaud.vincent@gmail.com

## Rename 
- eval_nls_folk
- eval_nls_∇fobj!
- eval_nls_∇fobj!

eval -> compute

## Factorize ResultTypes

Our AbstractNLSResult, AbstractNLSBCResult, AbstractQuadResult result
structure all have soem common methods like `converged()`,
`iteration_count()` or `solution()`

-> one must define and use a **Base** type with these **common**
functions.

### Generalize to some other types?

-> have a look a configuration type if a possible factorization is
possbile.

