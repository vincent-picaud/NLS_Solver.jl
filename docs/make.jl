using NLS_Solver
using Documenter

makedocs(;
    modules=[NLS_Solver],
    authors="Picaud Vincent <picaud.vincent@gmail.com> and contributors",
    repo="https://github.com/vincent-picaud/NLS_Solver.jl/blob/{commit}{path}#L{line}",
    sitename="NLS_Solver.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vincent-picaud.github.io/NLS_Solver.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vincent-picaud/NLS_Solver.jl",
)
