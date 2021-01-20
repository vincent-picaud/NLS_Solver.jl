using Kunisch
using Documenter

makedocs(;
    modules=[Kunisch],
    authors="Picaud Vincent <picaud.vincent@gmail.com> and contributors",
    repo="https://github.com/vincent-picaud/Kunisch.jl/blob/{commit}{path}#L{line}",
    sitename="Kunisch.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vincent-picaud.github.io/Kunisch.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vincent-picaud/Kunisch.jl",
)
