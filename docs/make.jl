using Documenter, JAMS

makedocs(sitename = "JAMS.jl", remotes = nothing,repo = Remotes.GitHub("sanderkammeraat", "JAMS.jl"),
format = Documenter.HTML(
        repolink = "https://github.com/sanderkammeraat/JAMS.jl",
        edit_link = "main",
        footer = "Part of this documentation is written with the help of Claude (Opus 5.5) and reviewed by the author. If you find something unclear or incorrect, please open a Github issue. Powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) and the [Julia Programming Language](https://julialang.org/).",
    ),
pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Overview of JAMS" => "overview.md",
        "Particles" => "particles.md",
        "Forces" => "forces.md",
        "DOF evolvers" => "dofevolvers.md",
        "Live plotting" => "live_plotting.md",
        "Save functions" => "save_functions.md",
        "Initial conditions" => "initial_conditions.md",
        "Extending JAMS" => "extending.md",
        "API reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/sanderkammeraat/JAMS.jl.git",devbranch = "main"
)