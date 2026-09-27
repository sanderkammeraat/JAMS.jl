using Documenter, JAMS

makedocs(sitename = "JAMS.jl", remotes = nothing,repo = Remotes.GitHub("sanderkammeraat", "JAMS.jl"),
format = Documenter.HTML(
        repolink = "https://github.com/sanderkammeraat/JAMS.jl",
        edit_link = "main",
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