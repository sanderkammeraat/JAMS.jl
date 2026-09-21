using Documenter, JAMS

makedocs(sitename = "JAMS.jl", remotes = nothing,repo = Remotes.GitHub("sanderkammeraat", "JAMS.jl"),
format = Documenter.HTML(
        repolink = "https://github.com/sanderkammeraat/JAMS.jl",
        edit_link = "main", 
    )
)

deploydocs(
    repo = "github.com/sanderkammeraat/JAMS.jl.git",devbranch = "main"
)