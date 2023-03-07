using Documenter, Pesto

makedocs(
     sitename="Pesto.jl",
     pages = [
          "Installation" => "install.md",
          "Rate analysis" => "analysis.md", 
          "index.md",
             ]
)

deploydocs(
    repo = "github.com/kopperud/Pesto.jl.git",
)


