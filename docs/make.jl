using Documenter, Diversification

makedocs(
     sitename="Diversification",
     pages = [
          "Installation" => "install.md",
          "Rate analysis" => "analysis.md", 
          "index.md",
             ]
)

deploydocs(
    repo = "github.com/kopperud/Diversification.jl.git",
)


