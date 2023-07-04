using Documenter, Pesto

makedocs(
     sitename="Pesto.jl",
     authors = "Bjørn Tore Kopperud and Sebastian Höhna",
     modules = [Pesto],
     pages = [
          "Home" => "home.md",
          "Installation" => "install.md",
          "Analyses" => [
            "Simple analysis" => "analysis/simple.md", 
            "Extended analysis" => "analysis/extended.md",
            "Number of shifts" => "analysis/shifts.md",
            "Tip rates" => "analysis/tiprates.md",
          ],
          "Functions" => "index.md",
             ]
)

deploydocs(
    repo = "github.com/kopperud/Pesto.jl.git",
)


