using Documenter, Pesto, Makie, CairoMakie

CairoMakie.activate!(type="svg")

makedocs(
     sitename="Pesto.jl",
     authors = "Bjørn Tore Kopperud and Sebastian Höhna",
     modules = [
      Pesto,
      isdefined(Base, :get_extension) ? Base.get_extension(Pesto, :PestoMakieExt) :
      Pesto.PestoMakieExt
      ],
     pages = [
          "Home" => "index.md",
          "Installation" => "install.md",
          "Analyses" => [
            "Simple analysis" => "analysis/simple.md", 
            "Extended analysis" => "analysis/extended.md",
            "Number of shifts" => "analysis/shifts.md",
            "Significant shifts" => "analysis/bayes_factors.md",
            "Tip rates" => "analysis/tiprates.md",
            "Plot with ggtree" => "plotting/ggtree.md"
          ],
          "Functions" => "functions.md",
             ]
)

deploydocs(
    repo = "github.com/kopperud/Pesto.jl.git",
)


