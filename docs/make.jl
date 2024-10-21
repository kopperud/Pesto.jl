using Documenter, Pesto, Makie, CairoMakie, ColorSchemes, DataFrames 

CairoMakie.activate!(type="svg")

modules = Module[Pesto]
y = Base.get_extension(Pesto, :PestoMakieExt)
if !isnothing(y)
    push!(modules, y)
end

makedocs(
     sitename="Pesto.jl",
     authors = "Bjørn Tore Kopperud and Sebastian Höhna",
     modules = modules,
     pages = [
          "Home" => "index.md",
          "Installation" => "install.md",
          "Analyses" => [
            "Simple analysis" => "analysis/simple.md", 
            "Two-step analysis" => "analysis/twostep.md",
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


