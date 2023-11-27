module Pesto

import OrdinaryDiffEq
import RCall
import CSV
import Distributions
import Optim
import LinearAlgebra
import DataFrames
import FastGaussQuadrature
import LoopVectorization
import ForwardDiff
import RecipesBase
import ColorSchemes

## the rest
include("datatypes.jl")
include("utils.jl")
include("polytomy.jl")
include("shiftbins.jl")
include("display.jl")

## constant birth-death process
include("bd_constant/constantBDP.jl")

## birth-death-shift process
include("bd_shift/postorder.jl")
include("bd_shift/postorder_nosave.jl")
include("bd_shift/preorder.jl")
include("bd_shift/ODE.jl")
include("bd_shift/logLroot.jl")
include("bd_shift/tree_rates.jl")
include("bd_shift/backwards_forwards.jl")
include("bd_shift/birth_death_shift.jl")
include("bd_shift/extinction.jl")
include("bd_shift/nshifts.jl")
include("bd_shift/analysis.jl")
include("bd_shift/optimize.jl")
include("bd_shift/shift_probability.jl")
include("bd_shift/tip_rates.jl")

## input-output
include("io/writenewick.jl")
include("io/readnewick.jl")
include("io/readtree.jl")

## rcall
include("rcall/rconvert.jl")

## plot
include("plot/plot_tree.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

## precompile
#PrecompileTools.@setup_workload begin
#    bears_tree = readtree(path("bears.tre"))
#
#end


end
