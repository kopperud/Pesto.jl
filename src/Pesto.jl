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
import ProgressMeter
import SparseArrays
import SciMLBase


## the rest
include("datatypes.jl")
include("treechecks.jl")
include("shiftbins.jl")
include("display.jl")

## pointer based tree
include("tree/types.jl")
include("tree/show.jl")
include("tree/construct.jl")
include("tree/print_indices.jl")
include("tree/height.jl")
include("tree/length.jl")
include("tree/utils.jl")


## utils
include("utils/utils.jl")
include("utils/logistic.jl")

## surface
include("loglsurface/surface.jl")

## constant birth-death process
include("bd_constant/constantBDP.jl")

## birth-death-shift process
include("bd_shift/postorder.jl")
include("bd_shift/postorder_sync.jl")
include("bd_shift/postorder_async.jl")
include("bd_shift/preorder.jl")
include("bd_shift/ODE.jl")
include("bd_shift/logLroot.jl")
include("bd_shift/tree_rates.jl")
include("bd_shift/backwards_forwards.jl")
include("bd_shift/birth_death_shift.jl")
include("bd_shift/extinction.jl")
include("bd_shift/nshifts.jl")
include("bd_shift/analysis.jl")
include("bd_shift/optimize_eta.jl")  ## optimize eta only
include("bd_shift/optimize_bd.jl") ## birth-death
include("bd_shift/optimize_bd_separate_shift.jl") ## birth-death with α and β, different shift rates
include("bd_shift/optimize_fbd.jl") ### fossilized birth-death
include("bd_shift/shift_probability.jl")
include("bd_shift/tip_rates.jl")
include("bd_shift/magnitude.jl")
include("bd_shift/shifts_along_branch.jl")

## input-output
include("io/writenewick.jl")
include("io/readnewick.jl")
include("io/readtree.jl")

## rcall
include("rcall/rconvert.jl")

## plot
include("plot/plot_tree.jl")

## custom exceptions
include("errors.jl")

## try some new interface
include("interface.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

## precompile
#PrecompileTools.@setup_workload begin
#    bears_tree = readtree(path("bears.tre"))
#
#end


end
