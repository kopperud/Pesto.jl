module Diversification

import OrdinaryDiffEq
import RCall
import CSV
import ProgressMeter
import Statistics
import DataStructures
import Distributions
import Optim
import ForwardDiff
import LinearAlgebra

include("datatypes.jl")
include("optimize.jl")
include("utils.jl")
include("diversitree_bisse.jl")
include("postorder.jl")
include("postorder_nosave.jl")
include("postorder_chunk.jl")
include("preorder.jl")
include("ODE.jl")
include("logLroot.jl")
include("constantBDP.jl")
include("tree_rates.jl")
include("backwards_forwards.jl")
include("birth_death_shift.jl")
include("extinction.jl")
include("writenewick.jl")
include("nshifts.jl")



# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)
#path(x::String) = joinpath(@__DIR__, x)

# Write your package code here.
end
