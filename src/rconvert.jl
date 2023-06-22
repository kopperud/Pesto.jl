import RCall: sexp, protect, unprotect, setclass!, RClass, sexpclass

RCall.sexpclass(::SSEdata) = RClass{:phylo}

function sexp(::Type{RClass{:phylo}}, f::SSEdata)
    phy = protect(sexp(Dict(
        "edge" => f.edges, 
        "Nnode" => f.Nnode,
        "tip.label" => f.tiplab,
        "edge.length" => f.branch_lengths
        )))
    setclass!(phy, sexp("phylo"))
    unprotect(1)
    phy
end

# function sexp(::Type{RClass{:treedata}}, f::SSEresult)
#     td = protect(sexp(
#         Dict(
#             "phy" => Dict(
#                 "edge" => f.phy["edge"], 
#                 "Nnode" => f.phy["Nnode"],
#                 "tip.label" => f.phy["tip.label"],
#                 "edge.length" => f.phy["edge.length"]
#                 ),
#             "data" => f.rates
#         )))
#     #setclass!(td, sexp("treedata"))
#     RCall.R"""
#         td <- new("treedata", data = )
#     """
#     unprotect(1)

#   #  x = rcopy(reval("""
#         #setClass("Foo", representation(x = "numeric"))
#    #     foo <- new("Foo", x = 20)
#    # """))
#     td
# end
