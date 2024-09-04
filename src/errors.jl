export PolytomyError, NegativeBranchError

#struct NonUltrametricError <: Exception end
struct PolytomyError <: Exception end

struct NegativeBranchError <: Exception end

Base.showerror(io::IO, e::PolytomyError) = print(io, "Your tree contains at least one hard polytomy. Pesto does not work with tress that have polytomies.")

Base.showerror(io::IO, e::NegativeBranchError) = print(io, "Your tree contains at least one branch with negative branch length, preventing Pesto from running.")
