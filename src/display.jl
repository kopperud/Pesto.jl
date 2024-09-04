#import Base.Multimedia.display


function Base.Multimedia.display(phy::phylo)
    n = phy.Nnode
    ntip = length(phy.tip_label)
    println("Rooted phylogenetic tree with $ntip tips and $n internal nodes.\n")
    println("Tip labels:")
    print("\t")
    for i in 1:minimum((ntip, 6))
        print(phy.tip_label[i])
        print(", ")
    end
    print("...")
end

function Base.Multimedia.display(data::SSEdata)
    n = data.Nnode
    ntip = length(data.tiplab)
    sampling_probability = data.sampling_probability

    println("Rooted phylogenetic tree with $ntip tips and $n internal nodes, and sampling fraction sampling_probability = $sampling_probability.\n")
    println("Tip labels:")
    print("\t")
    for i in 1:minimum((ntip, 6))
        print(data.tiplab[i])
        print(", ")
    end
    print("...")
end
