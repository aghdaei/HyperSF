using SparseArrays
using LinearAlgebra
using Clustering
using Distances
using Metis
using Laplacians


include("HyperSF.jl")

Input = "../data/ibm01.hgr"

ar_coarse = HyperSF(Input)

# The incidence matrix corresponding to the coarsened hypergraph
H = INC3(ar_coarse)


## Generate the output hypergraph in hMetis format
# First line: #hyperedges, #nodes
Whgr("Output.hgr", ar_coarse)

