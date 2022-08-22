using SparseArrays
using LinearAlgebra
using Clustering
using Distances
using Metis
using Laplacians
using StatsBase

include("HyperSF.jl")

filename = "ibm01.hgr"

## L controls the coarsening ratio by applying L-levels of k-mean clustering
L = 4

## R adjusts the ratio of selected clusters (low-quality clusters)
# for applying the flow-based technique
R = .1

## IDN is the index cluster that is assigned to each node
IDN = HyperSF(filename, L, R)

