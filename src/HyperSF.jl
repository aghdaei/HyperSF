include("Functions.jl")
include("../include/HyperLocal.jl")
include("../include/Helper_Functions.jl")
include("../include/maxflow.jl")

function HyperSF(Inp)

    ## Reading the input hypergraph
    ar = ReadInp(Inp)

    for in_loop = 1:1

        ## constructing the incidence matrix corresponding to the original hypergraph
        H = INC(ar)

        ## find the number of nodes in hypergraph
        mx = mxF(ar)

        ## geenrate a binary array (Flag) for all the nodes in hypergraph
        flag = falses(mx)

        ## converting hypergraph to simple graph using Star expansion
        MM = Star(ar)

        ## parameters 
        # changing these parameters allow you to generate different rresults
        szT = 10
        
        grow_lvl = 2

        grownum = 100
        ## end of parameters

        ## generating 10 random vectors
        RD = (rand(Float64, size(MM,1), 10) .- 0.5).*2

        ## smoothing steps
        k = 100

        ## Filtering random vectors by applying k smoothing steps
        SV = Filter(RD, k, MM, mx)

        ## Computing hyperedge scores
        Hscore = HSC(ar ,SV)

        ## clustering the nodes within each hyperedge applying k-mean
        idx_new = HEC(ar, SV)

        fdz = findall(x->x==0, idx_new)

        mx_idx = maximum(idx_new)

        vec1 = collect(mx_idx + 1 : mx_idx + length(fdz))

        idx_new[fdz] = vec1

        ###### Clusters dictionary #########
        dict = Dict{Any, Any}()

        count = 0
        for jj =1:length(idx_new)

            vals = get!(Vector{Int}, dict, idx_new[jj])

            push!(vals, jj)

        end # for jj

        KS = collect(keys(dict))

        VL = collect(values(dict))

        ##############  end of dictionary clusters ######################

        ##### Sort the clusters according to their sizes

        szCL = zeros(Int, length(VL))

        for ii = 1:length(VL)
            szCL[ii] = length(VL[ii])
        end

        spC1 = findall(x->x< szT, szCL)
        #####

        flag2 = falses(length(idx_new))

        idx2 = zeros(Int, length(idx_new))

        val = 1

        ## Flow-baed method
        for ii = 1:length(spC1)

            seedN = VL[spC1[ii]]

            ## expanding the network around the seed nodes
            nd = copy(seedN)

            HE = Any[]

            for ee = 1:grow_lvl

                HE = Any[]

                for dd = 1:length(nd)

                    append!(HE, NH[nd[dd]])

                end

                HE = sort(unique(HE))

                new_nd = Any[]

                for mm=1:length(HE)

                    nnd = ar[HE[mm]]

                    append!(new_nd, nnd)

                end #end of mm

                nd = sort(unique(new_nd))

            end #end of grow_lvl

            seedN2 = findall(x->in(x, seedN), nd)

            IM = H[HE, nd]
            #println("size(IM): ", size(IM))

            IMt = sparse(IM')

            ## d is node degree
            d = vec(sum(IM,dims=1))

            epsilon = 1.0

            delta = 1.0

            order = round.(Int, sum(IM, dims = 2))
            order = order[:,1]


            OneHop = get_immediate_neighbors(IM,IMt,seedN2)
            Rmore = BestNeighbors(IM,d,seedN2,OneHop,grownum)
            R = union(Rmore,seedN2)

            Rs = FD3(x->in(x,seedN2),R)  #Force

            S, lcond = HyperLocal(IM,IMt,order,d,R,epsilon,delta,Rs,true)

            volA = sum(d)

            S_org = nd[S]

            ## finding the not discovered nodes
            fgs = flag2[S_org]

            S_org2 = S_org[.!fgs]

            if length(S_org2) > 10

                S_org2 = S_org2[1:10]
            end


            if length(S_org2) > 0

                S3 = findall(x->in(x, S_org2), nd)

                condR,volR, cutR = tl_cond(IM,S3,d,delta,volA,order)

                idx2[S_org2] .= val

                flag2[S_org2] .= 1

                val+=1

            end#end of if S_org2


        end #for ii


        ## mappers
        fdz2 = findall(x->x==0, idx2)
        fdnz2 = findall(x->x!=0, idx2)

        dict2 = Dict{Any, Any}()

        count = 0
        for jj =1:length(fdz2)

            vals = get!(Vector{Int}, dict2, idx_new[fdz2[jj]])

            push!(vals, fdz2[jj])

        end # for jj

        KS2 = collect(keys(dict2))

        VL2 = collect(values(dict2))

        mx_idx2 = maximum(idx2[fdnz2])

        ###### Clusters dictionary - final clusters #########
        dict3 = Dict{Any, Any}()

        count = 0
        for jj =1:length(idx2)

            vals = get!(Vector{Int}, dict3, idx2[jj])

            push!(vals, jj)

        end # for jj

        KS3 = collect(keys(dict3))

        VL3 = collect(values(dict3))

        global szCL3 = zeros(Int, length(VL3))

        for ii = 1:length(VL3)
            szCL3[ii] = length(VL3[ii])
        end

        ##############  end of dictionary clusters ######################
        val2 = 1

        for ii = 1:length(VL2)

            ps1 = VL2[ii]

            idx2[ps1] .= val2 + mx_idx2

            val2 +=1

        end # for ii

        ar_new = Any[]

        for i =1:length(ar)

            nd = ar[i]

            new_nd = idx2[nd]

            push!(ar_new, sort(unique(new_nd)))


        end #end of for i

        ar_new = unique(ar_new)

        ### removing hyperedges with cardinality of 1
        HH = INC(ar_new)
        ss = sum(HH, dims=2)
        fd1 = findall(x->x==1, ss[:,1])
        deleteat!(ar_new, fd1)

        ar = ar_new


    end 

    return ar
end 
