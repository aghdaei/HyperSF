include("Functions.jl")
include("../include/HyperLocal.jl")
include("../include/Helper_Functions.jl")
include("../include/maxflow.jl")


function HyperSF(Inp, L, R)

    ar = ReadInp(Inp)
    ar_org = copy(ar)
    mx = mxF(ar_org)

    ## Decomposition
    ar_new, idx_mat, SV = decomposition(ar, L)

    idx = idx_mat[end]
    ## Computing the effective resistance diameter of clusters and sort them
    dict = Dict{Any, Any}()
    @inbounds for jj =1:length(idx)

        vals = get!(Vector{Int}, dict, idx[jj])

        push!(vals, jj)

    end # for jj

    KS = collect(keys(dict))

    VL = collect(values(dict))

    score = zeros(eltype(SV), length(VL))
    @inbounds for i =1:length(VL)
        nodes = VL[i]
        for j in axes(SV, 2)
            mx, mn = -Inf, +Inf
            for node in nodes
                x = SV[node, j]
                mx = ifelse(x > mx, x, mx)
                mn = ifelse(x < mn, x, mn)
            end
            score[i] += (mx - mn)^2
        end
    end

    fdP = sortperm(score)
    RT = round(Int, R * length(fdP))
    fdC = fdP[1: RT]

    CL1 = KS[fdC]
    ## end of Computing the effective resistance diameter of clusters

    ## flow-based
    mx = mxF(ar_new)
    grow_lvl=3
    NH = HyperNodes(ar_new)
    H = INC(ar_new)
    flag2 = falses(mx)
    val = 1
    idx2 = zeros(Int, mx)
    grownum=3
    cond_thresh=1

    for ii = 1:length(fdC)

        seedN = CL1[ii]

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

                nnd = ar_new[HE[mm]]

                append!(new_nd, nnd)

            end #end of mm

            nd = sort(unique(new_nd))

        end #end of grow_lvl

        seedN2 = findall(x->in(x, seedN), nd)

        IM = H[HE, nd]

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
        Rs = findall(x->in(x,seedN2),R)  #Force

        S, lcond = HyperLocal(IM,IMt,order,d,R,epsilon,delta,Rs,true)

        volA = sum(d)

        S_org = nd[S]

        ## finding the not discovered nodes
        fgs = flag2[S_org]

        S_org2 = S_org[.!fgs]

        if length(S_org2) > 50

            S_org2 = S_org2[1:50]
        end


        if length(S_org2) > 0

            S3 = findall(x->in(x, S_org2), nd)

            condR,volR, cutR = tl_cond(IM,S3,d,delta,volA,order)

            if condR < cond_thresh

                idx2[S_org2] .= val

                flag2[S_org2] .= 1

                val+=1

            end

        end#end of if S_org2


    end #for ii


    ## indexing the isolated nodes
    fdz = findall(x-> x==0, idx2)

    fdnz = findall(x-> x!=0, idx2)

    V = vec(val:val+length(fdz)-1)

    idx2[fdz] = V

    push!(idx_mat, idx2)

    ## Mapper

    idx1 = 1:maximum(idx_mat[end])

    @inbounds for ii = length(idx_mat):-1:1

        idx1 = idx1[idx_mat[ii]]

    end # for ii


    ## global conductance

    ar = ar_org

    NH = HyperNodes(ar)

    H = INC(ar)

    order = vec(round.(Int, sum(H, dims = 2)))

    vol_V = vec(round.(Int, sum(H, dims = 1)))

    global Cvec = zeros(Float64, 0)


    dict = Dict{Any, Any}()

    count = 0
    @inbounds for jj =1:length(idx1)

        vals = get!(Vector{Int}, dict, idx1[jj])

        push!(vals, jj)

    end # for jj

    KS = collect(keys(dict))

    VL = collect(values(dict))

    global szCL_HE = zeros(Int, length(VL))

    ndV = collect(1:mxF(ar))

    @inbounds for ii = 1:length(KS)

        S = VL[ii]

        szCL_HE[ii] = length(S)

        ct = tl_cut(H, vec(S), 1.0, order)

        vol1 = sum(vol_V[S])

        S_hat = deleteat!(ndV, S)

        vol2 = sum(vol_V[S_hat])

        cnd = ct / vol1

        push!(Cvec, cnd)

        ndV = collect(1:mxF(ar))


    end # for ii


    NR = (mx_org - maximum(idx1)) / mx_org *100

    #ER = (LN_org - length(ar_new)) / LN_org *100

    println("NR = ", NR)

    #println("ER = ", ER)

    println("average conductance: ", round(mean(Cvec), digits=2))


    #Whgr("ibm05_t1.hgr", ar_new)

    return idx1


end # function
