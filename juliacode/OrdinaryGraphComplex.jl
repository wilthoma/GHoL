
# python compatibility
len(lst) = length(lst)

ordinaryDataDirOdd = DATA_DIR + "/ordinarydata/oddedge/"
ordinaryDataDirEven = DATA_DIR + "/ordinarydata/evenedge/"
imgBaseDir = "img/"

type OrdinaryGraphVectorSpace <: GraphVectorSpace{SmallGraph}
  nVertices :: Int
  nLoops :: Int
  evenEdges :: Bool
end

function get_file_name(self::OrdinaryGraphVectorSpace)
  dataDir = self.evenEdges ? ordinaryDataDirEven : ordinaryDataDirOdd
  s = @sprintf "gra%d_%d.g6" self.nVertices self.nLoops
  return string(dataDir, s)
end

function get_svg_dir(self::OrdinaryGraphVectorSpace)
  dataDir = self.evenEdges ? ordinaryDataDirEven : ordinaryDataDirOdd
  s = @sprintf "imgs%d_%d/" self.nVertices self.nLoops
  return string(dataDir, imgBaseDir, s)
end


function get_color_counts(self::OrdinaryGraphVectorSpace)
        return Void # no coloring
end

function get_fstring(self::OrdinaryGraphVectorSpace)
        return "z"
end

"""Check whether the T graded subspace can in principle be non-empty.
   For each such T, a file is expected. Otherwise the corresponding
   graded component is considered not computed yet."""
function is_valid(self::OrdinaryGraphVectorSpace)
        nEdges = self.nLoops + self.nVertices -1
        # at least trivalent, and simple
        return  (3*self.nVertices <= 2*nEdges) && self.nVertices > 0 && self.nLoops >= 0 && nEdges <= self.nVertices*(self.nVertices-1)/2
end

"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs(self::OrdinaryGraphVectorSpace)
        # generate List of unmarked graphs
        println(self.nVertices)
    LL = listG(self.nVertices, self.nLoops)
    return LL
end

"""For G a graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
   Here vertex j becomes vertex p[j] in the new graph."""
function get_perm_sign(self::OrdinaryGraphVectorSpace, G::SmallGraph, p)
    nVert, nLoops, evenEdges = (self.nVertices, self.nLoops, self.evenEdges)
    nEdges = nLoops + nVert - 1
    if evenEdges
        # The sign is (induced sign on vertices) * (induced sign edge orientations)
        sgn = permSign(p)
        for (u,v) in edges(G)
            # we assume the edge is always directed from the larger to smaller index
            if (u < v && p[u] > p[v]) || (u > v && p[u] < p[v])
                sgn *= -1
            end
        end
        return sgn
    else
        # The sign is (induced sign of the edge permutation)
        # we assume the edges are always lex ordered
        # for the computation we use that edges(G) returns the edges in lex ordering
        pp = collect(1:nEdges)
        G2 = permuteGraph(G,p)
        edge_idx2 = -ones(Int64, nVert, nVert)
        for (j, e) = enumerate(edges(G2))
            u,v = e
            edge_idx2[u,v] = j
            edge_idx2[v,u] = j
        end
        for (j, e) = enumerate(edges(G))
            u,v = e
            pp[j] = edge_idx2[p[u],p[v]]
        end

        #println(pp)
        return permSign(pp)
    end
end

"""Converts the graph to a graphviz dot format string.
   This method is used only for visualization, not for computation."""
function get_dot(self::OrdinaryGraphVectorSpace, G)
    return render_to_dot(G)
end


# -----  Contraction operator --------

type ContractDOrdinary <: GraphOperator{SmallGraph,SmallGraph}
  # source
  nVertices :: Int
  nLoops :: Int
  evenEdges :: Bool
end

"""Retrieve the file name of the file storing the matrix list for the operator."""
function get_file_name(self::ContractDOrdinary)
  dataDir = self.evenEdges ? ordinaryDataDirEven : ordinaryDataDirOdd
  s = @sprintf "contractD%d_%d.txt" self.nVertices self.nLoops
  return string(dataDir, s)
end

function get_unique_file_name(self::ContractDOrdinary)
  prefix = self.evenEdges ? "ordinary_even_" : "ordinary_odd_"
  s = @sprintf "contractD%d_%d.sms" self.nVertices self.nLoops
  return string(prefix, s)
end

"""Returns target graph vector space."""
function get_target(self::ContractDOrdinary)
  return OrdinaryGraphVectorSpace(self.nVertices-1, self.nLoops, self.evenEdges)
end

"""
Returns the GraphVectorSpace on which this operator acts
"""
function get_source(self::ContractDOrdinary)
  return OrdinaryGraphVectorSpace(self.nVertices, self.nLoops, self.evenEdges)
end

"""For G a graph returns a list of pairs (GG, x),
   such that (operator)(G) = sum x GG.
"""
function operate_on(self::ContractDOrdinary, G::SmallGraph)
    vs = get_source(self)

    ret=[]
    for (u,v) = edges(G)
        # permute u,v to position 1 and 2
        p = collect(1:self.nVertices)
        p[1] = u
        p[2] = v
        idx = 3
        for j = 1:self.nVertices
            if j == u || j== v
                continue
            else
                p[idx] = j
                idx +=1
            end
        end
        pp = invPermutation(p)
        #println(pp)
        sgn = get_perm_sign(vs,G, pp)
        GG = permuteGraph(G,pp)

        # now delete the first edge
        remove_edge!(GG,1,2) #.... done later
        # ... and join 0 and 1 and fix labelings
        p = vcat([1],collect(1:self.nVertices-1))
        if !self.evenEdges
          sgn *= get_perm_sign(vs,GG, p)  # TODO: no good to call get_perm_sign with non-permutation if evenEdges
        end
        GG = permuteGraph(GG,p)

        # finally delete the last vertex
        # remove_vertex(GG, nVert) ... by creating a new graph
        GGG = small_graph(self.nVertices-1)
        for (uu,vv) = edges(GG)
          if uu != vv
            add_edge!(GGG, uu, vv)
          end
        end

        push!(ret, (GGG, sgn))
    end
    return ret
end

function get_work_estimate(self::OrdinaryGraphVectorSpace)
  # give estimate of number of graphs
  nEdges = self.nLoops + self.nVertices -1
  n = self.nVertices
  return binomial(div(n*(n-1),2), nEdges) / factorial(n)
end

function get_work_estimate(self::ContractDOrdinary)
  # give estimate of number of graphs
  vs = get_source(self)
  nEdges = vs.nLoops + vs.nVertices -1
  return get_work_estimate(vs) * nEdges
end

function dispListCoverageOrdinary(nDisplay=0)
  data = []
  nVR = collect(3:12)
  nLR=collect(2:15)

  for evenEdges in [true, false]
    curdata = Array{Any}(length(nVR), length(nLR))
    for (i,v) in enumerate(nVR)
      for (j,l) in enumerate(nLR)
        vs = OrdinaryGraphVectorSpace(v,l,evenEdges)
        dim = getDimension(vs)
        deco = is_valid(vs) ? "" : "class=redcell"
        curdata[i,j] = Dict("data"=>dim, "style"=>deco)
      end
    end
    push!(data, curdata)
  end

  dispTables(["Even edges (vertices\\loops)", "Odd edges (vertices\\loops)"], Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)

end

function dispCohomologyOrdinary(nDisplay=0)
  data = []
  nVR = collect(3:8)
  nLR=collect(2:8)

  for evenEdges in [true, false]
    curdata = Array{Any}(length(nVR), length(nLR))
    for (i,v) in enumerate(nVR)
      for (j,l) in enumerate(nLR)
        vs = OrdinaryGraphVectorSpace(v,l,evenEdges)
        op = ContractDOrdinary(v,l,evenEdges)
        op2 = ContractDOrdinary(v+1,l,evenEdges)

        dim = get_cohomology(op, op2)
        deco = is_valid(vs) ? "" : "class=redcell"
        curdata[i,j] = Dict("data"=>dim, "style"=>deco)
      end
    end
    push!(data, curdata)
  end

  dispTables(["Cohomology Even edges (vertices\\loops)", "Cohomology Odd edges (vertices\\loops)"], Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)

end

"""
  Shows a table of the number of nonzero entries in operator.
  Or -1 if operator could not be loaded.
"""
function dispOperatorCoverageOrdinary(nDisplay=0)
  data = []
  nVR = collect(3:12)
  nLR=collect(2:15)

  for evenEdges in [true, false]
    curdata = Array{Any}(length(nVR), length(nLR))
    for (i,v) in enumerate(nVR)
      for (j,l) in enumerate(nLR)
        vs = OrdinaryGraphVectorSpace(v,l,evenEdges)
        theop=ContractDOrdinary(v,l,evenEdges)
        nrEntries = -1
        try
          D = load_matrix(theop)
          if D==[]
            nrEntries = 0
          else
            nrEntries = nnz(D)
          end
        catch
        end
        deco = is_valid_op(theop) ? "" : "class=redcell"

        curdata[i,j] = Dict("data"=>nrEntries, "style"=>deco)
      end
    end
    push!(data, curdata)
  end

  dispTables(["Even edges (vertices\\loops)", "Odd edges (vertices\\loops)"], Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)

end


#@implements OrdinaryGraphVectorSpace <: GraphComplex
