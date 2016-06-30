
"""
    This file contains the implementation of the GraphVectorSpace interface for connected, >= trivalent graphs
    with marked edges, such that removing the edges makes the graph int a forest.

    There is a four-grading:
       (nVert, nLoops, nMarks, evenEdges)
    with nVert, nLoops, nMarks integers and evenEdges boolean
""";

using Memoize

markedDataDirOdd = "markeddata/oddedge/"
markedDataDirEven = "markeddata/evenedge/"
imgBaseDir = "img/"

# ----- helper methods
"""
Converts a marked multigraph to HackG format, i.e., adds one two-valent
vertex in the middle of each edge, with the marked edges standing at the end.
Tadpoles (if any) are encoded by a (necessarily marked) univalent vertex,
attached to the vertex at which the tadpole is.
"""
function multiGtoHackG_marked(G::MarkedSmallMultiGraph)
    # same as before, but the extra vertices corresponding to the marked edges are now at the end
    nEdges = num_edges(G) # only unmarked
    nMarks = num_marked_edges(G)
    nVerts = num_vertices(G)
    retG = small_graph(nVerts + nEdges + nMarks )
    for (i, (u,v)) in enumerate(edges(G))
        add_edge!(retG,u,nVerts+i)
        add_edge!(retG,v,nVerts+i)
    end
    for (i, (u,v)) in enumerate(marked_edges(G))
        add_edge!(retG,u,nVerts+nEdges+i)
        add_edge!(retG,v,nVerts+nEdges+i)
    end
    retG
end

"""
The reverse conversion from HackG format to marked multigraph.
Note that the data of the number of marks is not part of the HackG graph and
has to be provided.
"""
function hackGtoMultiG_marked(G::SmallGraph, nMarks)
    n = num_vertices(G)
    degs = [out_degree(v,G) for v in 1:n]
    trivalentVertices = find(d->d>=3, degs)
    invTrivalentVertices=zeros(Int, n)
    for (i,v) in enumerate(trivalentVertices)
        invTrivalentVertices[v] = i
    end

    nTrueVertices = length(trivalentVertices)
    nTotalEdges = n-nTrueVertices
    retG = marked_small_multi_graph(nTrueVertices)
    for u=1:n-nMarks
        if degs[u] ==2
            a,b = out_neighbors(u, G)
            println( "   $a $b")
            add_edge!(retG,invTrivalentVertices[a],invTrivalentVertices[b])
        end # note: there are no non-marked tadpoles, so no need to add them here
    end
    for u=n-nMarks+1:n
      if degs[u] == 2
        a,b = out_neighbors(u, G)
        add_marked_edge!(retG, invTrivalentVertices[a],invTrivalentVertices[b])
      elseif degs[u]==1 # tadpole
        a = out_neighbors(u, G)
        add_marked_edge!(retG,invTrivalentVertices[a], invTrivalentVertices[a])
      else
        error("marked vertices in HackG must have valence 1 or 2")
      end
    end
    return retG
end

#@memoize
nChooseMSetsMemoizer = Dict{Tuple{Int,Int}, Vector{Any}}()
function nChooseMSets(n::Int,m::Int)
    # lists all n-lists of bools, of which exactly m are True
    k = (n,m)
    global nChooseMSetsMemoizer
    if haskey(nChooseMSetsMemoizer, k)
      return copy(nChooseMSetsMemoizer[k])
    else
      ret = Any[]
      if m==0
          ret = Any[ [false for i in 1:n] ]
      elseif m==n
          ret = Any[ [true for i in 1:n] ]
      elseif m>n || m<0
          ret = Any[]
      else
          ret = vcat( [vcat(S,[true]) for S in nChooseMSets(n-1,m-1)] ,  [vcat(S,[false]) for S in nChooseMSets(n-1,m)]   )
      end
      nChooseMSetsMemoizer[k] = copy(ret)
      return ret
    end
end

#@memoize
function nChooseMSets2(n::Int,m::Int)
    # lists all n-lists of bools, of which exactly m are True
    if m==0
        return Any[ [false for i in 1:n] ]
    elseif m==n
        return Any[ [true for i in 1:n] ]
    elseif m>n || m<0
        return Any[]
    else
        return vcat( [vcat(S,[true]) for S in nChooseMSets(n-1,m-1)] ,  [vcat(S,[false]) for S in nChooseMSets(n-1,m)]   )
    end
end

#@profile
"""
Produces a marked multigraph by marking the edges of G according to l.
"""
function markEdges(G::SmallMultiGraph, l::Vector{Bool})
    GG = marked_small_multi_graph(num_vertices(G))
    for (i , (u,v)) in enumerate(edges(G))
       if l[i]
          add_marked_edge!(GG,u,v)
       else
          add_edge!(GG,u,v)
       end
    end
    return GG
end

# uses suboptimal O(n^2) algorithm, however without recursion
function myIsForest(G, excludedEdges)
    n=num_vertices(G)
    color_tbl = [ u => i for (i,u) in enumerate(vertices(G))]
    for (k,EE) in enumerate(edges(G))
        if !excludedEdges[k]
            u,v = EE
            c1 = color_tbl[u]
            c2 = color_tbl[v]
            if c1 == c2
                return false
            else
                new_color = min(c1, c2)
                for j = 1:n
                    if (color_tbl[j]==c1 || color_tbl[j] == c2)
                        color_tbl[j] = new_color
                    end
                end
            end
        end
    end
    return true
end

function isCutForest(G,l)
    # checks whether removing edges j s.t. l[j]=true makes G into a forest
    #GG = myCopyMultiGraph(G)
    #GG.remove_edges_from([e for i,e in enumerate(GG.edges(keys=True)) if l[i]])
    #aa = nt.is_forest(GG)
    if length(l) != num_edges(G)
      error("isCutForest called with incorrect parameters")
    end
    bb = myIsForest(G, l)
    #if aa != bb:
    #    print "Error!!!!!"
    return bb
end

function listMarkingsOfG(G::SmallMultiGraph, nrMarkings::Int)
    if nrMarkings < 0
        return []
    else
        S = nChooseMSets(num_edges(G),nrMarkings)
        filter!(l -> isCutForest(G,l), S)
        return [markEdges(G, l) for l in S]
    end
end


# ----- interface implementations

type MarkedEdgeGraphVectorSpace <: GraphVectorSpace{SmallGraph}
  nVertices :: Int
  nLoops :: Int
  nMarks :: Int
  evenEdges :: Bool
end

"""Retrieve the file name of the file storing the matrix list for the T graded component."""
function get_file_name(self::MarkedEdgeGraphVectorSpace)
        dataDir = self.evenEdges ? markedDataDirEven : markedDataDirOdd
        s = @sprintf "gra%d_%d_%d.g6" self.nVertices self.nLoops self.nMarks
        return string(dataDir, s)
end

function get_svg_dir(self::MarkedEdgeGraphVectorSpace)
    dataDir = self.evenEdges ? markedDataDirEven : markedDataDirOdd
    s = @sprintf "imgs%d_%d_%d/" self.nVertices self.nLoops self.nMarks
    return string(dataDir, imgBaseDir, s)
end

"""Internally, each graph is considered colored. This method returns a
   vector [a,b,c,...] such that the first a vertices are in the first color,
   the second b in the second etc."""
function get_color_counts(self::MarkedEdgeGraphVectorSpace)
        nEdges = self.nLoops + self.nVertices -1
        return color_counts_to_nauty_fmt( [self.nVertices+nEdges-self.nMarks, self.nMarks] )
end

function get_fstring(self::MarkedEdgeGraphVectorSpace)
  nEdges = self.nLoops + self.nVertices -1
        return "a"^(self.nVertices+nEdges-self.nMarks)
end

"""Check whether the T graded subspace can in principle be non-empty.
   For each such T, a file is expected. Otherwise the corresponding
   graded component is considered not computed yet."""
function is_valid(self::MarkedEdgeGraphVectorSpace)
        nEdges = self.nLoops + self.nVertices -1
        # at least trivalent, marks break graph into forest, and each edge can have at most one mark
        return  (3*self.nVertices <= 2*nEdges) &&  (self.nMarks >= self.nLoops) && (self.nMarks <=nEdges)
end

"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs_nomulti(self::MarkedEdgeGraphVectorSpace)
        # generate List of unmarked graphs
        LL = listG(self.nVertices, self.nLoops)
        LL = [small_graph_to_multigraph(G) for G in LL]
        LL = [listMarkingsOfG(G, self.nMarks) for G in LL]
        LL = vcat(LL...)
        return [multiGtoHackG_marked(G) for G in LL]
end

function tadpolify(G, nTadpoles)
  ret = tadpolify_int(G, nTadpoles,1)
  #println("Tadpoling $nTadpoles: $(length(ret)) created")
  return ret
end
function tadpolify_int(G, nTadpoles, nStartVertex)
  if nTadpoles==0
    return Any[G]
  end
  n = num_vertices(G)
  if nStartVertex == n
    GG = copy(G)
    GG[n,n] = nTadpoles
    return Any[GG]
  end
  ret = Any[]
  for j=0:nTadpoles
    pre = tadpolify_int(G, nTadpoles-j, nStartVertex+1)
    for GG in pre
      GGG = copy(GG)
      GGG[nStartVertex,nStartVertex] = j
      push!(ret, GGG)
    end
  end
  return ret
end


"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs(self::MarkedEdgeGraphVectorSpace)
        # generate List of unmarked graphs
        LL=Any[]
        nEdges = self.nVertices + self.nLoops - 1
        if self.evenEdges # no tadpoles
            LL = listMultiG(self.nVertices, self.nLoops, false)
        else
          for nTadpoles = 0:(nEdges-self.nVertices+1)
            LLL = listMultiGex(self.nVertices, self.nLoops-nTadpoles)
            LLL2 = [tadpolify(G, nTadpoles) for G in LLL]
            LLL3 = vcat(LLL2...)
            filter!(g-> (minimum_out_degree(g) >=3), LLL3)
            append!(LL, LLL3)
            #println("   __ $nTadpoles $(length(LLL)) $(length(LLL2)) $(length(LLL3))" )
          end
          for g in LL
            println("   Min outdegree $(minimum_out_degree(g))")
          end
        end
        #dispGraphList(LL)
        #error("stop here")

        LL = [listMarkingsOfG(G, self.nMarks) for G in LL]
        LL = vcat(LL...)
        #dispGraphList(LL, nDisplay = 2)
        return [multiGtoHackG_marked(G) for G in LL]
end

"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs_notp(self::MarkedEdgeGraphVectorSpace)
        # generate List of unmarked graphs
        LL = listMultiG(self.nVertices, self.nLoops, false)
        LL = [listMarkingsOfG(G, self.nMarks) for G in LL]
        LL = vcat(LL...)
        return [multiGtoHackG_marked(G) for G in LL]
end

"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs_orig(self::MarkedEdgeGraphVectorSpace)
        # generate List of unmarked graphs
        LL = listMultiG(self.nVertices, self.nLoops)
        LL = [listMarkingsOfG(G, self.nMarks) for G in LL]
        LL = vcat(LL...)
        return [multiGtoHackG_marked(G) for G in LL]
end

"""For G a (HackG-)graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
   Here vertex j in G becomes vertex p[j] in the new (relabeled) graph."""
function get_perm_sign(self::MarkedEdgeGraphVectorSpace, G::SmallGraph, p)
        nVert, nLoops, nMarks, evenEdges = (self.nVertices, self.nLoops, self.nMarks, self.evenEdges)
        nEdges = div(num_edges(G), 2) # two edges in HackG encode one edge
        nHackVerts = num_vertices(G)
        nVerts = nHackVerts - nEdges

        if evenEdges
            # The sign is (induced sign on vertices) * (induced sign on marked edges) * (induced sign edge orientations)
            markedEdges = collect(nHackVerts-nMarks+1 : nHackVerts)
            sgnMarkedEdges = permSign(inducedPerm(p, markedEdges))
            trueVertices = find(v-> out_degree(v,G) >=3, vertices(G))
            sgnVertices = permSign(inducedPerm(p, trueVertices))
            allEdges = find(v-> out_degree(v,G) == 2, vertices(G))
            sgnDirections=1
            for v in allEdges
                nb = out_neighbors(v,G)
                # we assume the edge is always directed from the larger to smaller index
                if (nb[1] < nb[2] && p[nb[1]] > p[nb[2]]) || (nb[1] > nb[2] && p[nb[1]] < p[nb[2]])
                    sgnDirections *= -1
                end
            end
            return sgnMarkedEdges * sgnVertices * sgnDirections
        else
            # The sign is (induced sign nonmarked edge permutations)
            nonMarkedEdges = find(v-> out_degree(v,G) == 2 && v<=nHackVerts-nMarks, vertices(G)) # mind that edges are encoded by 2-valent verts in HackG
            return permSign(inducedPerm(p, nonMarkedEdges))
        end
end

"""Converts the graph to a graphviz dot format string.
   This method is used only for visualization, not for computation."""
function get_dot(self::MarkedEdgeGraphVectorSpace, G)
        GG = hackGtoMultiG_marked(G, self.nMarks)
        ret = "graph {\n"
        for (u,v) in edges(GG)
            ret = string(ret, "$u -- $v;\n")
        end
        for (u,v) in marked_edges(GG)
            ret = string(ret, "$u -- $v[color=\"black:invis:black\"];\n")
        end
        ret = string(ret, "}")
        return ret
end

function get_work_estimate(self::MarkedEdgeGraphVectorSpace)
  # give estimate of number of graphs
  nEdges = self.nLoops + self.nVertices -1
  n = self.nVertices
  return binomial(div(n*(n-1),2), nEdges) * binomial(nEdges,self.nMarks) / factorial(n)
end

type MarkD <: GraphOperator{SmallGraph,SmallGraph}
  nVertices :: Int
  nLoops :: Int
  nMarks :: Int
  evenEdges :: Bool
end

function get_work_estimate(self::MarkD)
  # give estimate of number of graphs
  vs = get_source(self)
  nEdges = vs.nLoops + vs.nVertices -1
  return get_work_estimate(vs) * nEdges
end

"""Retrieve the file name of the file storing the matrix list for the T graded component.
   Here and below T is the grading of the source."""
function get_file_name(self::MarkD)
  dataDir = self.evenEdges ? markedDataDirEven : markedDataDirOdd
  s = @sprintf "markD%d_%d_%d.txt" self.nVertices self.nLoops self.nMarks
  return string(dataDir, s)
end

"""The operator maps the T graded piece to the T' graded piece.
   This functions returns T' given T."""
function get_target(self::MarkD)
   return MarkedEdgeGraphVectorSpace(self.nVertices, self.nLoops, self.nMarks+1, self.evenEdges)
end

"""
Returns the GraphVectorSpace on which this operator acts
"""
function get_source(self::MarkD)
  return MarkedEdgeGraphVectorSpace(self.nVertices, self.nLoops, self.nMarks, self.evenEdges)
end

"""For G a (HackG-format) graph, returns a list of pairs (GG, x),
   such that (operator)(G) = sum x GG."""
function operate_on(self::MarkD, G)
    nVert,nLoops,nMarks,evenEdges = (self.nVertices,self.nLoops, self.nMarks, self.evenEdges)
    vs = get_source(self)
    nHackVerts = num_vertices(G)
    degs = [out_degree(v,G) for v in vertices(G)]
    edge_vertices = find(d -> d==2, degs)
    #nonMarkedEdges = [v for (j,v) in enumerate(vertices(G)) if G.degree(v) == 2 and j<nHackVerts-nMarks]
    ret=[]
    for v in edge_vertices
      if v<=nHackVerts-nMarks
        # just permute the newly marked (edge-)vertex to the end of the list of non-marked stuff
        p = collect(1:nHackVerts)
        p[v] = nHackVerts - nMarks  # the last non-marked position
        p[nHackVerts - nMarks] = v
        sgn = get_perm_sign(vs, G, p)
        push!(ret, (permuteGraph(G, p), sgn) )
      end
    end
    return ret
end


type ContractD <: GraphOperator{SmallGraph,SmallGraph}
  nVertices :: Int
  nLoops :: Int
  nMarks :: Int
  evenEdges :: Bool
end

"""Retrieve the file name of the file storing the matrix list for the T graded component.
   Here and below T is the grading of the source."""
function get_file_name(self::ContractD)
  dataDir = self.evenEdges ? markedDataDirEven : markedDataDirOdd
  s = @sprintf "contractD%d_%d_%d.txt" self.nVertices self.nLoops self.nMarks
  return string(dataDir, s)
end

"""The operator maps the T graded piece to the T' graded piece.
   This functions returns T' given T."""
function get_target(self::ContractD)
        return MarkedEdgeGraphVectorSpace(self.nVertices-1, self.nLoops, self.nMarks, self.evenEdges)
end

"""
Returns the GraphVectorSpace on which this operator acts
"""
function get_source(self::ContractD)
  return MarkedEdgeGraphVectorSpace(self.nVertices, self.nLoops, self.nMarks, self.evenEdges)
end

"""For G a graph in HackG format, returns a list of pairs (GG, x),
   such that (contractionOperator)(G) = sum x GG."""
function operate_on(self::ContractD, G)
        nVert,nLoops,nMarks,evenEdges = (self.nVertices, self.nLoops,self.nMarks,self.evenEdges)
        vs = get_source(self)
        nHackVerts = num_vertices(G)
        degs = [out_degree(v,G) for v in vertices(G)]
        edge_vertices = find(d -> d==2, degs)
        #nonMarkedEdges = [v for j,v in enumerate(G.nodes()) if G.degree(v) == 2 and j<nHackVerts-nMarks]
        ret=[]
        for v in edge_vertices
          if v <= nHackVerts-nMarks # it is a non-marked edge
            nb = out_neighbors(v, G)
            a=nb[1]
            b=nb[2]
            # permute the contract-edge-vertex to the first position (1), its neighbors to position 2 and 3
            p = collect(1:nHackVerts)
            p[1] = v
            p[2]=a
            p[3]=b
            idx = 4
            for j =1:nHackVerts
                if j == v || j== a || j== b
                    continue
                else
                    p[idx] = j
                    idx +=1
                end
            end
            pp = invPermutation(p)
            sgn = get_perm_sign(vs, G, pp)
            GG = permuteGraph(G,pp)

            # now delete vertex 0
            remove_edge!(GG, 1,2)
            remove_edge!(GG, 1,3)
            #GG.remove_node(0) -> later
            # ... and join 2 and 3 and fix labelings
            GGG = small_graph(nHackVerts-2)
            for (u,v) in edges(GG)
              uu = (u<=3?1:u-2)
              vv = (v<=3?1:v-2)
              add_edge!(GGG, uu,vv)
            end

            push!(ret, (GGG, sgn))
          end
        end
        return ret
end

function get_work_estimate(self::ContractD)
  # give estimate of number of graphs
  vs = get_source(self)
  nEdges = vs.nLoops + vs.nVertices -1
  return get_work_estimate(vs) * nEdges
end

function dispListCoverageMarked(nDisplay=0)
  data = []
  titles = []
  xHeaders = []
  yHeaders = []
  nVR = collect(3:12)  # vertex range
  nLR=collect(2:8)
  nMR=collect(2:15)

  for evenEdges in [true, false]
      for (k,l) in enumerate(nLR)
        curdata = Array{Any}(length(nVR), length(nMR))
        for (i,v) in enumerate(nVR)
          for (j,m) in enumerate(nMR)
            vs = MarkedEdgeGraphVectorSpace(v,l,m,evenEdges)
            dim = getDimension(vs)
            deco = is_valid(vs) ? "" : "class=redcell"
            curdata[i,j] = Dict("data"=>dim, "style"=>deco)
          end
        end
        push!(data, curdata)
        evo = evenEdges ? "Even " : "Odd "
        push!(titles, "$evo edges, $l loops (vertices \\ marks)")
        push!(xHeaders, nMR)
        push!(yHeaders, nVR)
      end
  end

  dispTables(titles, xHeaders, yHeaders, data, nDisplay=nDisplay)

end

function dispMarkingOperatorCoverageMarked(nDisplay=0)
  data = []
  titles = []
  xHeaders = []
  yHeaders = []
  nVR = collect(3:12)  # vertex range
  nLR=collect(2:8)
  nMR=collect(2:15)

  for evenEdges in [true, false]
      for (k,l) in enumerate(nLR)
        curdata = Array{Any}(length(nVR), length(nMR))
        for (i,v) in enumerate(nVR)
          for (j,m) in enumerate(nMR)
            vs = MarkedEdgeGraphVectorSpace(v,l,m,evenEdges)
            nrEntries = -1
            try
              D = load_matrix(MarkD(v,l,m,evenEdges))
              if D==[]
                nrEntries = 0
              else
                nrEntries = nnz(D)
              end
            catch
            end
            deco = is_valid(vs) ? "" : "class=redcell"

            curdata[i,j] = Dict("data"=>nrEntries, "style"=>deco)
          end
        end
        push!(data, curdata)
        evo = evenEdges ? "Even " : "Odd "
        push!(titles, "$evo edges, $l loops (vertices \\ marks)")
        push!(xHeaders, nMR)
        push!(yHeaders, nVR)
      end
  end

  dispTables(titles, xHeaders, yHeaders, data, nDisplay=nDisplay)

end

function dispContractOperatorCoverageMarked(nDisplay=0)
  data = []
  titles = []
  xHeaders = []
  yHeaders = []
  nVR = collect(3:12)  # vertex range
  nLR=collect(2:8)
  nMR=collect(2:15)

  for evenEdges in [true, false]
      for (k,l) in enumerate(nLR)
        curdata = Array{Any}(length(nVR), length(nMR))
        for (i,v) in enumerate(nVR)
          for (j,m) in enumerate(nMR)
            vs = MarkedEdgeGraphVectorSpace(v,l,m,evenEdges)
            nrEntries = -1
            try
              D = load_matrix(ContractD(v,l,m,evenEdges))
              if D==[]
                nrEntries = 0
              else
                nrEntries = nnz(D)
              end
            catch
            end
            deco = is_valid(vs) ? "" : "class=redcell"

            curdata[i,j] = Dict("data"=>nrEntries, "style"=>deco)
          end
        end
        push!(data, curdata)
        evo = evenEdges ? "Even " : "Odd "
        push!(titles, "$evo edges, $l loops (vertices \\ marks)")
        push!(xHeaders, nMR)
        push!(yHeaders, nVR)
      end
  end

  dispTables(titles, xHeaders, yHeaders, data, nDisplay=nDisplay)

end

function dispContractDCohomologyMarked(nDisplay=0)
  data = []
  titles = []
  xHeaders = []
  yHeaders = []
  nVR = collect(3:6)  # vertex range
  nLR=collect(2:6)
  nMR=collect(2:6)

  for evenEdges in [true, false]
      for (k,l) in enumerate(nLR)
        curdata = Array{Any}(length(nVR), length(nMR))
        for (i,v) in enumerate(nVR)
          for (j,m) in enumerate(nMR)
            vs = MarkedEdgeGraphVectorSpace(v,l,m,evenEdges)
            op = ContractD(v,l,m,evenEdges)
            op2 = ContractD(v+1,l,m,evenEdges)
            println("Computing $v $l $m ...")
            dim = get_cohomology(op, op2)

            deco = is_valid(vs) ? "" : "class=redcell"

            curdata[i,j] = Dict("data"=>dim, "style"=>deco)
          end
        end
        push!(data, curdata)
        evo = evenEdges ? "Even " : "Odd "
        push!(titles, "$evo edges, $l loops (vertices \\ marks)")
        push!(xHeaders, nMR)
        push!(yHeaders, nVR)
      end
  end

  dispTables(titles, xHeaders, yHeaders, data, nDisplay=nDisplay)

end


function dispMarkDCohomologyMarked(nDisplay=0)
  data = []
  titles = []
  xHeaders = []
  yHeaders = []
  nVR = collect(3:6)  # vertex range
  nLR=collect(2:6)
  nMR=collect(2:6)

  for evenEdges in [true, false]
      for (k,l) in enumerate(nLR)
        curdata = Array{Any}(length(nVR), length(nMR))
        for (i,v) in enumerate(nVR)
          for (j,m) in enumerate(nMR)
            vs = MarkedEdgeGraphVectorSpace(v,l,m,evenEdges)
            op = MarkD(v,l,m,evenEdges)
            op2 = MarkD(v,l,m-1,evenEdges)
            println("Computing $v $l $m ...")
            dim = get_cohomology(op, op2)

            deco = is_valid(vs) ? "" : "class=redcell"

            curdata[i,j] = Dict("data"=>dim, "style"=>deco)
          end
        end
        push!(data, curdata)
        evo = evenEdges ? "Even " : "Odd "
        push!(titles, "$evo edges, $l loops (vertices \\ marks)")
        push!(xHeaders, nMR)
        push!(yHeaders, nVR)
      end
  end

  dispTables(titles, xHeaders, yHeaders, data, nDisplay=nDisplay)

end
