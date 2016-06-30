
ribbonDataDirOdd = "ribbondata/oddedges/"
ribbonDataDirEven = "ribbondata/evenedges/"
imgBaseDir = "img/"

  type RibbonGraphVectorSpace <: GraphVectorSpace{SmallDiGraph}
    nVertices :: Int
    nGenus :: Int
    nPunctures :: Int
    evenEdges :: Bool
  end

  function get_file_name(self::RibbonGraphVectorSpace)
    dataDir = self.evenEdges ? ribbonDataDirEven : ribbonDataDirOdd
    s = @sprintf "gra%d_%d_%d.d6" self.nVertices self.nGenus self.nPunctures
    return string(dataDir, s)
  end

  function get_svg_dir(self::RibbonGraphVectorSpace)
    dataDir = self.evenEdges ? ordinaryDataDirEven : ordinaryDataDirOdd
    s = @sprintf "imgs%d_%d_%d/" self.nVertices self.nGenus self.nPunctures
    return string(dataDir, imgBaseDir, s)
  end


  function get_color_counts(self::RibbonGraphVectorSpace)
          return Void # no coloring
  end

  """Check whether the T graded subspace can in principle be non-empty.
     For each such T, a file is expected. Otherwise the corresponding
     graded component is considered not computed yet."""
  function is_valid(self::RibbonGraphVectorSpace)

          nEdges = 2*self.nGenus + self.nVertices + self.nPunctures-2
          nLoops = nEdges -nVertices +1
          # at least trivalent, and simple
          return  (3*self.nVertices <= 2*nEdges) && self.nVertices > 0 && nLoops >= 0 && nEdges >= 0
  end

"""
  Computes the number of punctures in the ribbon graph provided
"""
function compute_punctures(G::SmallDiGraph)
  #dispGraphList([G])
  m = G.g
  visited = [false for u in vertices(G)]
  ret = 0
  for u in vertices(G)
    if !visited[u]
      ret +=1
      v = u
      while !visited[v]
        visited[v] = true
        # find next vertex by going along double edge, and then single edge
        w,ww = out_neighbors(v,G)
        if m[ww,v]
          w = ww
        end
        v,vv = out_neighbors(w,G)
        if m[v,w]
          v = vv
        end
      end
    end
  end
  return ret
end





"""
  Finds the indices of directed cycles in the HackG graph, i.e., the vertices in the
  ribbon graph.
  Returns an array of lists of indices, starting with the first index, and in lex ordering
"""
function find_vertex_loops(G::SmallDiGraph)
  m = G.g
  visited = [false for u in vertices(G)]
  ret = Any[]
  for u in vertices(G)
    if !visited[u]
      curCycle = Int[]
      v = u
      while !visited[v]
        visited[v] = true
        push!(curCycle,v)
        # find next vertex by going along single edge
        w,ww = out_neighbors(v,G)
        v = m[w,v]?ww:w
      end
      push!(ret, curCycle)
    end
  end
  ret
end



"""
  The input G is a graph whose vertices are organized in groups specified by vranges.
  The function computes all graphs obtained by cyclically connecting the vertices of each group
  (in all possible ways). The resulting graphs are pushed to ret.
"""
  function create_cycles(G::SmallDiGraph, vranges, ret::Vector{SmallDiGraph})#::SmallDiGraph[]
    println(vranges)

    if length(vranges)>0
      vr = vranges[1]
      deg = length(vr)
      for p in permutations(deg-1)
        pp = [p, deg]
        gg = copy(G)
        for i=1:deg
          add_edge!(gg, vr[pp[i]], vr[pp[i==deg?1:i+1]])
        end
        if length(vranges) == 1
          push!(ret, gg)
        else
          create_cycles(gg, vranges[2:end], ret)
        end
      end
    end
  end

  """Produces a set of graphs whose isomorphism classes span the T component of the vector space.
     (Not necessarily freely!)"""
  function get_generating_graphs(self::RibbonGraphVectorSpace)
      # generate List of unmarked graphs
      nEdges = 2 *self.nGenus + self.nVertices + self.nPunctures-2
      nLoops = nEdges -self.nVertices +1
      LL = listMultiG(self.nVertices, nLoops)
      ret = SmallDiGraph[]
      println("$(length(LL)) plain graphs generated, adding cyclic orders of stars")
      for G in LL
        # generate base graph with all (bidirectional) edges
        basegraph = small_digraph(2 *nEdges)
        vals = [out_degree(v,G) for v in vertices(G)]
        counts = [sum(vals[1:j-1]) for j in 1:self.nVertices+1]
        vranges = [counts[j]+1:counts[j+1] for j in 1:self.nVertices]
        degs = [0 for j in 1:self.nVertices]
        for (j, (u,v)) in enumerate(edges(G))
          degs[u] += 1
          degs[v] += 1
          add_edge!(basegraph, counts[u]+degs[u], counts[v]+degs[v])
          add_edge!(basegraph, counts[v]+degs[v], counts[u]+degs[u])
        end

        create_cycles(basegraph, vranges, ret)
      end
      println("$(length(ret)) cyclified graphs generated")
      filter!(x -> compute_punctures(x)==self.nPunctures, ret)
      return ret
  end

  """For G a graph and p a permutation of the vertices, returns the sign induced by the relabelling by p.
     Here vertex j becomes vertex p[j] in the new graph.
     For ribbon graphs we use sign conventions where the punctures are always even.
     The edges and vertices follow the same convention as for the ordinary complexes.
  """
  function get_perm_sign(self::RibbonGraphVectorSpace, G::SmallDiGraph, p)
      #error("Not implemented yet")
      nEdges = 2 *self.nGenus + self.nVertices + self.nPunctures-2
      nLoops = nEdges -nVertices +1
      nVertices = self.nVertices
      if self.evenEdges
          # The sign is (induced sign on vertices) * (induced sign edge orientations)
          pp = collect(1:nVertices)
          G2 = permuteDiGraph(G,p)
          vert_idx = zeros(Int,nVertices)
          for (i,c) in enumerate(find_vertex_loops(G2))
            for cc in c
              vert_idx[cc] = i
            end
          end
          cycles = find_vertex_loops(G)
          for (i,c) in enumerate(cycles)
            pp[i] = vert_idx[p[first(c)]]
          end
          sgn = permSign(pp)
          sgn *= permSign(p) # sign for edges... note that each HackG vertex is one half-edge
          return sgn
      else
          # The sign is (induced sign of the edge permutation)
          # we assume the edges are always lex ordered
          # for the computation we use that edges(G) returns the edges in lex ordering
          pp = collect(1:nEdges)
          G2 = permuteDiGraph(G,p)
          edge_idx2 = -ones(Int64, nVert, nVert)
          for (j, e) = enumerate(biedges(G2))
              u,v = e
              edge_idx2[u,v] = j
              edge_idx2[v,u] = j
          end
          for (j, e) = enumerate(biedges(G))
              u,v = e
              pp[j] = edge_idx2[p[u],p[v]]
          end

          #println(pp)
          return permSign(pp)
      end
  end


  """Converts the graph to a graphviz dot format string.
     This method is used only for visualization, not for computation."""
  function get_dot(self::RibbonGraphVectorSpace, G)
      ret = "digraph {\n"
      for e = edges(G)
          u,v = e
          if !G.g[v,u]
            ret = string(ret, "$u -> $v;\n")
          elseif u<v
            ret = string(ret, "$u ->[dir=none] $v;\n")
          end
      end
      ret = string(ret, "}")
      return ret
  end

  # -----  Contraction operator --------

  type ContractDRibbon <: GraphOperator{SmallDiGraph,SmallDiGraph}
    # source
    nVertices :: Int
    nGenus :: Int
    nPunctures :: Int
    evenEdges :: Bool
  end

  """Retrieve the file name of the file storing the matrix list for the operator."""
  function get_file_name(self::ContractDRibbon)
    dataDir = self.evenEdges ? ribbonDataDirEven : ribbonDataDirOdd
    s = @sprintf "contractD%d_%d_%d.txt" self.nVertices self.nGenus self.nPunctures
    return string(dataDir, s)
  end

  """Returns target graph vector space."""
  function get_target(self::ContractDRibbon)
    return RibbonGraphVectorSpace(self.nVertices-1, self.nGenus, self.nPunctures, self.evenEdges)
  end

  """
  Returns the GraphVectorSpace on which this operator acts
  """
  function get_source(self::ContractDRibbon)
    return RibbonGraphVectorSpace(self.nVertices, self.nGenus, self.nPunctures, self.evenEdges)
  end

  """For G a graph returns a list of pairs (GG, x),
     such that (operator)(G) = sum x GG.
  """
  function operate_on(self::ContractDRibbon, G::SmallDiGraph)
      vs = get_source(self)

      ret=[]
      for (u,v) = biedges(G)
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
          GG = permuteDiGraph(G,pp)

          # now delete the first biedge
          remove_biedge!(GG,1,2)
          a1 = first(out_neighbors(1, GG))
          b1 = first(in_neighbors(1, GG))
          a2 = first(out_neighbors(2, GG))
          b2 = first(in_neighbors(2, GG))
          remove_edge!(GG, 1,a1)
          remove_edge!(GG, b1,1)
          remove_edge!(GG, 2,a2)
          remove_edge!(GG, b2,2)
          add_edge!(GG, b2,a1)
          add_edge!(GG, b1,a2)

          # finally delete the first two vertices
          # remove_vertex(GG, nVert) ... by creating a new graph
          GGG = small_graph(self.nVertices-2)
          for (uu,vv) = edges(GG)
              add_edge!(GGG, uu, vv)
          end

          push!(ret, (GGG, sgn))
      end
      return ret
  end

  function get_work_estimate(self::RibbonGraphVectorSpace)
    # give estimate of number of graphs
    nEdges = self.nLoops + self.nVertices -1
    n = self.nVertices
    return binomial(div(n*(n-1),2), nEdges) / factorial(n)
  end

  function get_work_estimate(self::ContractDRibbon)
    # give estimate of number of graphs
    vs = get_source(self)
    nEdges = vs.nLoops + vs.nVertices -1
    return get_work_estimate(vs) * nEdges
  end

  function dispListCoverageRibbon(nDisplay=0)
    error("Not implemented")
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

  function dispCohomologyRibbon(nDisplay=0)
    error("Not implemented")
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
  function dispOperatorCoverageRibbon(nDisplay=0)
    error("Not implemented")
    data = []
    nVR = collect(3:12)
    nLR=collect(2:15)

    for evenEdges in [true, false]
      curdata = Array{Any}(length(nVR), length(nLR))
      for (i,v) in enumerate(nVR)
        for (j,l) in enumerate(nLR)
          vs = OrdinaryGraphVectorSpace(v,l,evenEdges)
          nrEntries = -1
          try
            D = load_matrix(ContractDOrdinary(v,l,evenEdges))
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
    end

    dispTables(["Even edges (vertices\\loops)", "Odd edges (vertices\\loops)"], Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)

  end
