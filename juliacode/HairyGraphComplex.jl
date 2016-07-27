#######################################################################
# Provides a representation of the hairy graph complexes (4 cases)
# For efficiency we work with graphs without multiple edges, and
# without multiple hairs.  
# The tripod and 2-hedgehog graphs are hence not represented and should
# be "added by hand".
#######################################################################


# here the translated copies of files reside, letters OE=odd edges, even hair(-vertices)
hairyDataDirOO = joinpath(DATA_DIR, "hairy/OO/")
hairyDataDirOE = joinpath(DATA_DIR, "hairy/OE/")
hairyDataDirEE = joinpath(DATA_DIR, "hairy/EE/")
hairyDataDirEO = joinpath(DATA_DIR, "hairy/EO/")


imgBaseDir = "img/"

type HairyGraphVectorSpace <: GraphVectorSpace{SmallGraph}
  nVertices :: Int # this excludes hair vertices
  nLoops :: Int
  nHairs :: Int
  evenEdges :: Bool
  evenHairs :: Bool # this refers to the vertex at the end of the hair, not the whole hair
end

function get_dataDir(self::HairyGraphVectorSpace)
    return self.evenEdges ? ( self.evenHairs ? hairyDataDirWrapperEE : hairyDataDirWrapperEO) : ( self.evenHairs ? hairyDataDirWrapperOE : hairyDataDirWrapperOO)
end

function get_file_name(self::HairyGraphVectorSpace)
  dataDir = get_dataDir(self)
  s = @sprintf "gra%d_%d_%d.g6" self.nVertices self.nLoops self.nHairs
  return string(dataDir, s)
end



function get_svg_dir(self::HairyGraphVectorSpace)
  dataDir = get_dataDir(self)
  s = @sprintf "imgs%d_%d_%d/" self.nVertices self.nLoops self.nHairs
  return string(dataDir, imgBaseDir, s)
end


function get_color_counts(self::HairyGraphVectorSpace)
    # obviously, all internal vertices are in color 1, the hair vertices are in color 2
    return color_counts_to_nauty_fmt( [self.nVertices, self.nHairs] )
end

function get_fstring(self::HairyGraphVectorSpace)
    #    error("This routine should not be called, this is only a wrapper")
    return "a"^(self.nVertices)
end

"""Check whether the T graded subspace can in principle be non-empty.
   For each such T, a file is expected. Otherwise the corresponding
   graded component is considered not computed yet."""
function is_valid(self::HairyGraphVectorSpaceWrapper)
        nEdges = self.nLoops + self.nVertices -1
        # at least trivalent
        l = (3*self.nVertices <= 2*nEdges + self.nHairs)
        # all numbers positive
        l = l && self.nVertices > 0 && self.nLoops >= 0 && self.nHairs >=0
        # Can have at most a full graph
        l = l && nEdges <= self.nVertices*(self.nVertices-1)/2
        # can have at most one hair per vertex
        l = l && self.nVertices >= self.nHairs
        return l
end

"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs(self::HairyGraphVectorSpaceWrapper)
    nEdges = self.nLoops + self.nVertices -1 +self.nHairs
    tempFile = get_temp_file_name()
    run(`$genbg -czl -d3:1 -D$(nEdges):1 $(self.nVertices) $(self.nHairs) $nEdges:$nEdges $tempFile`)
    lll=readAllLines(tempFile)
    return [parse_graph6(l) for l in lll]
end

"""For G a graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
   Here vertex j becomes vertex p[j] in the new graph."""
function get_perm_sign(self::HairyGraphVectorSpaceWrapper, G::SmallGraph, p)
    # the sign is the same as the corresponding sign in the 
    # ordinary graph complex, apart from an extra contribution from the hair-vertices
    ovs = OrdinaryGraphVectorSpace(self.nVertices+self.nHairs, self.nLoops, self.evenEdges)
    sign = get_perm_sign(ovs, G, p)

    # compute the extra contribution from hairs if necessary
    if self.evenHairs == self.evenEdges
        hairp = inducedperm(p, [self.nVertices+1 : self.nVertices+self.nHairs])
        sign = sign * permSign(hairp)
    end

    return sign
end

"""Converts the graph to a graphviz dot format string.
   This method is used only for visualization, not for computation."""
function get_dot(self::OrdinaryGraphVectorSpaceWrapper, G)
    # no custom rendering necessary
    return render_to_dot(G) 
end


# -----  Contraction operator --------

type ContractDHairy <: GraphOperator{SmallGraph,SmallGraph}
  # source
  nVertices :: Int
  nLoops :: Int
  nHairs :: Int
  evenEdges :: Bool
  evenHairs :: Bool
end


"""Returns target graph vector space."""
function get_target(self::ContractDHairy)
  return HairyGraphVectorSpace(self.nVertices-1, self.nLoops, self.nHairs, self.evenEdges, self.evenHairs)
end

"""
Returns the GraphVectorSpace on which this operator acts
"""
function get_source(self::ContractDHairy)
  return HairyGraphVectorSpace(self.nVertices, self.nLoops, self.nHairs, self.evenEdges, self.evenHairs)
end

"""Retrieve the file name of the file storing the matrix list for the operator."""
function get_file_name(self::ContractDHairy)
  dataDir = get_dataDir(get_source(self))
  s = @sprintf "contractD%d_%d_%d.txt" self.nVertices self.nLoops self.nHairs
  return joinpath(dataDir, s)
end

function get_unique_file_name(self::ContractDHairy)
  prefix = "hairy_" * (self.evenEdges ? "E":"O") * (self.evenHairs ? "E":"O") * "_"
  s = @sprintf "contractD%d_%d.sms" self.nVertices self.nLoops
  return string(prefix, s)
end


"""For G a graph returns a list of pairs (GG, x),
   such that (operator)(G) = sum x GG.
"""
function operate_on(self::ContractDOrdinaryWrapper, G::SmallGraph)
    vs = get_source(self)

    ret=[]
    for (u,v) = edges(G)
        # only edges not connecting to a hair vertex can be contracted
        if u>self.nVertices || v>self.nVertices
            continue
        end 

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

function get_work_estimate(self::HairyGraphVectorSpace)
  # give estimate of number of graphs
  nEdges = self.nLoops + self.nVertices -1
  n = self.nVertices
  return binomial(div(n*(n-1),2), nEdges) / factorial(n)
end

function get_work_estimate(self::ContractDHairy)
  # give estimate of number of graphs
  vs = get_source(self)
  nEdges = vs.nLoops + vs.nVertices -1
  return get_work_estimate(vs) * nEdges
end





#----- visualization ---------------------
function dispListCoverageHairy(nDisplay=0)
  data = []
  nVR = collect(3:24)
  nLR=collect(2:15)
  nHR=collect(1:15)
  captions=[]

  for evenEdges in [true, false]
   for evenHairs in [true, false]
   for (k,h) in enumerate(nHR) 
    curdata = Array{Any}(length(nVR), length(nLR))
    for (i,v) in enumerate(nVR)
      for (j,l) in enumerate(nLR)        
          vs = HairyGraphVectorSpace(v,l,h,evenEdges, evenHairs)
          dim = getDimension(vs)
          deco = is_valid(vs) ? "" : "class=redcell"
          curdata[i,j] = Dict("data"=>dim, "style"=>deco)
      end
    end
    push!(data, curdata)
    push!(captions, "$h-Hairy"*(evenEdges?"E":"O")*(evenHairs?"E":"O")*" (vertices\\loops)" ) 
    end
   end
  end

  dispTables(captions, Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)

end


"""
  Shows a table of the file size of operator.
  Or -1 if operator could not be loaded.
  The number in brackets is the rank, or -1 if not computed
"""
function dispOperatorCoverageHairy(nDisplay=0)
  data = []
  nVR = collect(3:24)
  nLR=collect(2:15)
  nHR=collect(1:15)
  captions=[]

  for evenEdges in [true, false]
   for evenHairs in [true, false]
   for (k,h) in enumerate(nHR) 
    curdata = Array{Any}(length(nVR), length(nLR))
    for (i,v) in enumerate(nVR)
      for (j,l) in enumerate(nLR)        
          theop=ContractDHairy(v,l,h,evenEdges, evenHairs)
          fsize = -1
          try
            fff = get_file_name(theop)
            if isfile(fff)
              fsize=filesize(fff)
            end 
          catch
          end

        rnk = readRank(theop)        
        deco = is_valid_op(theop) ? "" : "class=redcell"
        curdata[i,j] = Dict("data"=>"$fsize ($rnk)", "style"=>deco)
      end
    end
    push!(data, curdata)
    push!(captions, "ContractD $h-Hairy"*(evenEdges?"E":"O")*(evenHairs?"E":"O")*" (vertices\\loops)" ) 
    end
   end
   end

    dispTables(captions, Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)

end

function dispCohomologyHairy(nDisplay=0)
  data = []
  nVR = collect(3:24)
  nLR=collect(2:15)
  nHR=collect(1:15)
  captions=[]

  for evenEdges in [true, false]
   for evenHairs in [true, false]
     for (k,h) in enumerate(nHR) 
        curdata = Array{Any}(length(nVR), length(nLR))
        for (i,v) in enumerate(nVR)
            for (j,l) in enumerate(nLR)        
                vs = HairyGraphVectorSpace(v,l,h,evenEdges,evenHairs)
                op = ContractDHairy(v,l,h,evenEdges,evenHairs)
                op2 = ContractDHairy(v+1,l,h,evenEdges,evenHairs)

                dim = get_cohomology_by_rank(op, op2)
                deco = is_valid(vs) ? "" : "class=redcell"
                curdata[i,j] = Dict("data"=>dim, "style"=>deco)
            end
        end
        push!(data, curdata)
        push!(captions, "Cohomology $h-Hairy"*(evenEdges?"E":"O")*(evenHairs?"E":"O")*" (vertices\\loops)" )
     end 
    end
   end
  
    dispTables(captions, Any[nLR, nLR], Any[nVR,nVR], data, nDisplay=nDisplay)


end

