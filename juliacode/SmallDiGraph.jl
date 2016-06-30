

type SmallDiGraph
  g::BitMatrix
end

function small_digraph(n)
  g = BitMatrix(n,n)
  for i=1:n
    for j=1:n
      g[i,j] = false
    end
  end
  SmallDiGraph(g)
end

function Base.copy(G::SmallDiGraph)
  gg = copy(G.g)
  return SmallDiGraph(gg)
end

function add_edge!(g::SmallDiGraph, u,v)
  g.g[u,v] = true
end

function remove_edge!(g::SmallDiGraph, u,v)
  g.g[u,v] = false
end

num_vertices(g::SmallDiGraph) = size(g.g,1)
vertices(g::SmallDiGraph) = 1:num_vertices(g)
function out_neighbors(v, g::SmallDiGraph)
  ret = Int[]
  for j=1:num_vertices(g)
    if g.g[v,j]
      push!(ret,j)
    end
  end
  ret
end

"""
  Returns edges in lex ordering.
"""
function edges(g::SmallDiGraph)
  ret = Tuple{Int64, Int64}[]
  n = num_vertices(g)
  for i=1:n
    for j=1:n
      if g.g[i,j]
        push!(ret, (i,j))
      end
    end
  end
  return ret
end

"""
  Returns bidirectional edges in lex ordering.
  Tadpoles are ignored
"""
function biedges(g::SmallDiGraph)
  ret = Tuple{Int64, Int64}[]
  n = num_vertices(g)
  for i=1:n-1
    for j=i+1:n
      if g.g[i,j] && g.g[j,i]
        push!(ret, (i,j))
      end
    end
  end
  return ret
end


#function num_edges(g::SmallGraph)
#  return div(sum(sum(g)), 2)
#end

"""
TODO: currently supports only small graphs (<=63 vertices)
"""
function parse_digraph6(s::AbstractString)
  #println(s)
  n = Int(s[1]-63)
  nn = n*n
  nBlocks = ceil(Integer, nn / 6)
  nnn = nBlocks * 6
  bA = Vector{Bool}(nnn)
  for i=1:nnn
    bA[i]=false
  end
  for j=1:nBlocks
    a = Int(s[j+1])-63
    bA[j*6-5] = (a & 32) != 0
    bA[j*6-4] = (a & 16) != 0
    bA[j*6-3] = (a & 8) != 0
    bA[j*6-2] = (a & 4) != 0
    bA[j*6-1] = (a & 2) != 0
    bA[j*6  ] = (a & 1) != 0
  end

  G = small_digraph(n)
  countb = 0
  for j=1:n
    for i=1:n
      countb += 1
      if bA[countb]
        add_edge!(G,i,j)
      end
    end
  end
  #println(bA)
  #println(edges(G))
  G
end

"""
TODO: currently supports only small graphs (<=63 vertices)
"""

function generate_digraph6(G::SmallDiGraph)
      k = 6
      x = 0
      n = num_vertices(G)
      BIAS6 = 63
      p=IOBuffer()
      print(p, Char(63+n))
      for j = 1:n
          for i = 1:n
              x <<= 1;
              if (G.g[i,j])
                x |= 1
              end
              k -=1
              if k == 0
                  print(p, Char( BIAS6 + x))
                  k = 6
                  x = 0
              end
          end
      end

      if k != 6
         print(p, Char( BIAS6 + (x << k)))
      end
      takebuf_string(p)
end



"""
Creates a new simple graph isomorphic to G but with vertex labels permuted according to p.
Edges at vertex j in the old graph get connected to as p[j] in the new.
Non-bijective p are treated accordingly.
:param G: input graph
:param p: vertex permutation
:return: graph with vertex numbers permuted accordingly
"""
function permuteDiGraph(G::SmallDiGraph,p::Permutation)
    GG = small_digraph(num_vertices(G))
    for (u,v) in edges(G)
        add_edge!(GG,p[u],p[v])
    end
    return GG
end



to_string(G::SmallDiGraph) = generate_digraph6(G)
from_string(dummy::Type{SmallDiGraph}, s) = parse_digraph6(s)


function render_to_dot(G::SmallDiGraph)
    ret = "digraph {\n"
    for e = edges(G)
        u,v = e
        if !G.g[v,u]
          ret = string(ret, "$u -> $v;\n")
        elseif u<=v
          ret = string(ret, "$u -> $v [dir=none];\n")
        end
    end
    ret = string(ret, "}")
    return ret
end

# -------- generators ----------
""" Creates a wheel graph with n spokes (i.e., n+1 vertices).
    The central vertex has number 1, the others follow in cyclic order
"""
function wheel_digraph_outpointing(n)
  g = small_digraph(n+1)
  for i=1:n
    add_edge!(g, 1,i+1)
    add_edge!(g, i+1,i==n?2:i+2)
  end
  g
end

function loop_digraph(n)
  g = small_graph(n)
  for i=1:n
    add_edge!(g, i,i==n?1:i+1)
  end
  g
end
