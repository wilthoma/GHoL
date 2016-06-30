

typealias SmallGraph BitMatrix

typealias Permutation Union{Vector{Int32}, Vector{Int}}

function small_graph(n)
  ret = BitMatrix(n,n)
  for i=1:n
    for j=1:n
      ret[i,j] = false
    end
  end
  ret
end

function add_edge!(g::SmallGraph, u,v)
  g[u,v] = true
  g[v,u] = true
end

function remove_edge!(g::SmallGraph, u,v)
  g[u,v] = false
  g[v,u] = false
end

num_vertices(g::SmallGraph) = size(g,1)
vertices(g::SmallGraph) = 1:num_vertices(g)
out_degree2(v::Int,g::SmallGraph) = sum(g[:,v])
function out_degree(v::Int,g::SmallGraph)
  ret = 0
  for i=1:num_vertices(g)
    if g[v,i]
      ret +=1
    end
  end
  return ret
end
out_neighbors(v::Int, g::SmallGraph) = find(u->g[u,v], vertices(g))

"""
  Returns edges in lex ordering.
"""
function edges22(g::SmallGraph) # experimental version of edges
  function __edges()
    n = num_vertices(g)
    for i=1:n-1
      for j=i+1:n
        if g[i,j]
          produce( (i,j) )
        end
      end
    end
  end
  return Task(__edges)
end
function edges(g::SmallGraph)

  ret = Tuple{Int64, Int64}[]
  n = num_vertices(g)
  for i=1:n
    for j=i:n
      if g[i,j]
        push!(ret, (i,j))
      end
    end
  end
  return ret
end

num_edges(g) = length(edges(g))

#function num_edges(g::SmallGraph)
#  return div(sum(sum(g)), 2)
#end

"""
TODO: currently supports only small graphs (<=63 vertices)
"""
function parse_graph6(s::AbstractString)
  #println(s)
  n = Int(s[1]-63)
  nn = div(n*(n-1),2)
  nBlocks = ceil(Integer, nn / 6)
  nnn = nBlocks * 6
  bA = BitVector(nnn)
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

  G = small_graph(n)
  countb = 0
  for j=2:n
    for i=1:j-1
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
function generate_graph6_reference(G::SmallGraph)
  #println(edges(G))
  #println(G)
  n = num_vertices(G)
  nn = div(n*(n-1),2)
  nBlocks = ceil(Integer, nn / 6)
  nnn = nBlocks * 6
  bA = BitVector(nnn)
  A = G
  count = 0
  for j=2:n
    for i=1:j-1
      count += 1
      bA[count] = A[i,j]
    end
  end

  ret = Char(63+n)
  for i=1:nBlocks
    ic = 32*bA[6*i-5] + 16* bA[6*i-4]+ 8* bA[6*i-3]+ 4* bA[6*i-2]+ 2* bA[6*i-1] + bA[6*i]
    ret = string(ret, Char(63 + ic ))
  end
  #println(ret)
  ret
end

function generate_graph6(G::SmallGraph)
      k = 6
      x = 0
      n = num_vertices(G)
      BIAS6 = 63
      p=IOBuffer()
      print(p, Char(63+n))
      for j = 2:n
          for i = 1:j-1
              x <<= 1;
              if (G[i,j])
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


@inline function generate_graph6_3(G::SmallGraph)
      k = 6
      x = 0
      n = num_vertices(G)
      BIAS6 = 63
      p = Char(63+n)
      for j = 2:n
          for i = 1:j-1
              x <<= 1;
              #if (G[i,j])
                x |= G[i,j]
              #end
              k -=1
              if k == 0
                  p = string(p, Char( BIAS6 + x))
                  k = 6
                  x = 0
              end
          end
      end

      if k != 6
         p = string(p, Char( BIAS6 + (x << k)))
      end
      p
end

"""
Creates a new simple graph isomorphic to G but with vertex labels permuted according to p.
Edges at vertex j in the old graph get connected to as p[j] in the new.
Non-bijective p are treated accordingly.
:param G: input graph
:param p: vertex permutation
:return: graph with vertex numbers permuted accordingly
"""
function permuteGraph(G::SmallGraph,p::Permutation)
    GG = small_graph(num_vertices(G))
    for (u,v) in edges(G)
        add_edge!(GG,p[u],p[v])
    end
    return GG
end

function source(e)
  u,v = e
  u
end

function target(e)
  u,v = e
  v
end


# -------- generators ----------
""" Creates a wheel graph with n spokes (i.e., n+1 vertices).
    The central vertex has number 1, the others follow in cyclic order
"""
function wheel_graph(n)
  g = small_graph(n+1)
  for i=1:n
    add_edge!(g, 1,i+1)
    add_edge!(g, i+1,i==n?2:i+2)
  end
  g
end

function loop_graph(n)
  g = small_graph(n)
  for i=1:n
    add_edge!(g, i,i==n?1:i+1)
  end
  g
end

function line_graph(n)
  g = small_graph(n)
  for i=1:n-1
    add_edge!(g, i,i+1)
  end
  g
end

""" Creates a bipartite complete graph with 2n vertices, such that
    the first n are connected to the second n.
"""
function bipartite_complete(n)
  g = small_graph(2n)
  for i=1:n
    for j=1:n
      add_edge!(g,i,n+j)
    end
  end
end

to_string(G::SmallGraph) = generate_graph6(G)
from_string(dummy::Type{SmallGraph}, s) = parse_graph6(s)
