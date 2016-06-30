

typealias SmallMultiGraph Matrix{Int64}

# each matrix entry determines the number of edges
function small_multi_graph(n::Int)
  ret = Matrix{Int64}(n,n)
  ret :: SmallMultiGraph
  for i=1:n
    for j=1:n
      ret[i,j] = 0
    end
  end
  ret
end

function add_edge!(g::SmallMultiGraph, u::Int,v::Int, multiplicity::Int=1)
  newmul = g[u,v] + multiplicity
  g[u,v] = newmul
  g[v,u] = newmul
end

function remove_edge!(g::SmallMultiGraph, u::Int,v::Int, multiplicity::Int=1)
  newmul = g[u,v] - multiplicity
  if newmul<0
    newmul=0
  end
  g[u,v] = newmul
  g[v,u] = newmul
end

num_vertices(g::SmallMultiGraph) = size(g,1)
vertices(g::SmallMultiGraph) = 1:num_vertices(g)

out_degree2(v::Int,g::SmallMultiGraph) = sum(g[:,v])
function out_degree(v::Int,g::SmallMultiGraph)
  ret = 0
  for i=1:num_vertices(g)
      ret +=g[v,i]
  end
  ret += g[v,v] # a tadpole counts as valence 2
  return ret
end
in_degree(v,g) = out_degree(v,g)
out_neighbors(v::Int, g::SmallMultiGraph) = find(u->g[u,v]>0, vertices(g))

"""
Returns the minimum degree of a vertex in the graph g.
"""
function minimum_out_degree(g)
  degs = [out_degree(v,g) for v in vertices(g)]
  return min(degs...)
end

"""
  Returns edges in lex ordering.
"""
function edges(g::SmallMultiGraph)
  ret =Tuple{Int64, Int64}[]
  n = num_vertices(g)
  for i=1:n
    for j=i:n
      for k=1:g[i,j]
        push!(ret, (i,j))
      end
    end
  end
  return ret
end

num_edges(g) = length(edges(g))


# a graph with marked edges
type MarkedSmallMultiGraph
  g :: SmallMultiGraph
  markedg :: SmallMultiGraph
end

marked_small_multi_graph(n) = MarkedSmallMultiGraph(small_multi_graph(n), small_multi_graph(n))

function add_edge!(g::MarkedSmallMultiGraph, u,v; marked=false)
  if marked
    add_edge!(g.markedg,u,v)
  else
    add_edge!(g.g,u,v)
  end
end

add_marked_edge!(g::MarkedSmallMultiGraph, u,v) = add_edge!(g,u,v,marked=true)

edges(g::MarkedSmallMultiGraph) =edges(g.g)  # only the unmarked egdes!!
marked_edges(g::MarkedSmallMultiGraph) =edges(g.markedg)
num_edges(g) = length(edges(g))
num_marked_edges(g::MarkedSmallMultiGraph) = length(marked_edges(g))
num_vertices(g::MarkedSmallMultiGraph) = num_vertices(g.g)

function small_graph_to_multigraph(G::SmallGraph)
  n = num_vertices(G)
  GG = small_multi_graph(n)
  for i=1:n
    for j=1:n
      if G[i,j]
        GG[i,j]=1
      end
    end
  end
  GG
end

function render_to_dot(G::SmallMultiGraph)
  ret = "graph {\n"
  for (u,v) in edges(G)
      ret = string(ret, "$u -- $v;\n")
  end
  ret = string(ret, "}")
  return ret
end

function render_to_dot(G::MarkedSmallMultiGraph)
  ret = "graph {\n"
  for (u,v) in edges(G)
      ret = string(ret, "$u -- $v;\n")
  end
  for (u,v) in marked_edges(G)
      ret = string(ret, "$u -- $v[color=\"black:invis:black\"];\n")
  end
  ret = string(ret, "}")
  #println("  --- $ret")
  #println(G.markedg)
  return ret
end
