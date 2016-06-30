# nauty interface

# close and reopen to reload
if isdefined(:the_nauty_lib)
  Libdl.dlclose(the_nauty_lib)
end
the_nauty_lib = Libdl.dlopen("../nauty_wrapper/nautywrap.dylib")


"""
  Computes the canonical form of g, together with the permutation to normal form,
  and optionally with a list of generators of the automorphism group.
  The return value is a pair, (GG, lst), where GG is the graph in normal form and
  lst is a list of permutations, the last permutation is the canonizing permutation,
  the others are automorphism group generators.
  Note: The automorphisms are automorphisms of the parameter graph, not of the canonized graph!!!
"""
function fast_getcanon(g::SmallGraph, colors::Vector{Int32}, computeIsoms::Bool=false)
    gg = copy(g)
    n = num_vertices(g)
    #println(n)
    #println(num_edges(gg))
    buf = Array(Int32,500)
    t = ccall( :canonlabel, Int32, (Ptr{UInt64},Int32,Ptr{Int32}, Ptr{Int32},Int32,Int32), gg.chunks,n, colors!=Void ?colors: C_NULL, buf, computeIsoms,true)
    t :: Int32
    #ps = reshape(buf[1:(t+1)*n], n,t+1)
    #println(num_edges(gg))

    ps = Any[map(Int64, buf[n*(j-1)+1:n*j]+1) for j in 1:t+1]
    return (gg, ps)
end

"""
  Returns a pair (GG, p), where GG is the canonized version of G and
  p is the permutation to canonical form
"""
function get_canon(g::SmallGraph, colordata)
  gg = copy(g)
  n = num_vertices(g)
  #println(n)
  #println(num_edges(gg))
  buf = Array(Int32,500)
  t = ccall( :canonlabel, Int32, (Ptr{UInt64},Int32,Ptr{Int32}, Ptr{Int32},Int32, Int32), gg.chunks,n, colordata!=Void ?colordata: C_NULL, buf, false, true)
  t :: Int32
  #ps = reshape(buf[1:(t+1)*n], n,t+1)
  #println(num_edges(gg))
  #println("And t is $t")
  ps = map(Int64, buf[1:n]+1)
  return (gg, ps)
end

"""
  Returns a pair (GG, lst), where GG is the canonized version of G and
  lst is a list of automorphisms of G
"""
function get_canon_and_automorphisms(g::SmallGraph, colordata)
  gg = copy(g)
  n = num_vertices(g)
  #println(n)
  #println(num_edges(gg))
  buf = Array(Int32,500)
  t = ccall( :canonlabel, Int32, (Ptr{UInt64},Int32,Ptr{Int32}, Ptr{Int32},Int32, Int32), gg.chunks,n, colordata!=Void ?colordata: C_NULL, buf, true, false)
  t :: Int32
  #ps = reshape(buf[1:(t+1)*n], n,t+1)
  #println(num_edges(gg))

  ps = Any[map(Int64, buf[n*(j-1)+1:n*j]+1) for j in 1:t]
  return (gg, ps)
end

function get_canon(g::SmallDiGraph, colordata)
  gg = copy(g)
  n = num_vertices(g)
  #println(n)
  #println(num_edges(gg))
  buf = Array(Int32,500)
  t = ccall( :canonlabeldi, Int32, (Ptr{UInt64},Int32,Ptr{Int32}, Ptr{Int32},Int32, Int32), gg.g.chunks,n, colordata!=Void ?colordata: C_NULL, buf, false, true)
  t :: Int32
  #ps = reshape(buf[1:(t+1)*n], n,t+1)
  #println(num_edges(gg))

  ps = map(Int64, buf[1:n]+1)
  return (gg, ps)
end

"""
  Returns a pair (GG, lst), where GG is the canonized version of G and
  lst is a list of automorphisms of G
"""
function get_canon_and_automorphisms(g::SmallDiGraph, colordata)
  gg = copy(g)
  n = num_vertices(g)
  #println(n)
  #println(num_edges(gg))
  buf = Array(Int32,500)
  t = ccall( :canonlabeldi, Int32, (Ptr{UInt64},Int32,Ptr{Int32}, Ptr{Int32},Int32, Int32), gg.g.chunks,n, colordata!=Void ?colordata: C_NULL, buf, true, false)
  t :: Int32
  #ps = reshape(buf[1:(t+1)*n], n,t+1)
  #println(num_edges(gg))

  ps = Any[map(Int64, buf[n*(j-1)+1:n*j]+1) for j in 1:t]
  return (gg, ps)
end

function testingDi(g)
  gg = copy(g)
  n=num_vertices(g)
  ccall( :testingdi, Int32, (Ptr{UInt64},Int32), gg.g.chunks,n)
  return gg
end
