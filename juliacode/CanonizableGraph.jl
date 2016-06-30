
# The interface used for graphs
# at the moment, this is purely for documentation, and not enforced
# in order to have more flexibility for testing
abstract CanonizableGraph

function to_string(G::CanonizableGraph)
  error("Not implemented")
end

# called like: from_string(SmallGraph, "myg6code")
function from_string(dummy::Type{CanonizableGraph},s::AbstractString)
  error("Not implemented")
end

"""
  Returns a pair (GG, p), where GG is the canonized version of G and
  p is the permutation to canonical form
"""
function get_canon(G::CanonizableGraph, colordata)
  error("Not implemented")
end

"""
  Returns a pair (GG, lst), where GG is the canonized version of G and
  lst is a list of automorphisms of G
"""
function get_canon_and_automorphisms(G::CanonizableGraph, colordata)
  error("Not implemented")
end
