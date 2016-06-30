# Test code for directed graph
using Base.Test

function testNautyDi()
  g = wheel_digraph_outpointing(10)
  gg = testingDi(g)
  #dispGraphList(Any[g, gg])
  @test g.g == gg.g
  println("testNautyDi: Success")
end

function test_display_directed()
  g = wheel_digraph_outpointing(10)
  gg = copy(g)
  add_edge!(gg,2,1)
  add_edge!(gg,3,1)
  add_edge!(gg,4,1)
  add_edge!(gg,5,1)

  dispGraphList(Any[g, gg])
end

function test_tostring_directed()
  g = wheel_digraph_outpointing(30)
  s = to_string(g)
  g2 = from_string(SmallDiGraph, s)
  s2 = to_string(g2)
  @test s == s2
  @test g.g == g2.g
  println("tostring test: Success")
end

function test_automs_directed()
  n = 20 # number of spokes
  g = wheel_digraph_outpointing(20)
  gcan, automs = get_canon_and_automorphisms(g, Void)
  println( automs )

  gcan2, p = get_canon(g, Void)
  #gcan3 = permuteDiGraph(g, invPermutation(p))
  gcan3 = permuteDiGraph(g, p)
  #dispGraphList([g, gcan, gcan2, gcan3],nDisplay=1)
  @test gcan3.g == gcan2.g
  @test gcan3.g == gcan.g

  # add color
  colordata = color_counts_to_nauty_fmt( [1, n] )
  colordata2 = color_counts_to_nauty_fmt( [2, n-1] )

  gcan, automs = get_canon_and_automorphisms(g, colordata)
  @test length(automs) == 1
  gcan, automs = get_canon_and_automorphisms(g, colordata2)
  @test length(automs) == 0

  println("automs test: Success")
end

testNautyDi()
#test_display_directed()
test_tostring_directed()
test_automs_directed()
