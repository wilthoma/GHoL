using Base.Test


function testPermSignMarked()
    # HackG tetrahedron
    G = small_graph(10)
    add_edge!(G,1,5)
    add_edge!(G,1,6)
    add_edge!(G,1,7)
    add_edge!(G,2,5)
    add_edge!(G,2,8)
    add_edge!(G,2,10)
    add_edge!(G,3,6)
    add_edge!(G,3,8)
    add_edge!(G,3,9)
    add_edge!(G,4,7)
    add_edge!(G,4,9)
    add_edge!(G,4,10)
    p = [1,2,4,3,5,7,6,10,9,8]
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,0,false),G,p ) == 1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,3,false),G,p ) == -1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,0,true),G,p ) == 1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,3,true),G,p ) == -1
    p = [1,4,3,2,7,6,5,9,8,10]
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,0,false),G,p ) == 1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,1,false),G,p ) == 1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,0,true),G,p ) == 1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,1,true),G,p ) == 1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,3,true),G,p ) == -1
    @test get_perm_sign(MarkedEdgeGraphVectorSpace(4,3,3,false),G,p ) == -1

    println( "testPermSignOrdinary : success")
end

function testNChooseM()
  a = nChooseMSets(5,3)
  @test length(a)==10
  for i=1:length(a)
    @test length(a[i]) == 5
    @test sum(a[i]) == 3
  end
  println( "testNChooseM : success")
end

function testCreateAllMarked()
  L = [MarkedEdgeGraphVectorSpace(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  L=L[:]
  sort!(L, by= x->get_work_estimate(x))
  createListFiles(L, skipExisting=false)
end

function testCreateAllMarkedMarkingOps()
  L = [MarkD(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  L=L[:]
  sort!(L, by= x->get_work_estimate(x))
  createOperatorFiles(L, skipExisting=false)
end
function testCreateAllMarkedContractOps()
  L = [ContractD(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  L=L[:]
  sort!(L, by= x->get_work_estimate(x))
  createOperatorFiles(L, skipExisting=false)
end

function testSquareZeroMarkedContract()
  L = [ContractD(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  squareZeroTestGeneric(L)
end

function testSquareZeroMarkedMarking()
  L = [MarkD(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  squareZeroTestGeneric(L)
end

function testSquareZeroMarkedSingle()
  createListFiles([MarkedEdgeGraphVectorSpace(v,l,m,e) for v in 6:6, l in 4:4, m in 4:6, e in [true, false]])
  L = [MarkD(v,l,m,e) for v in 6:6, l in 4:4, m in 4:6, e in [true, false]]
  createOperatorFiles(L, skipExisting=false)
  squareZeroTestGeneric(L)
end

function testCommutatorsMarked()
  L1 = [MarkD(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  L2 = [ContractD(v,l,m,e) for v in 3:8, l in 3:8, m in 3:8, e in [true, false]]
  commuteTestGeneric(L1,L2, antiCommute=false)
end

function testTadpoledCreation()
  nTadpoles = 2
  nMarks=4
  LLL = listMultiGex(6, 4-nTadpoles)
  dispGraphList(LLL, nDisplay=0)
  G = LLL[9]
  println(G)
  LL2 = tadpolify(G, nTadpoles)
  dispGraphList(LL2, nDisplay=1)
  println(LL2[1])
  println(LL2[2])
  GG = LL2[1]
  #LL = [listMarkingsOfG(G, nMarks) for G in LL2]
  #LL = vcat(LL...)
  LL = listMarkingsOfG(GG, nMarks)
  dispGraphList(LL, nDisplay=2)
    println(LL[1])
end

testPermSignMarked()
testNChooseM()

testTadpoledCreation()

#testCommutatorsMarked()
#testCreateAllMarked()
#testCreateAllMarkedMarkingOps()
#testCreateAllMarkedContractOps()
#testSquareZeroMarkedMarking()
#testSquareZeroMarkedContract()

#testCreateAllMarked()
#testCreateAllMarkedMarkingOps()
#testCreateAllMarkedContractOps()

#@time createListFile(MarkedEdgeGraphVectorSpace(4,3,3,false))
#@time createListFile(MarkedEdgeGraphVectorSpace(6,5,6,false))
#@time createListFile(MarkedEdgeGraphVectorSpace(6,5,6,true))
#@time createListFile(MarkedEdgeGraphVectorSpace(4,3,3,true))
#@time createListFile(MarkedEdgeGraphVectorSpace(4,3,4,true))
#@time createListFile(MarkedEdgeGraphVectorSpace(4,3,4,false))
#createSvgFiles(MarkedEdgeGraphVectorSpace(4,3,3,true))
#createSvgFiles(MarkedEdgeGraphVectorSpace(4,3,4,true))
#createSvgFiles(MarkedEdgeGraphVectorSpace(4,3,4,false))

#createSvgFiles(MarkedEdgeGraphVectorSpace(4,3,3,false))
#createSvgFiles(MarkedEdgeGraphVectorSpace(6,5,6,false))
#createSvgFiles(MarkedEdgeGraphVectorSpace(6,5,6,true))

#dispListFile(MarkedEdgeGraphVectorSpace(4,3,3,true), nDisplay=0)
#dispListFile(MarkedEdgeGraphVectorSpace(4,3,4,true), nDisplay=1)
#dispListFile(MarkedEdgeGraphVectorSpace(4,3,4,false), nDisplay=2)
#dispListFile(MarkedEdgeGraphVectorSpace(4,3,3,false), nDisplay=0)
#dispListFile(MarkedEdgeGraphVectorSpace(6,5,6,false), nDisplay=1)
#dispListFile(MarkedEdgeGraphVectorSpace(6,5,6,true), nDisplay=2)
