using Base.Test

function testSquareZeroAll()
  L = [ContractDOrdinary(v,l,ee) for v in 2:10, l in 2:8, ee in [true, false]]
  L=L[:]
  squareZeroTestGeneric(L)
end

function testSquareZero()
    A = load_matrix(ContractDOrdinary(8,10,false))
    B = load_matrix(ContractDOrdinary(9,10,false))
    #print inducedPerm(perm, lst)
    @test norm(full(B*A)) < 1e-10
    println( "testSquareZero : success")
end

function testSquareZero2()
    createOperatorFile(ContractDOrdinary(7,6,false))
    createOperatorFile(ContractDOrdinary(8,6,false))

    A = load_matrix(ContractDOrdinary(7,6,false))
    B = load_matrix(ContractDOrdinary(8,6,false))
    #print inducedPerm(perm, lst)
    @test norm(full(B*A)) < 1e-10
    println( "testSquareZero : success")
end

function testSquareZero3()
    A = load_matrix(ContractDOrdinary(7,6,true))
    B = load_matrix(ContractDOrdinary(8,6,true))
    #print inducedPerm(perm, lst)
    @test norm(full(B*A)) < 1e-10
    println( "testSquareZero : success")
end

function testPermSignOrdinary()
  G = wheel_graph(5)
  p = [1,3,4,5,6,2]
  @test get_perm_sign(OrdinaryGraphVectorSpace(6,5,false),G,p ) == 1
  @test get_perm_sign(OrdinaryGraphVectorSpace(6,5,true),G,p ) == 1
  p = [1,2,6,5,4,3]
  @test get_perm_sign(OrdinaryGraphVectorSpace(6,5,false),G,p ) == 1
  @test get_perm_sign(OrdinaryGraphVectorSpace(6,5,true),G,p ) == -1
  p = [1,2,3,5,4,6]
  @test get_perm_sign(OrdinaryGraphVectorSpace(6,5,true),G,p ) == 1
  @test get_perm_sign(OrdinaryGraphVectorSpace(6,5,false),G,p ) == -1


  println( "testPermSignOrdinary : success")
end

function testCreateAllOrdinary()
  L = [OrdinaryGraphVectorSpace(v,l,ee) for v in 3:10, l in 3:8, ee in [true, false]]
  L=L[:]
  #createListFiles(L, timeout=20, skipExisting=false)
  createListFiles(L, skipExisting=false)
end

function testCreateAllOrdinaryOps()
  L = [ContractDOrdinary(v,l,ee) for v in 2:10, l in 2:8, ee in [true, false]]
  L=L[:]
  createOperatorFiles(L, skipExisting=false)
end


#testPermSignOrdinary()
testCreateAllOrdinary()
testCreateAllOrdinaryOps()
#testSquareZero()
#testSquareZero3()
#testSquareZero2()
testSquareZeroAll()

#test = OGC.OrdinaryGraphVectorSpace(8,10,true)

#println(GC.get_file_name(test))

#ggg = sc.listG(4,3)

#println(ggg)

#GC.createListFile(test)

#@time GraphComplexes.createListFile(OGC.OrdinaryGraphVectorSpace(8,10,false))
#@time GraphComplexes.createListFile(OGC.OrdinaryGraphVectorSpace(7,10,false))
#@time createListFile(OrdinaryGraphVectorSpace(9,10,false))
#@time createListFile(OrdinaryGraphVectorSpace(9,12,false))

#@time createOperatorFile(ContractDOrdinary(8,10,false))
#@time createOperatorFile(ContractDOrdinary(9,10,false))
#createSvgFiles(OrdinaryGraphVectorSpace(7,5,false))
#createSvgFiles(OrdinaryGraphVectorSpace(6,5,false))
#@time createListFile(OrdinaryGraphVectorSpace(7,5,false))
#@time createListFile(OrdinaryGraphVectorSpace(6,5,false))
#@time createOperatorFile(ContractDOrdinary(7,5,false))
#A=load_matrix(ContractDOrdinary(7,5,false))
#disp_array(A,xVS=OrdinaryGraphVectorSpace(6,5,false),yVS=OrdinaryGraphVectorSpace(7,5,false))

#testSquareZero()
#println(A)
#println(B)
