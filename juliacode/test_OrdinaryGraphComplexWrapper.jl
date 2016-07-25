function testCreateAllOrdinaryWrapper()
  L = [OrdinaryGraphVectorSpaceWrapper(v,l,ee) for v in 3:20, l in 3:15, ee in [true, false]]
  L=L[:]
  #createListFiles(L, timeout=20, skipExisting=false)
  #createListFiles(L, skipExisting=true)
  for LL in L
    createListFile(LL,importOnly=true, skipExisting=true)
  end
end

function testCreateAllOrdinaryWrapper2()
  L = [OrdinaryGraphVectorSpaceWrapper(v,l,ee) for v in 3:17, l in 3:10, ee in [true, false]]
  L=L[:]
  #createListFiles(L, timeout=20, skipExisting=false)
  #createListFiles(L, skipExisting=true)
  for LL in L
    createListFile(LL,importOnly=false, skipExisting=true)
  end
end

function testCreateAllOrdinaryWrapperOps()
  L = [ContractDOrdinaryWrapper(v,l,ee) for v in 2:20, l in 2:15, ee in [true, false]]
  L=L[:]
  for LL in L
    createOperatorFile(LL,importOnly=true, skipExisting=true)
  end
end

#testCreateAllOrdinaryWrapper()
#testCreateAllOrdinaryWrapperOps()

testCreateAllOrdinaryWrapper2()