function testCreateAllOrdinaryWrapper()
  L = [OrdinaryGraphVectorSpaceWrapper(v,l,ee) for v in 3:11, l in 3:11, ee in [true, false]]
  L=L[:]
  #createListFiles(L, timeout=20, skipExisting=false)
  #createListFiles(L, skipExisting=true)
  for LL in L
    createListFile(L,importOnly=true, skipExisting=true)
  end
end

function testCreateAllOrdinaryOps()
  L = [ContractDOrdinaryWrapper(v,l,ee) for v in 2:10, l in 2:8, ee in [true, false]]
  L=L[:]
  for LL in L
    createOperatorFile(L,importOnly=true, skipExisting=true)
  end
end