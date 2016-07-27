function testCreateAllHairyWrapper()
  L = [HairyGraphVectorSpaceWrapper(v,l,h,ee,eh) for v in 3:15, l in 3:9, h in 1:8, ee in [true, false], eh in [true, false]]
  L=L[:]
  for LL in L
    createListFile(LL, skipExisting=true, importOnly=true)
  end
end

function testCreateAllHairyOpsWrapper()
  L = [ContractDHairyWrapper(v,l,h,ee,eh) for v in 3:15, l in 3:8, h in 1:8, ee in [true, false], eh in [true, false]]
  L=L[:]
  for LL in L
    createOperatorFile(LL, skipExisting=true, importOnly=true)
  end
end

function testSquareZeroAllHairyWrapper()
  L = [ContractDHairyWrapper(v,l,h,ee,eh) for v in 3:15, l in 3:8, h in 1:8, ee in [true, false], eh in [true, false]]
  L=L[:]
  squareZeroTestGeneric(L)
end

function computeHairyRanksWrapper()
  L = [ContractDHairyWrapper(v,l,h,ee,eh) for v in 3:8, l in 3:8, h in 1:5, ee in [true, false], eh in [true, false]]
  L=L[:]
  computeRanks(L, skipExisting=true)
end

#testCreateAllHairyWrapper()

#testCreateAllHairyOpsWrapper()
testSquareZeroAllHairyWrapper()

#computeHairyRanksWrapper()