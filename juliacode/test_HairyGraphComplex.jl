function testCreateAllHairy()
  L = [HairyGraphVectorSpace(v,l,h,ee,eh) for v in 3:12, l in 3:8, h in 1:5, ee in [true, false], eh in [true, false]]
  L=L[:]
  #createListFiles(L, timeout=20, skipExisting=false)
  createListFiles(L, skipExisting=true)
end

function testCreateAllHairyOps()
  L = [ContractDHairy(v,l,h,ee,eh) for v in 3:8, l in 3:8, h in 1:5, ee in [true, false], eh in [true, false]]
  L=L[:]
  createOperatorFiles(L, skipExisting=true)
end

function testSquareZeroAll()
  L = [ContractDHairy(v,l,h,ee,eh) for v in 3:8, l in 3:8, h in 1:5, ee in [true, false], eh in [true, false]]
  L=L[:]
  squareZeroTestGeneric(L)
end

function computeHairyRanks()
  L = [ContractDHairy(v,l,h,ee,eh) for v in 3:8, l in 3:8, h in 1:5, ee in [true, false], eh in [true, false]]
  L=L[:]
  computeRanks(L, skipExisting=true)
end

#testCreateAllHairy()

#testCreateAllHairyOps()
#testSquareZeroAll()

computeHairyRanks()