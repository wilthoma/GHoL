function computeOrdinaryRanks()
  L = [ContractDOrdinary(v,l,ee) for v in 2:8, l in 2:6, ee in [true, false]]
  L=L[:]
  computeRanks(L, skipExisting=true)
end

function scheduleOrdinaryRanks()
  L = [ContractDOrdinary(v,l,ee) for v in 2:8, l in 2:6, ee in [true, false]]
  L=L[:]
  scheduleForRankComputations(L, skipExisting=true)
end

computeOrdinaryRanks()
scheduleOrdinaryRanks()