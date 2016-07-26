function computeOrdinaryRanks()
  L = [ContractDOrdinary(v,l,ee) for v in 2:8, l in 2:6, ee in [true, false]]
  L=L[:]
  computeRanks(L, skipExisting=true)
end

function computeOrdinaryRanksWrapper()
  L = [ContractDOrdinaryWrapper(v,l,ee) for v in 2:20, l in 2:8, ee in [true, false]]
  L=L[:]
  computeRanks(L, skipExisting=true)
end

function scheduleOrdinaryRanks()
  L = [ContractDOrdinary(v,l,ee) for v in 2:8, l in 2:8, ee in [true, false]]
  L=L[:]
  scheduleForRankComputations(L, skipExisting=true)
end

function scheduleAllOrdinaryRanks()
  L = [ContractDOrdinary(v,l,ee) for v in 2:15, l in 2:10, ee in [true, false]]
  L=L[:]
  scheduleForRankComputations(L, skipExisting=true)
end

function scheduleAllOrdinaryRanksWrapper()
  L = [ContractDOrdinaryWrapper(v,l,ee) for v in 2:20, l in 2:11, ee in [true, false]]
  L=L[:]
  scheduleForRankComputations(L, skipExisting=true)
end

#computeOrdinaryRanks()
#scheduleOrdinaryRanks()

#scheduleAllOrdinaryRanksWrapper()
computeOrdinaryRanksWrapper()