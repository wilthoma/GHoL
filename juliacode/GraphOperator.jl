
# S and T should be the canonizable graph types for source and target vector space
abstract GraphOperator{S, T}


# ------ interface GraphOperator ---------
function get_file_name{S, T}(self::GraphOperator{S,T})
     """Retrieve the file name (and path) of the file storing the matrix."""
     error("Not implemented")
end

function get_unique_file_name{S, T}(self::GraphOperator{S,T})
     """Retrieve a unique file name for the matrix.
        This filename is used when interchanging files with other computers
     """
     error("Not implemented")
end

function get_source{S, T}(self::GraphOperator{S,T})
     """Returns the GraphVectorSpace{S} on which the operator acts."""
     error("Not implemented")
end

function get_target{S, T}(self::GraphOperator{S,T})
          """Returns the GraphVectorSpace{T} in which the operator takes values."""
     error("Not implemented")
end

"""For G::S a graph in the domain, returns a list of pairs (GG, x), GG::T graph
   in the target, x a number,
   such that (operator)(G) = sum x GG."""
function operate_on{S, T}(self::GraphOperator{S,T}, G::S)
     error("Not implemented")
end
# ----- interface end ------


"""
  Provides a rough estimate of the amount of work needed to create the operator file.
  (In arbitrary units)
"""
function get_work_estimate{S, T}(vs::GraphOperator{S,T})
  return 999999999999999
end




function is_valid_op{S, T}(op::GraphOperator{S,T})
  vs = get_source(op)
  tvs = get_target(op)
  return is_valid(vs) && is_valid(tvs)
end


"""
Creates the matrix file that holds the operator.
The corresponding list files for source and target
must exist when calling this function.
"""
function createOperatorFile{S,T}(self::GraphOperator{S,T})

        vs = get_source(self)
        outFile = get_file_name(self)
        inListFile = get_file_name(vs)
        tvs = get_target(self)
        tgtListFile = get_file_name(tvs)

        colorData = get_color_counts(tvs)

        if !isfile(inListFile)
            println( "Cannot create operator file: First create list file $inListFile")
            return
        end
        if !isfile(tgtListFile)
            println("Cannot create operator file: First create list file $tgtListFile")
            return
        end

        println( "Creating File $outFile ..." )

        sss = readAllLines(inListFile)
        lst = [from_string(S, s) for s in sss] # list of source graphs

        src_count = length(lst)
        ll = readAllLines(tgtListFile)
        tgt_count = length(ll)

        println( "List files read ($src_count, $tgt_count graphs)..." )
        #println("Test $(length(sss)) $(length(lst)) __ $(length(ll))")
       # mat = op.computeEdgeMarkingDifferential(ggg, tgtListFile, nMarks, evenEdges)

        count = 0 # counts current index in dummy file
        entries = [[] for G in lst] # will hold matrix

        if src_count == 0 || tgt_count == 0
            # create empty file and return
            open(outFile,"w") do f
            end
            println("Wrote empty file.")
            return
        end

        # lookup g6 -> index in target vector space
        lookup = Dict{String,Int}( s => j for (j,s) in enumerate(ll))

        open(outFile,"w") do f
          for (i,G) in enumerate(lst)
            the_image = operate_on(self, G)
            for (GG,prefactor) in the_image
                # canonize and look up
                GGcan, to_canonical_p = get_canon(GG, colorData)

                GGcang6 = to_string(GGcan)
                #println("$GGcang6 <- $(to_string(GG)): $to_canonical_p  ___ $(invPermutation( to_canonical_p ))")
                if haskey(lookup, GGcang6)
                  sgn = get_perm_sign(tvs, GGcan, to_canonical_p)
                  write(f, "$i $(lookup[GGcang6]) $(sgn * prefactor)\n" )
                end
            end
          end
          # write matrix size
          write(f, "$(length(lst)) $tgt_count 0\n")
          write(f, "$(length(lst)) $tgt_count 0\n")
        end

        println("done")
end

"""
Creates the matrix file that holds the operator.
The corresponding list files for source and target
must exist when calling this function.
"""
function createOperatorFile_ref{S, T}(self::GraphOperator{S,T})

        vs = get_source(self)
        outFile = get_file_name(self)
        inListFile = get_file_name(vs)
        tvs = get_target(self)
        tgtListFile = get_file_name(tvs)

        if !isfile(inListFile)
            println( "Cannot create operator file: First create list file $inListFile")
            return
        end
        if !isfile(tgtListFile)
            println("Cannot create operator file: First create list file $tgtListFile")
            return
        end

        println( "Creating File $outFile ..." )

        sss = readAllLines(inListFile)
        lst = [parse_graph6(s) for s in sss] # list of source graphs

        src_count = length(lst)
        ll = readAllLines(tgtListFile)
        tgt_count = length(ll)
        #println("Test $(length(sss)) $(length(lst)) __ $(length(ll))")
       # mat = op.computeEdgeMarkingDifferential(ggg, tgtListFile, nMarks, evenEdges)

        count = 0 # counts current index in dummy file
        entries = [[] for G in lst] # will hold matrix

        if src_count == 0 || tgt_count == 0
            # create empty file and return
            open(outFile,"w") do f
            end
            return
        end

        # otherwise, compute differential
        tempFile1 = get_temp_file_name()
        tempFile2 = get_temp_file_name()
        tempFile3 = get_temp_file_name()
        tgtFString = get_fstring(tvs) #"z" # get_color_string_from_counts(get_color_counts(tvs))

        open(tempFile1, "w") do f
            for (i,G) = enumerate(lst)
                the_image = operate_on(self, G)
                cur_entries = []
                for (G,prefactor) = the_image
                    write(f, generate_graph6(G)+"\n")
                    count += 1
                    push!(cur_entries,  (count, prefactor) )
                end
                entries[i] = cur_entries
            end
        end

        # canonize target graphs
        perm_lookup = getCanonPermAndAutoms(tempFile1, tgtFString, outFile=tempFile2)
        canonical_g6 = readAllLines(tempFile2) # canonized form of outputs
        canonical_G = [parse_graph6(s) for s in canonical_g6]

        # ... and create lookup table to find the correct indices
        lookup = [ s => j for (j,s) in enumerate(ll)]
#testestlst = readAllLines(tempFile1)
        # fill matrix
        open(outFile,"w") do f
            for i =1:length(lst)
                for p =entries[i]
                    idx, sgn = p
                    #to_canonical_p = invPermutation(perm_lookup[idx][-1]) # last permutation is the one to canonical form
                    to_canonical_p = perm_lookup[idx][end] # last permutation is the one to canonical form
                    can_g6 = canonical_g6[idx]
                    can_G=canonical_G[idx]
                    #println("$can_g6 <- $(testestlst[idx]): $to_canonical_p  ___ ")

                    if haskey(lookup, can_g6)
                        tgt_idx = lookup[can_g6]
                        sgn2 = get_perm_sign(tvs, can_G,to_canonical_p)
                        write(f, "$i $tgt_idx $(sgn * sgn2)\n" )
                    end
                end
            end
            # add a zero entry to fix the size (twice to ensure it is read as 2d array)
            write(f, "$(length(lst)) $tgt_count 0\n")
            write(f, "$(length(lst)) $tgt_count 0\n")
        end

        println("done")
end

function load_matrix{S, T}(self::GraphOperator{S,T})
    inFile = get_file_name(self)
    vs = get_source(self)
    tvs = get_target(self)
    if (!is_valid(vs)) || (!is_valid(tvs))
        return []
    end

    if !isfile(inFile)
        println( "Error: Cannot load matrix: No such file $inFile")
        error("Cannot load file")
    else
        A=[]
        try
          #A = readdlm(inFile)
          A=read_matrix_file_plain(inFile)
        catch y
          # readdlm fails on empty files, corresponding to empty matrices
          println(self)
          println(y)
          return []
        end

        #println(A)
        # unfortunately, scipy cannot handle size zero sparse matrices
        #return sparse(round(Int,A[:,1]),round(Int,A[:,2]),round(Int128,A[:,3])//1)
        return A
    end
end


"""
Driver function for the creation of operator matrix files.
"""
function createOperatorFiles(lst; timeout=0, skipExisting=true)
  for op = lst
      vs = get_source(op)
      tvs = get_target(op)
      if (!is_valid(vs)) || (!is_valid(tvs))
          continue
      end
      outFile = get_file_name(op)
      if skipExisting && isfile(outFile)
          println("Skipping $outFile")
          continue
      end

      if timeout <= 0
          createOperatorFile(op)
      else
          #println("Running $op...")
          run_code_with_timeout(createOperatorFile , Any[op], timeout)
      end
  end
end

"""
  Creates all operator and list files necessary for the given list of operator files.
"""
function createListAndOperatorFiles(oplist; skipExisting= true)
  vslst = vcat( [get_source(op) for op in oplist], [get_target(op) for op in oplist] )
  vslst = myunique(vslst, string)
  createListFiles(vslst,skipExisting=skipExisting)
  createOperatorFiles(oplist,skipExisting=skipExisting)
end


function getDiracOperator{U,S,T}(operatorD::GraphOperator{S,T},operatorDD::GraphOperator{U,S})
  vsD = get_source(operatorD)
  tvsD = get_target(operatorD)
  vsDD = get_source(operatorD)
  tvsDD = get_target(operatorDD) # we assume tvsDD  is equal to vsD !
  Dfile = get_file_name(operatorD)
  DDfile = get_file_name(operatorDD)

  D=[]
  DD=[]
  try
      D = load_matrix(operatorD)
      DD = load_matrix(operatorDD)
  catch
      if returnOnlyDimension
        return -1
      else
        return Void
      end
  end

  A=[]
  if D != [] && DD != []
      A = vcat(D',DD)
  elseif D != [] && DD ==[]
      A = D'
  elseif D == [] && DD != []
      A = DD
  else
      # if both matrices are empty, then the full space is in the kernel
      A=[]
  end


  return A

end

"""
  Computes the cohomology, i.e., ker(D)/im(DD)
"""
function get_cohomology{U,S, T}(operatorD::GraphOperator{S,T},operatorDD::GraphOperator{U,S}, returnOnlyDimension=true)
      vsD = get_source(operatorD)
      tvsD = get_target(operatorD)
      vsDD = get_source(operatorD)
      tvsDD = get_target(operatorDD) # we assume tvsDD  is equal to vsD !
      Dfile = get_file_name(operatorD)
      DDfile = get_file_name(operatorDD)

      D=[]
      DD=[]
      try
          D = load_matrix(operatorD)
          DD = load_matrix(operatorDD)
      catch
          if returnOnlyDimension
            return -1
          else
            return Void
          end
      end

      A=[]
      if D != [] && DD != []
          A = vcat(D',DD)
      elseif D != [] && DD ==[]
          A = D'
      elseif D == [] && DD != []
          A = DD
      else
          # if both matrices are empty, then the full space is in the kernel
          dim = getDimension(vsD)
          if returnOnlyDimension
              return dim
          else
              if dim >=0
                  return eye(dim)
              else
                  return Void
              end
          end
      end


      #NN = nullspace(full(A))
      if returnOnlyDimension
          return sparse_nullspace_dim(A) 
      else
         NN = nullspace(full(A))
          return NN
      end
end


"""
For a given list (or array) of operators, checks whether all the possible compositions are zero
"""
function squareZeroTestGeneric(oplst)
  succ = [] # holds pairs for which test was successful
  fail = [] # failed pairs
  triv = [] # pairs for which test trivially succeeded because at least one operator is the empty matrix
  inc = [] # pairs for which operator matrices are missing
  lst = oplst[:]
  for op1 in lst
    vs1 = get_source(op1)
    for op2 in lst
      tvs2 = get_target(op2)
      if vs1 == tvs2
        # A composable pair is found
        p = (op1, op2)
        if !(is_valid_op(op1) && is_valid_op(op2))
          push!(triv, p)
        else
          D=[]
          DD=[]
          try
            #println(op1)
            D=load_matrix(op1)
            #println("dne1")
            #println(op2)
            DD=load_matrix(op2)
            #println("dne")
          catch
            #println("cannot load ...")
            push!(inc, p)
            continue
          end
          if D==[] || DD==[]
            push!(triv, p)
          else
            if mynorm(DD*D) < 1e-10
              push!(succ, p)
            else
              push!(fail,p)
            end
          end
        end
      end
    end
  end

  # Display results
  println("Square zero test results:")
  println("Success: $(length(succ)), Trivial success $(length(triv)) pairs")
  println("Inconclusive (data not available): $(length(inc))")
  println("Failed: $(length(fail))")
  for (op1,op2) in fail
    println("     $op1 $op2")
  end
end

"""
For a given list (or array) of operators, checks whether all the possible compositions are zero
"""
function commuteTestGeneric(oplst1, oplst2; antiCommute = false)
  succ = [] # holds pairs for which test was successful
  fail = [] # failed pairs
  triv = [] # pairs for which test trivially succeeded because at least one operator is the empty matrix
  inc = [] # pairs for which operator matrices are missing
  lst1 = oplst1[:]
  lst2 = oplst2[:]
  slst1 = [string(get_source(op)) for op in lst1]
  tlst1 = [string(get_target(op)) for op in lst1]

  slst2 = [string(get_source(op)) for op in lst2]
  tlst2 = [string(get_target(op)) for op in lst2]

  println("Commute test: list1 _ $(length(lst1)) entries, list 2 _ $(length(lst2)) entries")
  for (i1a,op1a) in enumerate(lst1)
    for (i2a,op2a) in enumerate(lst2)
      if slst1[i1a] == slst2[i2a]
        for (i1b,op1b) in enumerate(lst1)
          if tlst2[i2a] == slst1[i1b]
            for (i2b,op2b) in enumerate(lst2)
             if tlst1[i1a] == slst2[i2b] && tlst2[i2b] == tlst1[i1b]
                # we found a composable quadruple
                println("quadruple found")
                p = (op1a, op1b, op2a, op2b)

                if !(is_valid_op(op1a) && is_valid_op(op2b)) && !(is_valid_op(op1b) && is_valid_op(op2a))
                  push!(triv, p)
                elseif (is_valid_op(op1a) && is_valid_op(op2b)) && !(is_valid_op(op1b) && is_valid_op(op2a))
                  D=[]
                  DD=[]
                  try
                    D=load_matrix(op1a)
                    DD=load_matrix(op2b)
                  catch
                    #println("cannot load")
                    push!(inc, p)
                    continue
                  end
                  if D==[] || DD==[]
                    push!(triv, p)
                  else
                    if mynorm(D*DD) < 1e-10
                      push!(succ, p)
                    else
                      push!(fail,p)
                    end
                  end
                elseif !(is_valid_op(op1a) && is_valid_op(op2b)) && (is_valid_op(op1b) && is_valid_op(op2a))
                  D=[]
                  DD=[]
                  try
                    D=load_matrix(op2a)
                    DD=load_matrix(op1b)
                  catch
                    #println("cannot load")
                    push!(inc, p)
                    continue
                  end
                  if D==[] || DD==[]
                    push!(triv, p)
                  else
                    if mynorm(D*DD) < 1e-10
                      push!(succ, p)
                    else
                      push!(fail,p)
                    end
                  end
                else
                  D1a=[]
                  D1b=[]
                  D2a=[]
                  D2b=[]
                  try
                    D1a=load_matrix(op1a)
                    D1b=load_matrix(op1b)
                    D2a=load_matrix(op2a)
                    D2b=load_matrix(op2b)
                  catch
                    #println("cannot load")
                    push!(inc, p)
                    continue
                  end
                  if (D1a==[] || D2b==[]) && (D2a==[] || D1b==[])
                    push!(triv, p)
                  elseif !(D1a==[] || D2b==[]) && (D2a==[] || D1b==[])
                    if mynorm(D1a*D2b) < 1e-10
                      push!(succ, p)
                    else
                      push!(fail,p)
                    end
                  elseif (D1a==[] || D2b==[]) && !(D2a==[] || D1b==[])
                    if mynorm(D2a*D1b) < 1e-10
                      push!(succ, p)
                    else
                      push!(fail,p)
                    end
                  else
                    if mynorm(D1a*D2b + (antiCommute?1:-1)*D2a*D1b) < 1e-10
                      push!(succ, p)
                    else
                      push!(fail,p)
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end

  # Display results
  println("Commute test results:")
  println("Success: $(length(succ)), Trivial success $(length(triv)) pairs")
  println("Inconclusive (data not available): $(length(inc))")
  println("Failed: $(length(fail))")
  for (op1,op2,op3,op4) in fail
    println("     $op1 $op2 $op3 $op4")
  end
end

function get_rank_file{S, T}(self::GraphOperator{S,T})
    return get_file_name(self) * ".rank.txt"
end

"""
  readRank
  Tries to read the rank from file <MATRIXFILENAME>.rank if present.
  
  Returns -1 if file is not present.
  Returns 0 if parameters are not in valid range
"""
function readRank{S, T}(self::GraphOperator{S,T})
    if !(is_valid_op(self))
      return 0
    end

    RankFile1 = get_rank_file(self)
    RankFile2 = joinpath(EXCHANGE_DIR_RANK , get_unique_file_name(self) + ".rank")
    if isfile(RankFile1)
       c = readall(RankFile1)
       return parse(Int, c)
    elseif isfile(RankFile2)
       c = readAllLines(RankFile2)
       # Copy the rank file from the exchange dir to the "correct place" in the local filesystem
       cp(RankFile2, RankFile1)
       return parse(Int, c)
    else
       return -1
    end

end

"""
  computeRank
  Computes the rank in julia, returns the result and writes the result to file
  <MATRIXFILENAME>.rank for later use
"""
function computeRank{S, T}(self::GraphOperator{S,T})
    if !(is_valid_op(self))
      return 0
    end

    RankFile = get_rank_file(self)
    r = 0
    try
      A = load_matrix(self)
      if A != []
        r = rank(full(A))
      end
    catch
      # if file cannot be loaded -> indicate failure
      return -1
    end    

    outfile = open(RankFile, "w")
    write(outfile, "$r")
    close(outfile)

    return r

end



"""
  scheduleForRankComputation
  Copies the matrix file to the exchange dir, so that the rank can be computed externally.
  The intended workflow is that the exchange directory is mirrored on the external compute host.
  There the rank is computed, and the rank file created, and mirrored back to our machine.
  Then reading the rank will copy the rank file to the correct directory.
"""
function scheduleForRankComputation{S, T}(self::GraphOperator{S,T}; skipExisting= false)
    if !(is_valid_op(self))
      #println(self)
      return
    end
    ExchangeMatrixFile = joinpath(EXCHANGE_DIR_RANK , get_unique_file_name(self))
    if skipExisting && isfile(ExchangeMatrixFile)
      println("Skipping $ExchangeMatrixFile...")
      return
    end
    matfile = get_file_name(self)
    if !isfile(matfile)
      println(matfile*": file not present, skipping.")
      return
    else
      #println(matfile*": file present")
    end
    A = load_matrix(self)
    if A != []
      # first create folder if necessary
      outDir = dirname(ExchangeMatrixFile)
      if !isdir(outDir)
          mkpath(outDir)
      end
      #write_matrix_file_sms(A, ExchangeMatrixFile)
      write_matrix_file_sms(round(Int,A), ExchangeMatrixFile)
      println("Wrote $ExchangeMatrixFile.")
    else
      # don't schedule empty files, just report 0 rank
      RankFile = get_rank_file(self)
      outfile = open(RankFile, "w")
      write(outfile, "0")
      close(outfile)
      println("Found empty matrix, wrote $RankFile.")
    end
    # TODO: it is inefficient to load matrix and then write again,... 
    # rather the standard format should be .sms and one just copies the file
end

function scheduleForRankComputations(oplist; skipExisting= true)
    for op in oplist
      scheduleForRankComputation(op, skipExisting=skipExisting)
    end
end

function computeRanks(oplist; skipExisting= true)
    for op in oplist
      RankFile = get_rank_file(op)
      if (skipExisting && isfile(RankFile))
        println("Skipping $RankFile...")
      else
        println("Computing rank, output to $RankFile...")
        computeRank(op)
      end
    end
end


"""
  Computes the homology as dim(V) -rank(D) -rank(DD)
  Returns 0 is insufficient information is present (e.g., rank files missing)
"""
function get_cohomology_by_rank{U,S, T}(operatorD::GraphOperator{S,T},operatorDD::GraphOperator{U,S})
  vs = get_source(operatorD)
  if !is_valid(vs) 
    return 0
  end
  d = getDimension(vs)
  if d<0 
    return -1 # list file not yet computed
  end
  rD = readRank(operatorD)
  rDD = readRank(operatorDD)
  if rD<0 || rDD <0
    return -1
  end
  return d - rD - rDD
end