

function test_fastcanon()
        self = OrdinaryGraphVectorSpace(7,7,false)

        #println("welt")
        L = get_generating_graphs(self)
        fString = "z"
        #println("hallo")
        # remove duplicates
        tempFile1 = get_temp_file_name()
        tempFile2 = get_temp_file_name()
        tempFile3 = get_temp_file_name()
        tempFile4 = get_temp_file_name()
        open(tempFile1,"w") do f
          for G=L
            write(f,string(generate_graph6(G),"\n"))
          end
        end
        println("Running $shortg -qf$fString $tempFile1 $tempFile2")
        run(`$shortg -qf$fString $tempFile1 $tempFile2`)

        # remove graphs with odd automorphisms
        automList = getCanonPermAndAutoms(tempFile2,fString, outFile=tempFile3);

        automListXX = getCanonPermAndAutoms(tempFile1,fString, outFile=tempFile4);
        sss = readAllLines(tempFile4)

        graphSet = Set{AbstractString}()
        for (i, s) in enumerate(sss)
          G = parse_graph6(s)
          GG, automs2 = fast_getcanon(L[i], Void, true)
          GG2, automs22 = fast_getcanon(GG, Void, true)
          s1 = generate_graph6(G)
          s2 = generate_graph6(GG)
          s22 = generate_graph6(GG2)
          s3 = generate_graph6(L[i])
          # should be same
          if s1 != s2
            println("Error: canonizations not the same---: $s1 vs $s2, orig: $s3,other $s22")
            println(edges(G))
            println(edges(GG))

          end

            canong6 = generate_graph6(GG)
            if !in(canong6,graphSet)
              if !hasOddAutomorphisms(self,L[i], automs2)
                push!(graphSet, canong6)
                println("$(length(graphSet)) : $canong6")
              end
            end

        end
        println("canon test done, found $(length(graphSet)) graphs")

        sss = readAllLines(tempFile3)

cnt = 0
        for (i, s) in enumerate(sss)
          #println("Graph $i:")
          G = parse_graph6(s)
          automs1 = automList[i]
          GG, automs2 = fast_getcanon(G, Void, true)

          s1 = generate_graph6(G)
          s2 = generate_graph6(GG)
          # should be same
          if s1 != s2
            println("Error: canonizations not the same: $s1 vs $s2")
            println(edges(G))
            println(edges(GG))

          end

          if hasOddAutomorphisms(self,G, automs1) != hasOddAutomorphisms(self, GG, automs2)
            println("Error: Odd automs different")
          end
          if !hasOddAutomorphisms(self, GG, automs2)
            cnt +=1
            println("$cnt : $s2")
          else
            println("        $s2")
          end


          if length(automs1) != length(automs2)
            println("Error: autom length not the same")
            println(automs1)
            println(automs2)
          else
            for (a1,a2) in zip(automs1, automs2)
              if a1 != a2
              println(a1)
              println(a2)
              println("---")
             end
            end
          end


        end

        println( "done" )
end


function testOperatorCreation()
# compares the
  op = ContractDOrdinary(10,9,false)
  outfile = get_file_name(op)
  @time createOperatorFile_ref(op)
  lst1 = readAllLines(outfile)
  @time createOperatorFile(op)
  lst2 = readAllLines(outfile)

  if lst1 != lst2
    println("Error: Mismatch")
  else
    println("Success: Operator files match")
  end

end

function testListFileCreation()
# compares the
  #vs = OrdinaryGraphVectorSpace(10,9,false)
  vs = MarkedEdgeGraphVectorSpace(5,4,4,false)

  outfile = get_file_name(vs)
  @time createListFile_ref(vs)
  lst1 = readAllLines(outfile)
  @time createListFile(vs)
  lst2 = readAllLines(outfile)

  # compare sets of graphs
  is_hit = [l => false for l in lst2]
  for l in lst1
    if is_hit[l]
      error("Duplicate detected")
      return
    else
      is_hit[l]=true # error if key doesn't exist
    end
  end

  if length(lst1) != length(lst2)
    println("Error: Mismatch")
  else
    println("Success: List files match (possibly up to reordering)")
  end

end

function testListFileCreation2()
# compares the
  #vs = OrdinaryGraphVectorSpace(10,9,false)
  vs = MarkedEdgeGraphVectorSpace(5,4,4,false)

  outfile = get_file_name(vs)
  @time createListFile(vs)
  lst1 = readAllLines(outfile)
  lst2 = readAllLines("../pycodes/"*outfile)

  # compare sets of graphs
  is_hit = [l => false for l in lst2]
  for l in lst1
    if is_hit[l]
      error("Duplicate detected")
      return
    else
      is_hit[l]=true # error if key doesn't exist
    end
  end

  if length(lst1) != length(lst2)
    println("Error: Mismatch : length $(length(lst1)) vs  $(length(lst2))")
  else
    println("Success: List files match (possibly up to reordering)")
  end

end

function testOperatorCreation2()
# compares the
  #op = ContractDOrdinary(8,7,false)
  op=ContractD(5,4,4,false)
  outfile = get_file_name(op)
  @time createOperatorFile(op)
  lst2 = readAllLines(outfile)
  @time createOperatorFile_ref(op)
  lst1 = readAllLines(outfile)

  # compare sets of graphs
  for (l,ll) in zip(lst1, lst2)
      if l != ll
        println("Error: $l     vs       $ll")
      else
        println("ok")
      end
  end

end

function testOperatorCreation3()
# compares the
  #op = ContractDOrdinary(8,7,false)
  op=ContractD(5,4,4,false)
  outfile = get_file_name(op)
  @time createOperatorFile(op)
  lst1 = readAllLines(outfile)
  lst2 = readAllLines("../pycodes/"*outfile)


  # compare sets of graphs
  for (l,ll) in zip(lst1, lst2)
      if l != ll
        println("Error: $l     vs       $ll")
      else
        println("ok")
      end
  end

end

#testListFileCreation()
testListFileCreation2()
#testOperatorCreation3()

#test_fastcanon()
