dir = "olddata/"

println("...")
for nrv=1:20
  for deg=-15:15
    nrE = div(3*nrv -3-deg,2)
    nLoops=nrE-nrv+1
    cFile = dir*"even/graev$(nrv)_$(deg).txt"
    outFile = dir*"even/gra$(nrv)_$(nLoops).sms"
    if isfile(cFile) && !isfile(outFile)
      siz = stat(cFile).size
      if siz < 2e9
        println("Reading $cFile...")
        try
          A = read_matrix_file_plain(cFile)

          if A != []
            println("Writing $outFile...")
            write_matrix_file_sms(round(Int,A),outFile)
          end
        catch y
          println(y)
          println("Couldn't read file")
        end
      end
    end
  end
end


for nrv=1:20
  for deg=-15:15
    nrE = 2*nrv -2-deg
    nLoops=nrE-nrv+1
    cFile = dir*"odd/gra$(nrv)_$(deg).txt"
    outFile = dir*"odd/gra$(nrv)_$(nLoops).sms"
    if isfile(cFile) && !isfile(outFile)
      siz = stat(cFile).size
      if siz < 2e9
        println("Reading $cFile...")
        try
          A = read_matrix_file_plain(cFile)
          if A != []
            println("Writing $outFile...")
            write_matrix_file_sms(round(Int,A),outFile)
          end
        catch y
          println(y)
          println("Couldn't read file")
        end
      end
    end
  end
end
