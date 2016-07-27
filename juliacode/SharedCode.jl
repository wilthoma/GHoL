
#export get_temp_file_name

import Base.+

# python compatibility
len(lst) = length(lst)
+(x::AbstractString, y::AbstractString)=x*y
argsort(lst) = sortperm(lst)

nautypath = "../nauty25r9_hacked/"
shortg = nautypath + "shortg"
countg = nautypath + "countg"
geng = nautypath + "geng"
genbg = nautypath + "genbg"
multig = nautypath + "multig"
labelg = nautypath + "labelg"
hackhack = "autominfo.txt"  # for now, automorphism information is written to this file

displayHtmlFile = "display%d.html"
displayHtmlFilePrefix = "display"



global nTempFile = -1;
function get_temp_file_name()
  global nTempFile
  nTempFile += 1
  return "temp/dummy$nTempFile.txt"
end


"""
  File copy ( cp(...) ) is broken in Julia in that only files <2G can  be copied.
  When the bug is fixed in Julia, this should be replaced by cp() again.
  ToDO: make platform independent. 
"""
function mycp(fromFile, toFile)
  println(`cp $fromFile $toFile`)
  run(`cp $fromFile $toFile`)
end

function parsePerm(s)
    aa = split(s)
    return [parse(Int,a)+1 for a in aa]
end

"""
Calls labelg, and returns a list of lists of permutations.
The last permutation thus returned is the one used to bring the graph to canonical form.
The other permutations are generators of the automorphism group.
"""
function getCanonPermAndAutoms(cFile,FString; outFile="dummyout_canonperm.txt")
    run(`$labelg -qf$FString $cFile $outFile`)
    # currently automs are written to separate file
    #println("Hallo welt")
    lll = readAllLines(hackhack)

    ret=[[] for l in lll]
    #print len(lll)
    for l = lll
        a = split(l,":")
        if len(a) == 2
            b = split(a[2],",")
            #print a[0] + " \n"

            ret[parse(Int,a[1])] = [parsePerm(s) for s in b]
        end
    end
    return ret
end

"""
Computes the inverse of a permutation
:param p: The permutation
"""
function invPermutation(p::Permutation)
    ret = p[:]
    for (i,j) = enumerate(p)
        ret[j]=i
    end
    return ret
end

"""
Computes the sign of a permutation
:param p: The permutation
:return: +1 or -1
"""
function permSign(p::Permutation)
    sgn = 1
    pp=p[:]
    ip = invPermutation(p)
    for i =1:length(p)
        if pp[i] != i
            sgn = sgn * (-1)
            pp[ip[i]] = pp[i]
            ip[pp[i]] = ip[i]
        end
    end
    return sgn
end

"""
Given a subset (list) of numbers lst, and a (larger) permutation p,
compute the induced permutation on elements of lst, i.e., how they change place relative to each other,
ignoring the other elements.
:param p: The (larger) permutation of N elements (numbers)
:param lst: a subset (list) of k<=N elements
:return: the permutation in S_k
"""
function inducedPerm(p::Permutation, lst::Permutation)
    plst = [p[l] for l in lst]
    return argsort(plst)
end

function readAllLines(cFile)
  lll = split(strip(readall(cFile)),"\n")
  filter!(s->s!="", lll)
  return lll
end

"""
creates a list of simple 1vi graphs with at least trivalent vertices
"""
function listG(nVertices, nLoops, onlyonevi=true)
  nEdges = nLoops + nVertices -1
  if (3*nVertices > 2*nEdges) || (nEdges > nVertices*(nVertices-1)/2 )
    # impossible
    return SmallGraph[]
  end
  tempFile = get_temp_file_name()
  run(`$geng $(onlyonevi?"-Cd3":"-cd3") $nVertices $nEdges:$nEdges $tempFile`)
  lll=readAllLines(tempFile)
  return [parse_graph6(l) for l in lll]
end

"""
parses the nauty multig T format to produce a multigraph
"""
function parse_graphT(s)
    a = split(s)
    nEdges = parse(Int,a[2])
    nVert = parse(Int,a[1])
    G = small_multi_graph(nVert)
    for i = 1:nEdges
        u = parse(Int,a[3*i]) +1
        v = parse(Int,a[3*i+1]) +1
        emultiplicity = parse(Int,a[3 * i + 2])
        add_edge!(G,u,v, emultiplicity)
    end
    return G
end

"""
Creates list of all 1-vertex irreducible connected at least trivalent multigraphs,
except for those with just 2 vertices.
"""
function listMultiG(nVert, nLoops, onlyonevi=true)
    # create single edge graphs
    nEdges = nLoops + nVert -1
    tempFile = get_temp_file_name()
    tempFile2 = get_temp_file_name()
    run(`$geng $(onlyonevi?"-Cd1":"-cd1") $nVert 0:$nEdges $tempFile`)
    # the graph with 2 vertices, one edge is missed... we have to add it manually
    if nVert==2
      s2vgraph = to_string(line_graph(2))
      open(tempFile,"a") do f
        write(f,s2vgraph*"\n")
      end
    end
    # make multigraphs
    run(`$multig -T -e$nEdges $tempFile $tempFile2`)
    lll = readAllLines(tempFile2)
    ggg = [parse_graphT(l) for l in lll]
    filter!(g-> (minimum_out_degree(g) >=3), ggg)
    return ggg
end

function listMultiGex(nVert, nLoops)
    # create single edge graphs
    nEdges = nLoops + nVert -1
    tempFile = get_temp_file_name()
    tempFile2 = get_temp_file_name()
    run(`$geng -cd1 $nVert 0:$nEdges $tempFile`)
    # make multigraphs
    run(`$multig -T -e$nEdges $tempFile $tempFile2`)
    lll = readAllLines(tempFile2)
    ggg = [parse_graphT(l) for l in lll]

    return ggg
end


#addprocs(1) # once on program startup to launch a dedicated computational worker
#require("my_computational_funcs.jl")   # load computational functions on all processes

function run_code_with_timeout(computational_func, args, timeout)
println("Starting process....")
  ps = procs()
  if length(ps) < 2
    addprocs(1)
    @everywhere include("startup.jl")
  end
  pindex = procs()[end]
  println("Using process $pindex")

  response = RemoteRef()
  @async put!(response, remotecall_fetch(pindex, computational_func, args...))  # Run computation on worker 2

  start=time()
  while !isready(response) && (time() - start) < timeout   # timeout of 30 seconds
    sleep(0.1)
  end

  if !isready(response)
     interrupt(pindex)      # interrupt the computation on 2
     println("timeout: interrupted call")
     #do_error_processing()
  else
     #do_response_processing(fetch(response))
     println("success")
     println(fetch(response))
  end
end


"""
Returns the first element in list satisfying the criterion f.
If none found, return Void
"""
function myfirst(lst, f)
  for l in lst
    if f(l)
      return l
    end
  end
  return Void
end

mynorm(C) = sum(sum(C.^2))

"""
  vcat, allowing for one or both of the matrices to be zero size
"""
@inline function myvcat(A,B)
  if prod(size(A))==0
    if prod(size(B))==0
      return []
    else
      return B
    end
  else
    if prod(size(B))==0
      return A
    else
      return vcat(A,B)
    end
  end
end

"""
  Produces a list from lst with duplicates removed, using the key function f.
  (I.e., x,y are duplicates if keyfun(x)==keyfun(y))
"""
function myunique(lst, keyfun)
  dd = [keyfun(l) => l for l in lst]
  return [b for (a,b) in dd]
end

@inline function myrank(A)
  if prod(size(A)) >0
    return rank(A)
  else
    return 0
  end
end

"""
  reads a sparse matrix file in plain format
  row col entry
  the max row and column entries determine the matrix size
"""
function read_matrix_file_plain(cFile)
  if filesize(cFile) == 0
    return []  # readdlm cannot process empty files
  end
  A = readdlm(cFile)
  if A==[]
    return []
  else
    # remove 0-index entries (those are to be ignored... this is a legacy thing)
    I = round(Int,A[:,1])
    J = round(Int,A[:,2])
    m=maximum(I)
    n=maximum(J)
    if m*n == 0
      return []
    end
    good =  ((I.>0) & (J .>0))
    return sparse(I[good],J[good], A[good,3],m,n)
  end
end

"""
  reads a sparse matrix file in sms format
  row col entry
"""
function read_matrix_file_sms(cFile)
  data, hdr = readdlm(cFile, header=true)
  m = parse(Int, hdr[1])
  n = parse(Int, hdr[2])
  if hdr[3]=="M"
    return sparse(round(Int,data[:,1]),round(Int,data[:,2]), round(Int, data[:,3]) )
  else
    return sparse(round(Int,data[:,1]),round(Int,data[:,2]), data[:,3] )
  end
end

"""
  writes a sparse matrix in sms file format
"""
function write_matrix_file_sms(A, cFile)
  typechar = eltype(A) <: Integer ? "M" : "R"
  m,n = size(A)
  AA = A.'  # mind that in sms format column index is the quickly changing one
  rows=rowvals(AA)
  vals=nonzeros(AA)
  open(cFile, "w") do f
    write(f, "$m $n $typechar\n")
    for i=1:m
      for j in nzrange(AA,i)
        write(f, "$i $(rows[j]) $(vals[j])\n")
      end
    end
    write(f,"0 0 0")
  end

end

function permutations(n::Int)
  if n==0
    return Vector{Int}[]
  elseif n== 1
    return Vector{Int}[ [1] ]
  else
    ps = permutations(n-1)
    #println(ps)
    ret = Vector{Int}[]
    for j=1:n
      pps = [vcat(p[1:j-1], [n], p[j:n-1])  for p in ps]
      append!(ret, pps )
    end
    return ret
  end
end
