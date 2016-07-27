
#using Interfaces
#using FileSystem
#using Graphs

import Base.==


#@interface GraphComplex2 begin
#    get_file_name(self::GraphComplex, T::Int)
#    get_svg_dir(self::GraphComplex, T::Int)
#    is_valid_t(self::GraphComplex, T::Int)
#end

# We assume that T mplements the interface CanonizableGraph, however, given
# lack of language support for interfaces this is not enforced
abstract GraphVectorSpace{T}


# this is the type of graph used for talking to nauty
#typealias BaseGraph SimpleGraph

# Interface GraphVectorSpace

function get_file_name{T}(self::GraphVectorSpace{T})
  error("Not implemented")
end

function get_svg_dir{T}(self::GraphVectorSpace{T})
  error("Not implemented")
end

function is_valid{T}(self::GraphVectorSpace{T})
  error("Not implemented")
end

"""Produces a set of graphs whose isomorphism classes span the T component of the vector space.
   (Not necessarily freely!)"""
function get_generating_graphs{T}(self::GraphVectorSpace{T})
  error("Not implemented")
end

"""Internally, each graph is considered colored. This method returns a
   vector indicating the coloring, cf. color_counts_to_nauty_fmt().
"""
function get_color_counts{T}(self::GraphVectorSpace{T})
  error("Not implemented")
end

"""For G a (HackG-) graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
   Here vertex j becomes vertex p[j] in the new graph."""
function get_perm_sign{T}(self::GraphVectorSpace{T}, G::T, p)
  error("Not implemented")
end

"""Converts the graph to a graphviz dot format string.
   This method is used only for visualization, not for computation."""
function get_dot{T}(self::GraphVectorSpace{T}, G::T)
  error("Not implemented")
end

# ----- interface GraphVectorSpace  end -------

=={T}(vs1::GraphVectorSpace{T}, vs2::GraphVectorSpace{T}) = string(vs1)==string(vs2)

"""
  Provides a rough estimate of the amount of work needed to create the graph list.
  (In arbitrary units)
"""
function get_work_estimate{T}(vs::GraphVectorSpace{T})
  return 999999999999999
end


"""
Determines whether the autom list contains odd automorphisms of G
:param vs:
:param G:
:param T:
:param automList:
:return:
"""
function hasOddAutomorphisms{T}(self::GraphVectorSpace{T}, G::T, automList)
  for p in automList
    if get_perm_sign(self, G, p) == -1
      return true
    end
  end

  return false
end

"""
  Converts a list of color counts to the format understood by nauty.
  E.g. [5, 4,1] ->[1,1,1,1,0,1,1,1,0,1]
"""
function color_counts_to_nauty_fmt(lst)
  ret=Int32[]
  for j in lst
    append!(ret, ones(j-1))
    push!(ret,0)
  end
  return ret
end


"""
Creates the list file for a single graph vector space.
:param self: the vector space of type GraphVectorSpace
:return:
"""
function createListFile{T}(self::GraphVectorSpace{T})
        outFile = get_file_name(self)
        outDir = dirname(outFile)
        colorData = get_color_counts(self)

        if !isdir(outDir)
            mkpath(outDir)
        end

        println( string( "Creating File ",outFile,"...") )
        #println("welt")
        L = get_generating_graphs(self)
        #println("hallo")
        # remove duplicates
        println("$(length(L)) graphs generated, removing duplicates and automorphs...")
        graphSet = Set{AbstractString}()
        for G in L
          canonG, isoms = get_canon_and_automorphisms(G,colorData)
          canong6 = to_string(canonG)
          if !in(canong6,graphSet)
            if !hasOddAutomorphisms(self,G, isoms)
              push!(graphSet, canong6)
            end
          end
        end

        println("$(length(graphSet)) graphs survived, writing to file...")
        # write output to file
        open(outFile, "w") do f
            for g6 in graphSet
              write(f, g6+"\n")
            end
        end

        println( "done" )
end

"""
Creates the list file for a single graph vector space.
:param self: the vector space of type GraphVectorSpace
"""
function createListFile_ref{T}(self::GraphVectorSpace{T})
        outFile = get_file_name(self)
        outDir = dirname(outFile)
        if ! isdir(outDir)
            mkpath(outDir)
        end

        println( string( "Creating File ",outFile,"...") )
        #println("welt")
        L = get_generating_graphs(self)
        #println("hallo")
        # remove duplicates
        fString = get_fstring(self)
        #fString = "z"
        tempFile1 = get_temp_file_name()
        tempFile2 = get_temp_file_name()
        tempFile3 = get_temp_file_name()
        open(tempFile1,"w") do f
          for G=L
            write(f,string(generate_graph6(G),"\n"))
          end
        end
        println("Running $shortg -qf$fString $tempFile1 $tempFile2")
        run(`$shortg -qf$fString $tempFile1 $tempFile2`)

        # remove graphs with odd automorphisms
        automList = getCanonPermAndAutoms(tempFile2,fString, outFile=tempFile3);

        sss = readAllLines(tempFile3)

        graphlist = [parse_graph6(s) for s in sss ]
        filter_arr = [ !hasOddAutomorphisms(self,G, automList[i]) for (i,G) in enumerate(graphlist)]
        graphlist = graphlist[filter_arr]

        # write output to file
        open(outFile, "w") do f
            for G in graphlist
              write(f, generate_graph6(G)+"\n")
            end
        end

        println( "done, $(length(graphlist)) graphs generated" )
end



"""
Returns the dimension of the space
Returns -1 if that piece has not been computed yet.
"""
function getDimension{T}(self::GraphVectorSpace{T})
   if !is_valid(self)
       return 0
   end
   fileName = get_file_name(self)
   #print fileName
   if !isfile(fileName)
       return -1 # not yet computed
   end
   return length(readAllLines(fileName))
end


"""
Driver function for the creation of graph list files.
:param lst: A list of GraphVectorSpaces corresponding to the files to be created
:param timeout: Computations taking longer than this (in s) are cancelled. Attention: this might leave a broken file!
:param skipExisting: Whether existing files should be skipped
:return:
"""
function createListFiles(lst; timeout=0, skipExisting=true)
        for vs = lst
            if !is_valid(vs)
                continue
            end
            outFile = get_file_name(vs)
            if skipExisting && isfile(outFile)
                println("Skipping $outFile")
                continue
            end

            if timeout <= 0
                createListFile(vs)
            else
                #println("Running $vs...")
                run_code_with_timeout(createListFile , Any[vs], timeout)
            end
        end
end
