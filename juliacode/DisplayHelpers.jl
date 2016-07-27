


"""
Creates svg files for all graphs.
"""
function createSvgFiles(self::GraphVectorSpace)
        lstFilename = get_file_name(self)
            if isfile(lstFilename)
                imgDir = get_svg_dir(self)
                if !isdir(imgDir)
                    mkpath(imgDir)
                end
                lll = readAllLines(lstFilename)
                ggg = [parse_graph6(s) for s in lll]
                for (idx,G) in enumerate(ggg)
                    cg="__dummy.gv"
                    open(cg,"w") do f
                        write(f,get_dot(self,G))
                    end
                    outFile = "$imgDir$idx.svg"
                    run(`neato -Nshape=point -Tsvg -o$outFile $cg`)
                end
                println( "$(length(ggg)) svgs for $lstFilename created" )
            else
                println( "Cannot create images for inexistent file: $lstFilename" )
            end
end

"""
Creates html file with list of graphs (in dot format already).
"""
function dispGraphListDot(lst; nDisplay=0, captionList=[])
  self = OrdinaryGraphVectorSpace(3,3,false)
  dispFile = "$(displayHtmlFilePrefix)$nDisplay.html"
  captions = (captionList==[] ? (1:length(lst)) : captionList )

      # generate html page
      open(dispFile,"w") do f
          prea="""<html>
          <head>
          <meta http-equiv="refresh" content="10">
          <title>graph list</title>
          <style>
          p {
              text-align:center;
          }
          div {
              display:inline-block;
          }
          img {
              height:200px;max-width:200px;width: expression(this.width > 200 ? 200: true);
          }
          svg {
              height:200px;max-width:200px;width: expression(this.width > 200 ? 200: true);
          }
          </style>
          </head>
          <body>
          """

          write(f,prea)
          for (idx,s) in enumerate(lst)
              #nx.write_graphml(G, dummyfile)
              #call(["graphml2gv", "-o "+dummyfilegv, dummyfile])
              write(f,"<div>" )
              cg="__dummy.gv"
              open(cg,"w") do fff
                  write(fff,s)
              end
              outFile = "__dummy.svg"
              run(`neato -Nshape=point -Tsvg -o$outFile $cg`)
              open(outFile) do ff
                write(f, readall(ff))
              end

              write(f,"<p>$(captions[idx]).</p></div>\n")
          end
          write(f,"</body></html>")

          println( "Display data for written to $dispFile")
        end
end

"""
Creates html file with list of graphs.
"""
function dispG6List(lst; nDisplay=0)
        self = OrdinaryGraphVectorSpace(3,3,false)
        dispGraphListDot([get_dot(self, parse_graph6(s)) for s in lst], nDisplay=nDisplay, captionList=lst)
end

"""
Creates html file with list of graphs.
Assumes that the function render_to_dot is defined for the graphs provided
"""
function dispGraphList(lst; nDisplay=0)
        dispGraphListDot([render_to_dot(g) for g in lst], nDisplay=nDisplay)
end


function dispListFile(self::GraphVectorSpace; nDisplay=0)
        """
        Creates html file with list of graphs.
        """
        inFile = get_file_name(self)
        imgDir = get_svg_dir(self)
        needSvgCreation = false
        dispFile = "$(displayHtmlFilePrefix)$nDisplay.html"
        if !isfile(inFile)
            println( "Cannot display inexistent file: $inFile" )
        else
            lll = readAllLines(inFile)
            ggg = [parse_graph6(s) for s in lll]

            # generate html page
            open(dispFile,"w") do f
                prea="""<html>
                <head>
                <meta http-equiv="refresh" content="10">
                <title>graph list</title>
                <style>
                p {
                    text-align:center;
                }
                div {
                    display:inline-block;
                }
                img {
                    height:200px;max-width:200px;width: expression(this.width > 200 ? 200: true);
                }
                svg {
                    height:200px;max-width:200px;width: expression(this.width > 200 ? 200: true);
                }
                </style>
                </head>
                <body>
                """

                write(f,prea)
                for (idx,G) in enumerate(ggg)
                    #nx.write_graphml(G, dummyfile)
                    #call(["graphml2gv", "-o "+dummyfilegv, dummyfile])
                    write(f,"<div><img src='$imgDir$idx.svg'>" )
                    write(f,"<p>$(idx).</p></div>\n")
                    if !isfile("$imgDir$idx.svg")
                      needSvgCreation = true
                    end
                end
                write(f,"</body></html>")

                println( "Display data for $inFile written to $dispFile")
            end
        end

        if needSvgCreation
          println(".... still need to create svg files ...")
          createSvgFiles(self)
        end
end


"""
Creates a html page containing the matrix A
Optionally, the axes may be labelled by pictures of graphs in the current graph vector space.
:param A: The array (2D-) Array to be displayed
:param xGraphVector:
:param yGraphVector:
:return:
"""
function disp_array(A; xVS=Void, yVS=Void, nDisplay=0)
    dispFile = "$displayHtmlFilePrefix$nDisplay.html"
    open(dispFile,"w") do f
        write(f, """<html>
        <head>
        <meta http-equiv="refresh" content="5">
        <title>graph list</title>
        <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">
        <style>
        p {
            text-align:center;
        }
        div {
            display:inline-block;
        }
        img {
            height:80px;max-width:80px;width: expression(this.width > 80 ? 80: true);
        }
        #datagrid {
          margin:auto;
          padding:30px;
        }
        </style>
        </head>
        <body>
        <div id="datagrid">
        <table class="pure-table">
        <thead><tr>
        <th></th>
        """ )
        m,n = size(A)
        # display header
        if xVS != Void
            imgDir = get_svg_dir(xVS)
            for i = 1:n
                write(f,"<th><img src='$imgDir$i.svg'></th>")
            end
        else
            for i =1:n
                write(f,"<th>$i</th>")
            end
        end
        write(f,"</tr></thead><tbody>")

        if yVS != Void
            imgDir = get_svg_dir(yVS)
        end

        for j =1:m
            write(f, "\n<tr>\n")
            if yVS != Void
                write(f,"<td><img src='$imgDir$j.svg'></td>")
            else
                write(f,"<td>$j</td>")
            end
            write(f,"\n")

            for i =1:n
                write(f,"<td>$(A[j,i])</td>")
            end
            write(f,"\n</tr>\n")
        end

        write(f,"</tbody></table></div></body></html>")
    end

    println( "Matrix displayed in $dispFile")
end


"""
Displays a 3D array in the form of multiple table
"""
function dispTables(tableHeaders, xHeaders, yHeaders, thedata; nDisplay=0)
    dispFile = "$displayHtmlFilePrefix$nDisplay.html"
    open(dispFile,"w") do f
        write(f,"""<html>
        <head>
        <meta http-equiv="refresh" content="5">
        <title>Dimensions Graph Complex</title>
        <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">
        <style>
        p {
            text-align:center;
        }
        div {
            display:inline-block;
        }
        td.redcell {
            background-color : red;
        }
        img {
            height:80px;max-width:80px;width: expression(this.width > 80 ? 80: true);
        }
        </style>
        </head>
        <body>
        """)

        for (i, tblHdr) in enumerate(tableHeaders)
            write(f,
                """
                <h2>$tblHdr</h2>
                <div id="datagrid">
                <table class="pure-table">
                <thead><tr>
                <th> </th>
                """ )
            for Xh in xHeaders[i]
                write(f,"<th>$Xh</th> ")
            end
            write(f,"</tr></thead><tbody>")
            data = thedata[i]
            yHdrs = yHeaders[i]
            for j=1:size(data,1)
                write(f,"<tr><td>$(yHdrs[j])</td>")
                for k=1:size(data,2)
                    dd = data[j,k]
                    txt = dd["data"]
                    style = dd["style"]
                    write(f,"<td $style>$txt</td>")
                end
                write(f,"</tr>\n")
            end
            write(f,"</tbody>\n</table>\n\n")
      end
    end
    println( "Table data displayed in $dispFile")
end
