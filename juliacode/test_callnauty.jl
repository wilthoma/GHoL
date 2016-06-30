
thelib=Libdl.dlopen("../nauty_wrapper/nautywrap.dylib")

function fast_generate_g6(g::SmallGraph)
  n = num_vertices(g)
  #ggg = BitArray(n,64*n)
  #for i=1:size(ggg,1), j=1:size(ggg,2)
#    ggg[i,j]=false
#  end
  #a = g.chunks
  t = ccall( :generateg6, Ptr{UInt8}, (Ptr{UInt64},Int32), g.chunks,n)
  #t = ccall( (:generateg6,"../nauty_wrapper/nautywrap.dylib"), Ptr{UInt8}, (Ptr{UInt64},Int32), g.chunks,n)
  return(bytestring(t))
end

function fast_getcanon2(g::SmallGraph)
    gg = copy(g)
    n = num_vertices(g)
    buf = Array(Int32,500)
    t = ccall( :canonlabel, Int32, (Ptr{UInt64},Int32,Ptr{Int32},Int32), gg.chunks,n, buf, true)
    #ps = reshape(buf[1:(t+1)*n], n,t+1)
    ps = Any[buf[n*(j-1)+1:n*j] for j in 1:t+1]
    #println(t)
    #println(ps')
    return (gg, ps)
end

#g = loop_graph(7)
g = wheel_graph(30)

ggg,ps = fast_getcanon(g)
println(length(ps))
println(fast_generate_g6(g))
println(fast_generate_g6(ggg))

#pp = getCanonPermAndAutoms( [g, g])
#println(pp)

function call_100ktimes(f, arg)
  for i=1:100000
    f(arg)
  end
end
a = ccall( :test_plus, Int32, (Int32,Int32), 5,7)
println(a)


n = num_vertices(g)
a = g.chunks

aa=copy(a)
ccall( :testing, Int32, (Ptr{UInt64},Int32), aa,n)
t = ccall( :generateg6, Ptr{UInt8}, (Ptr{UInt64},Int32), a,n)

println(bytestring(t))
println(generate_graph6(g))



println(fast_generate_g6(g))
println(generate_graph6(g))
#println(generate_graph6_2(g))
 #=
@time call_100ktimes(fast_generate_g6,g)
@time call_100ktimes(generate_graph6,g)
@time call_100ktimes(generate_graph6_reference,g)
@time call_100ktimes(generate_graph6_3,g)
=#

Libdl.dlclose(thelib)

#define DEFAULTOPTIONS_GRAPH(options) optionblk options = \
 #{0,FALSE,FALSE,FALSE,TRUE,FALSE,CONSOLWIDTH, \
#  NULL,NULL,NULL,NULL,NULL,NULL,100,0,1,0,&dispatch_graph,FALSE,NULL}

#=typedef struct optionstruct
{
    int getcanon;             /* make canong and canonlab? */
#define LABELONLY 2   /* new value UNIMPLEMENTED */
    boolean digraph;          /* multiple edges or loops? */
    boolean writeautoms;      /* write automorphisms? */
    boolean writemarkers;     /* write stats on pts fixed, etc.? */
    boolean defaultptn;       /* set lab,ptn,active for single cell? */
    boolean cartesian;        /* use cartesian rep for writing automs? */
    int linelength;           /* max chars/line (excl. '\n') for output */
    FILE *outfile;            /* file for output, if any */
    void (*userrefproc)       /* replacement for usual refine procedure */
         (graph*,int*,int*,int,int*,int*,set*,int*,int,int);
    void (*userautomproc)     /* procedure called for each automorphism */
         (int,int*,int*,int,int,int);
    void (*userlevelproc)     /* procedure called for each level */
         (int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
    void (*usernodeproc)      /* procedure called for each node */
         (graph*,int*,int*,int,int,int,int,int,int);
    void (*invarproc)         /* procedure to compute vertex-invariant */
         (graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
    int tc_level;             /* max level for smart target cell choosing */
    int mininvarlevel;        /* min level for invariant computation */
    int maxinvarlevel;        /* max level for invariant computation */
    int invararg;             /* value passed to (*invarproc)() */
    dispatchvec *dispatch;    /* vector of object-specific routines */
    boolean schreier;         /* use random schreier method */
    void *extra_options;      /* arbitrary extra options */
    =#
