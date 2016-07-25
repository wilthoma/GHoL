// Produces graph differential
#define USAGE \
"gradiff [-d#deg] [-p#v] [-e] [-l] n\n\n"

#define HELPTEXT \
" Counts the number of graphs without odd symmetries.\n\
\n\
      n     : the number of vertices\n\
     -d#deg : the number of edges.\n\
     -c     : only connected graphs\n\
     -e     : even edges, degree=3*vertices -2*edges-3 (as in CS theory), modulo graphs with double edges\n"


#include <iostream>
#include <string>
#include <stdlib.h>

#include "nauty/gtools.h"
#include "tgra.h"

using namespace std;

extern "C" {
extern int GENG_MAIN(int, char**);
void OUTPROC(FILE *, graph *, int n);
int PRUNE(graph *g,int n,int maxn);
int PREPRUNE(graph *g,int n,int maxn); 
}

static unsigned long counter;
static int gnmaxvertsofhival = -1;   // #verts of valence > 3 allowed
static boolean glevenedges = FALSE;    // used for computing CS graphs
static boolean glConnectedOnly=FALSE;
static boolean glListOnly=TRUE;


boolean isAdmissible(graph* g, int n)
{    
    // only gnmaxvertsofhival verticess of valence higher than three
  if (gnmaxvertsofhival <0) return TRUE;
    int vcount =0;
    for (int i=0;i<n;i++) {
        int ecount =0;
        for (int j=0;j<n;j++) 
            if (i!=j && TISELEMENT(g[i],j)) 
                ecount++;
        if (ecount>3) vcount++;
    }
    if (vcount>gnmaxvertsofhival) return FALSE;
    
    return TRUE;
}

void 
OUTPROC(FILE *outfile, graph *g, int n)
{
 /* This will be called for each graph. */
    if (!isAdmissible(g,n)) return;
    if (!isNonzeroGraph(g,n, glevenedges)) return;
    // write graph
/*    writeg6(outfile, g, 1, n);
    
 	setword temp[50];
	int sign;
	// compute differential
	if (!glListOnly) {
	    for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) if (TISELEMENT(g[i],j)) {
	        sign = contract(g,n,i,j, temp, glevenedges);
	        if (sign != 0) {
		        fprintf(outfile,"x %d ", sign);
		        writeg6(outfile, temp, 1, n-1);
	        }
	    }
	}   */
    
    ++counter;
}

int PRUNE(graph *g,int n,int maxn)
{
    if (isAdmissible(g, n))
        return 0;
    else return 1;    
}
int PREPRUNE(graph *g,int n,int maxn)
{
    if (isAdmissible(g, n))
        return 0;
    else return 1;
}

int main(int argc, char *argv[])
{
    boolean lbadargs= FALSE, lhasn=FALSE, lhasFile=FALSE;
    int n=0, d=0;
    string cFile;
    char *ccFile;

    // Parse command line arguments
    int nedges = 0;
    if (argc<2) lbadargs = TRUE;
    for (int i=1;i<argc && !lbadargs;i++) {
        if (argv[i][0] == 0) continue;
        if (argv[i][0]=='-') {
            if (argv[i][1]=='d') { nedges = atoi(argv[i] + 2); }
            else if (argv[i][1]=='c') { glConnectedOnly = true; }
            //else if (argv[i][1]=='p') { gnmaxvertsofhival = atoi(argv[i] + 2); }
            else if (argv[i][1]=='e') { glevenedges = TRUE; }           
            else lbadargs = TRUE; 
        } else {
            if (!lhasn) { n = atoi(argv[i]); lhasn=TRUE; }
            else if (!lhasFile) { cFile = argv[i]; ccFile = argv[i]; lhasFile = TRUE; }
            else { lbadargs = TRUE; }
        }
    } 
    
    // Check for valid argument values (TODO)
    if (!lhasn || n<2 || n>30 || nedges<1) lbadargs = TRUE;
    if (lbadargs) {
        cerr << USAGE << HELPTEXT;
        return -1;
    }
    
    if (FALSE)//(glevenedges)
    {
        nedges = (3*n-d-3)/2;
        if ((3*n-d-3) % 2 != 0)
        {
            cerr << "Impossible combination nr vertices / degree\n";
            return -1;
        }
    }
        
    cerr << "Computing graphs with n=" << n << " vertices and "<<nedges<<" edges\n";
    if (glevenedges) cerr << "Even edges\n";
    //if (lhasFile) cerr << "Data is written to "<< cFile << ".\n";
    //if (glListOnly) cerr << "List only (no differential computed).\n";
    //if (gnmaxvertsofhival>=0) 
    //    cerr << "Pruning all graphs with more than "<< gnmaxvertsofhival<<" vertices of degree higher than three.\n";
    

   int geng_argc;
    char *geng_argv[7], cn[40], b[40],cedg[40];
    sprintf(cn, "%d", n);
    sprintf(cedg, "%d:%d", nedges, nedges);
    
    strcpy(b,cedg);
//cerr << cedg << endl<< b << endl;
  /* Set up geng argument list.  The 0-th argument is the command name.
   * There must be a NULL at the end.  This example is for connected
   * bipartite graphs of order 10 and maximum degree at most 4. */

    geng_argv[0] = "geng";
    if (glConnectedOnly)
    	geng_argv[1] = "-cl";
    else
        geng_argv[1] = "-l";
    geng_argv[2] = "-d3";
    geng_argv[3] = cn;
    geng_argv[4] = b;
    geng_argv[5] = NULL;
    geng_argv[6] = NULL;
    geng_argc = 5;
    if (lhasFile) { geng_argv[5] = ccFile; geng_argc++; }

//cerr << geng_argv[4] << endl;

    counter = 0;
//    cerr << geng_argc << endl;
    for (int j=0;j<geng_argc;j++)
    {
    	cerr << geng_argv[j] << " ";
    }
    cerr << endl;
    
    GENG_MAIN(geng_argc,geng_argv);

    printf("Number of non-zero graphs = %lu.\n",counter);


  return 0;
}
