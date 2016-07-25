// Produces graph differential
#define USAGE \
"liegradiff [-d#deg] [-e] [-l] n [outfile]\n\n"

#define HELPTEXT \
" Generates The graph differential (edge contraction).\n\
\n\
      n     : the number of trivalent vertices, must be even\n\
     -d#deg : the degree, i.e., number of edge markings.\n\
     -l     : produces a list of graphs, but no differential\n\
     -e     : even edges, modulo graphs with double edges\n"


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
static boolean glListOnly=FALSE;
static int gnmaxvertsofhival = -1;   // #verts of valence > 3 allowed
static boolean glevenedges = FALSE;    // used for computing CS graphs

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
    writeg6(outfile, g, 1, n);
    
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
	}   
    
    ++counter;
}

int PRUNE(graph *g,int n,int maxn)
{
    return 0;
}
int PREPRUNE(graph *g,int n,int maxn)
{
    return 0;
}

int main(int argc, char *argv[])
{
    boolean lbadargs= FALSE, lhasn=FALSE, lhasFile=FALSE;
    int n=0, d=0;
    string cFile;
    char *ccFile;

    // Parse command line arguments
    if (argc<2) lbadargs = TRUE;
    for (int i=1;i<argc && !lbadargs;i++) {
        if (argv[i][0] == 0) continue;
        if (argv[i][0]=='-') {
            if (argv[i][1]=='l') { glListOnly = TRUE; }
            else if (argv[i][1]=='d') { d = atoi(argv[i] + 2); }
            else if (argv[i][1]=='e') { glevenedges = TRUE; }           
            else lbadargs = TRUE; 
        } else {
            if (!lhasn) { n = atoi(argv[i]); lhasn=TRUE; }
            else if (!lhasFile) { cFile = argv[i]; ccFile = argv[i]; lhasFile = TRUE; }
            else { lbadargs = TRUE; }
        }
    } 
    
    // Check for valid argument values (TODO)
    if (!lhasn || n<2 || n>30  ) lbadargs = TRUE;
    if (lbadargs) {
        cerr << USAGE << HELPTEXT;
        return -1;
    }
    
    int nedges = 3*n / 2; // marked edge = 2 edges, not counted here
    if (3*n / 2 < d || (n % 2 == 1))
    {
            cerr << "Impossible combination nr vertices / degree\n";
            return -1;
    }
        
    cerr << "Computing graphs with n=" << n << " and "<< d << " marked edges. ("<<nedges<<" edges)\n";
    if (glevenedges) cerr << "Even edges\n";
    if (lhasFile) cerr << "Data is written to "<< cFile << ".\n";
    if (glListOnly) cerr << "List only (no differential computed).\n";
 

   int geng_argc;
    char *geng_argv[7], cn[40], b[40],cedg[40];
    sprintf(cn, "%d", n+d);  // each marked edge amounts to +1 vertex and +1 edge
    sprintf(cedg, "%d:%d", nedges+d, nedges+d);
    
    strcpy(b,cedg);
//cerr << cedg << endl<< b << endl;
  /* Set up geng argument list.  The 0-th argument is the command name.
   * There must be a NULL at the end. */

    geng_argv[0] = "geng";
    geng_argv[1] = "-Cl"; // change to c for just connected graphs
    geng_argv[2] = "-d3";
    geng_argv[3] = "-D3";
    geng_argv[4] = cn;
    geng_argv[5] = b;
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
