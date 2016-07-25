// Produces graph differential
#define USAGE \
"gradiff [-d#deg] [-p#v] [-e] [-l] [-o] n [outfile_prefix]\n\n"

#define HELPTEXT \
" Generates the hairy graph differential (edge contraction).\n\
\n\
      n     : the number of vertices (including univalent)\n\
     -d#deg : the degree (#deg=2*n-#edges-2). The default is 0.\n\
     -l     : produces a list of graphs, but no differential\n\
     -p#v   : computes only graphs with no more than #v vertices of valence > 3\n\
     -e     : even edges, degree=3*vertices -2*edges-3 (as in CS theory), modulo graphs with double edges\n\
     -o     : opposite parity of hairs\n\
     -h     : output list of numbers of non-edge vertices\n"

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
static boolean glListNonEdgeVertexNumbers=FALSE;
static boolean glOppositeHairParity=FALSE;
static int gnmaxvertsofhival = -1;   // #verts of valence > 3 allowed
static boolean glevenedges = FALSE;    // used for computing CS graphs

boolean isAdmissible(graph* g, int n)
{    
    // kill graphs with vertices of valence 2
    for (int i=0;i<n;i++) {
        int ecount =0;
        for (int j=0;j<n;j++) 
            if (i!=j && TISELEMENT(g[i],j)) 
                ecount++;
        if (ecount==2) 
            return FALSE;
    }


       // find the univalent vertices
    boolean isUnivalent[50];
    for (int i=0;i<n;i++) {
        int ecount =0;
        for (int j=0;j<n;j++) 
            if (i!=j && TISELEMENT(g[i],j)) 
                ecount++;
        if (ecount == 1)
        {
            isUnivalent[i] = TRUE;
        }
        else isUnivalent[i] = FALSE;
    }
    
    // Kill graphs with two hairs at one vertex.
    for (int i=0;i<n;i++) {
        int haircount =0;
        for (int j=0;j<n;j++) 
            if (i!=j && TISELEMENT(g[i],j) && isUnivalent[j]) 
                haircount++;
        if (haircount > 1)
        {
            return FALSE;
        }
    }

    return TRUE;


    // only gnmaxvertsofhival verticess of valence higher than three
 /* if (gnmaxvertsofhival <0) return TRUE;
    int vcount =0;
    for (int i=0;i<n;i++) {
        int ecount =0;
        for (int j=0;j<n;j++) 
            if (i!=j && TISELEMENT(g[i],j)) 
                ecount++;
        if (ecount>3) vcount++;
    }
    if (vcount>gnmaxvertsofhival) return FALSE;
    
    return TRUE;*/
}

void 
OUTPROC(FILE *outfile, graph *g, int n)
{
 /* This will be called for each graph. */
    if (!isAdmissible(g,n)) return;
    
    graph *gg = g;
    int nn = n;
    setword newg[MAXN] ={0};
    graph gcan[MAXN];
    if (glOppositeHairParity)
    {
         /// replace the graph by one for which each hair is doubled (to a univalent and a bivalent vertex)        
         for (int i=0;i<nn;i++)
         {
            int ecount = 0;
            for (int j=0;j<n;j++)
               if (i != j && TISELEMENT(g[i],j))
               {
                  TADDELEMENT(newg[i], j);
                  ecount++;
               }
            
            // if univalent, add an extra edge to a new vertex
            if (ecount == 1)
            {
               TADDELEMENT(newg[i], nn);
               TADDELEMENT(newg[nn], i);
               nn++;
            }
         
         }
         
         // canonise graph
         if (!makecanon2(newg,gcan,nn,NULL, glevenedges))
             return;
         gg = gcan;
    } 
 
 
 
    if (!isNonzeroGraph(gg,nn, glevenedges)) return;
    
    
    // find the univalent and bivalent vertices -> hairs
    boolean isUnivalent[50];
    int nrUnivalent = 0;
    for (int i=0;i<nn;i++) {
        int ecount =0;
        for (int j=0;j<nn;j++) 
            if (i!=j && TISELEMENT(gg[i],j)) 
                ecount++;
        if (ecount == 1)
        {
            isUnivalent[i] = TRUE;
            nrUnivalent++;
        }
        else if (ecount == 2)
        {
            isUnivalent[i] = TRUE;
        }
        else isUnivalent[i] = FALSE;
    }


   if (!glListNonEdgeVertexNumbers)
   {

       // write graph
       fprintf(outfile,"g %d ", nrUnivalent);
       writeg6(outfile, gg, 1, nn);
       

    	setword temp[50];
    	boolean isHairyVertex[50]={0}; // whether a hair is attached to the vertex
    	
	   int sign;
	   // compute differential
       // contract only edges not incident to univalent vertices
       // if opposite parity one has to check that one does not create graphs with multiple hairs
	   if (!glListOnly) {
	       
	          // find all hairy vertices
	          for (int i=0;i<nn;i++)
	            if (isUnivalent[i])
	              for (int j=0;j<nn;j++) 
	                if (i != j && TISELEMENT(gg[i],j))
	                  isHairyVertex[j] = TRUE;
	       
	
	       for (int i=0;i<nn-1;i++) for (int j=i+1;j<nn;j++) 
               if (TISELEMENT(gg[i],j) && !isUnivalent[i] && !isUnivalent[j] && !(isHairyVertex[i] && isHairyVertex[j]) ) {              
       	        sign = contract(gg,nn,i,j, temp, glevenedges);
       	        if (sign != 0) {
       		        fprintf(outfile,"x %d ", sign);
       		        writeg6(outfile, temp, 1, nn-1);
       	        }
	           }
	   }   
    }
    else
    {
         // just compute number of non-edge vertices
           boolean isHairyVertex[50]={0}; // whether a hair is attached to the vertex
           int nrNonEdgeV=0;
           // find all hairy vertices
	          for (int i=0;i<nn;i++)
	          {
	            if (isUnivalent[i])
	              for (int j=0;j<nn;j++) 
	                if (i != j && TISELEMENT(gg[i],j))
	                  isHairyVertex[j] = TRUE;
             }
             for (int i=0;i<nn;i++) 
               if (!isUnivalent[i])
               {
                 int valence=0;
                 for (int j=0;j<nn;j++)
                   if (i != j && TISELEMENT(gg[i],j))
                     valence++;
                 
                 if (valence >= 4 || (valence == 3 && !isHairyVertex[i]))
                   nrNonEdgeV++;
               }
               
           
         fprintf(outfile,"%d %d\n", nrUnivalent, nrNonEdgeV);
    }
    ++counter;
}

int PRUNE(graph *g,int n,int maxn)
{
    return 0;
    if (isAdmissible(g, n))
        return 0;
    else return 1;    
}
int PREPRUNE(graph *g,int n,int maxn)
{
    return 0;
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
    if (argc<2) lbadargs = TRUE;
    for (int i=1;i<argc && !lbadargs;i++) {
        if (argv[i][0] == 0) continue;
        if (argv[i][0]=='-') {
            if (argv[i][1]=='l') { glListOnly = TRUE; }
            else if (argv[i][1]=='d') { d = atoi(argv[i] + 2); }
            else if (argv[i][1]=='p') { gnmaxvertsofhival = atoi(argv[i] + 2); }
            else if (argv[i][1]=='e') { glevenedges = TRUE; } 
            else if (argv[i][1]=='o') { glOppositeHairParity = TRUE; }
            else if (argv[i][1]=='h') { glListNonEdgeVertexNumbers = TRUE; }           
            else lbadargs = TRUE; 
        } else {
            if (!lhasn) { n = atoi(argv[i]); lhasn=TRUE; }
            else if (!lhasFile) { cFile = argv[i]; ccFile = argv[i]; lhasFile = TRUE; }
            else { lbadargs = TRUE; }
        }
    } 
    
    // Check for valid argument values (TODO)
    if (!lhasn || n<2 || n>30) lbadargs = TRUE;
    if (lbadargs) {
        cerr << USAGE << HELPTEXT;
        return -1;
    }
    
    int nedges = 2*n-d-2;
    if (glevenedges)
    {
        nedges = (3*n-d-3)/2;
        if ((3*n-d-3) % 2 != 0)
        {
            cerr << "Impossible combination nr vertices / degree\n";
            return -1;
        }
    }
        
    cerr << "Computing graphs with n=" << n << " vertices in degree "<< d << ". ("<<nedges<<" edges)\n";
    if (glevenedges) cerr << "Even edges\n";
    if (lhasFile) cerr << "Data is written to "<< cFile << ".\n";
    if (glListOnly) cerr << "List only (no differential computed).\n";
    if (gnmaxvertsofhival>=0) 
        cerr << "Pruning all graphs with more than "<< gnmaxvertsofhival<<" vertices of degree higher than three.\n";
    if (glOppositeHairParity)
        cerr << "Hairs with opposite parity.\n";    
    if (glListNonEdgeVertexNumbers)
        cerr << "Only output list of non-edge vertex numbers.\n";   

   int geng_argc;
    char *geng_argv[7], cn[40], b[40],cedg[40];
    sprintf(cn, "%d", n);
    sprintf(cedg, "%d:%d", nedges, nedges);
    
    strcpy(b,cedg);
//cerr << cedg << endl<< b << endl;
  /* Set up geng argument list.  The 0-th argument is the command name.
   * There must be a NULL at the end. */

    geng_argv[0] = "geng";
    geng_argv[1] = "-cl"; // change to c for just connected graphs
    geng_argv[2] = "-d1";
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
