#define ONE_WORD_SETS
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <stdio.h>
#include "nauty/gtools.h"
#include "nauty/nauty.h"
#include "nauty/naugroup.h"

#include "tgra.h"

using namespace std;
#define MMAXN 200

extern "C" {
void
userautomproct(int , permutation *, int *,  int , int , int );
void countsymm(permutation *, int);
}
static boolean glIsEven;
static boolean glevenedges;
static graph gtheg[MAXN];
static int ggrcount;

double& GETVAL(graphvect& v, const string gra)
{
	if (v.count(gra)<1) v[gra]=0;
	return v[gra];
}
int& GETVAL(graphmap& v, const string gra)
{
	if (v.count(gra)<1) v[gra]=0;
	return v[gra];
}

int nredges(graph *g, int n)
{
	int nre = 0;
	for (int i=0;i<n-1;i++) 
		for (int j=i+1;j<n;j++) 
			if (TISELEMENT(g[i],j)) 
				nre++;
	return nre;
}
void printperm(permutation* p, int n)
{
	int i;
	fprintf(stderr, "perm(");
	for(i=0; i<n;i++) fprintf(stderr, "%d", p[i]); fprintf(stderr, " ;%d)\n",n);
}
void invperm(permutation* perm, int n)
{
	int temp[MMAXN], i;
	for (i=0;i<n;i++) temp[i] = perm[i];
	for (i=0;i<n;i++) perm[temp[i]] = i;
}
boolean isperm(permutation* p, int n)
{
	boolean vistd[MMAXN];
	for (int i=0;i<n;i++) vistd[i]=FALSE;
	for (int i=0;i<n;i++) {
		if (p[i]>=n || p[i]<0) return FALSE;
		if (vistd[p[i]]) return FALSE;
		vistd[p[i]] = TRUE;
	}
	return TRUE;
}
int getpermsign(permutation* perm, int n)
{
	if (!isperm(perm, n)) {
		fprintf(stderr, "not a permutation:");
		printperm(perm,n);
		return 0;
	}
	boolean vistd[MMAXN]={0};
	for (int i=0;i<MMAXN;i++) vistd[i]=FALSE;
	int res = 1;
	for (int i=0;i<n;i++) if (!vistd[i]) {
	  res = -res;
	  int j=i;
	  do {
	    vistd[j] = TRUE;
	    j = perm[j];
	    res=-res;
	  } while (j!= i);
	}
	return res;
}
int getgpermsign(graph* g, int n, permutation* perm)
// g is the new graph, and perm[i] the label in g of vertex i of the old one
{
	if (!isperm(perm, n)) {
		fprintf(stderr, "not a permutation(ggps):");
		printperm(perm,n);
		return 0;
	}
	
	int inds[MAXN][MAXN] = {0}, nre=0, perm2[MMAXN];
	// loop through edges of the original graph
	for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) if (TISELEMENT(g[perm[i]],perm[j])) 
	{
	  inds[perm[i]][perm[j]] = nre;
	  inds[perm[j]][perm[i]] = nre++;
	}

	// loop through edges of the new graph (i.e., g)
	nre=0;
	for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) if (TISELEMENT(g[i],j)) {
	  perm2[nre++] = inds[i][j];
	}
//printperm(perm2,nre);
	return getpermsign(perm2, nre);
}

int getgpermsign_even(graph* g, int n, permutation* perm) // for even-edges-case
{
	if (!isperm(perm, n)) {
		fprintf(stderr, "not a permutation(ggps):");
		printperm(perm,n);
		return 0;
	}
	
	int sign_from_edges = 1;
	
	// loop through edges of the original graph
	for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) if (TISELEMENT(g[perm[i]],perm[j])) 
	{
	  // record sign change if edge orientation changed 
	  // (edge is oriented from the vertex with lower to that wth higher label)
	  if (perm[i] > perm[j])
	    sign_from_edges = - sign_from_edges;
	}
	
	// total sign is sign from edges * sign from vertices
	return getpermsign(perm, n) * sign_from_edges;
}

int graphperm(graph* g, int n, permutation* perm, graph* gout)
// perm[i] = new label of vertex i of g
{	
	int i,j;
	for (i=0;i<n; i++) gout[i] = 0;
	for (i=0;i<n-1;i++) for (j=i+1;j<n;j++) if (TISELEMENT(g[i],j)) {
	  TADDELEMENT(gout[perm[i]], perm[j]);
	  TADDELEMENT(gout[perm[j]], perm[i]);
	}

	return getgpermsign(gout, n, perm);
}


void
userautomproct(int count, permutation *p, int *orbits,
              int numorbits, int stabvertex, int n)
/* called by nauty;  operates on glIsEven */
{
//fprintf(stderr, "ain");
	// determine whether permutation is even
	int temp[MAXN],i,j;
	graph temp2[MAXN], *g=temp2;
	 
	for (i=0;i<n;i++) {temp[i] = p[i]; temp2[i]=gtheg[i]; }
	//invperm(temp,n);
//fprintf(stderr, "ain2");
//printperm(temp,n);
    int sign = 1;
    if (glevenedges)
        sign = getgpermsign_even(g, n, temp);
    else 
	    sign = getgpermsign(g, n, temp); //TODO:check if not invperm necessary
	    
	glIsEven = glIsEven && (sign==1);
//fprintf(stderr, "aout");
}

boolean
makecanon2(graph *g, graph *gcan, int n, int* perm, boolean evenedges)
/* gcan := canonise(g) */
{
    int lab[MAXN],ptn[MAXN],orbits[MAXN],i,j,r;
  //  set *gj, rr; /////
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk nauty_stats;
    setword workspace[50];

    options.getcanon = TRUE;
    options.userautomproc = userautomproct;  
        
    for (i=0;i<MAXN;i++) { lab[i]=0; ptn[i]=0; }
        
    for (i=0;i<n;i++) { 
        gtheg[i] = g[i]; 
    }

	glIsEven = TRUE;
	glevenedges = evenedges;
        nauty(g,lab,ptn,NULL,orbits,&options,&nauty_stats,
              workspace,50,1,n,gcan);

    if (perm) for (i=0;i<n;i++) perm[i]=lab[i];
    return glIsEven;
}

boolean isNonzeroGraph(graph *g, int n, boolean evenedges)
{
    graph gcan[MAXN];
    return makecanon2(g,gcan,n,NULL, evenedges);
}

boolean isNonzeroGraphOppositeHairParity(graph *g, int n, boolean evenedges)
{
    // find hairs
    boolean isUnivalent[MAXN];
    int nrUnivalent = 0, nrNonUnivalent=0;
    int oldindtonew[MAXN], newindtoold[MAXN];
    for (int i=0;i<n;i++) {
        int ecount =0;
        for (int j=0;j<n;j++) 
            if (i!=j && TISELEMENT(g[i],j)) 
                ecount++;
        if (ecount == 1)
        {
            isUnivalent[i] = TRUE;
            nrUnivalent++;
        }
        else 
        {
            isUnivalent[i] = FALSE;
            oldindtonew[i] = nrNonUnivalent;
            newindtoold[nrNonUnivalent] = i;
            nrNonUnivalent++;
        }
    }
    
    // attach an extra edge at the end of each hair to change symmetry
    setword gg[MAXN]={0};
	int i,j, N=n+nrUnivalent;
	
	for (i=0;i<N;i++) 
	  for (j=0;j<i;j++)
	    if (TISELEMENT(g[i], j))
	       TADDELEMENT(gg[i],j);
	
	/*{ 
	  for (j=0;j<n;j++) if (i!=j && !isUnivalent[j] && TISELEMENT(g[j], g[newindtoold[i]]))
	  {
	      TADDELEMENT(gg[i],oldindtonew[j]); 
	  }
 	}*/
	
	
    graph gcan[MAXN];
    return makecanon2(gg,gcan,N,NULL, evenedges);
}




/*int permsign(int* perm, int n)
{
	int i,j, res=1;
	boolean vistd[MAXN]={0};
	//for (i=0;i<n;i++) printf(" %d",perm[i]);
	//return 1;
	//for (i=0;i<n;i++) vistd[i]=FALSE;
	for (i=0;i<n;i++) if (!vistd[i]) {
	  j=i;
	  res=-res;
	  do {
	    vistd[j] = TRUE;
	    j = perm[j];
	  //  printf(" %d",j);
	    res=-res;
	  } while (j!=i);
	}
	return res;
}*/
int contract(graph* g, int n, int i, int j, graph* gout, boolean evenedges)
// takes graph g, contracts edge (ij) and puts the result in gout. The return value is the sign, or 0 in case of a zero graph (double edge etc)
{
	setword temp[50]={0};
	int inds[30][30]={0}, perm1[MAXN]={0}, perm2[MAXN]={0};
	int k,l,ii,newk, newl, sign, nre=0;
	int sign_from_vertex_order, sign_from_edge_orientations=1, sign_from_canonisation;
	if (!TISELEMENT(g[i],j)) return 0;

	// create the new graph (not yet canonised, stored in temp)
	// vertex j is deleted from the graph. The indices of vertices <j remain the same, 
	// those of vertices > j are reduced by one (and j "becomes" i)
	for (k=0;k<n-1;k++) for (l=k+1;l<n;l++) 
	if (TISELEMENT(g[k],l)) 
	{
	  newk = (k<j?k: (k==j?i:k-1));
	  newl = (l<j?l: (l==j?i:l-1));
	  if (TISELEMENT(temp[newk],newl)){ return 0;}
	  if (newk != newl) {
	    TADDELEMENT(temp[newk], newl);
	    TADDELEMENT(temp[newl], newk);
	    if (newk>newl)
	      sign_from_edge_orientations = - sign_from_edge_orientations;  // only used for even edges
	  }
	  inds[newk][newl]=inds[newl][newk]=nre++;      // only used for odd edges
	}

	// canonise
	if (!makecanon2(temp, gout, n-1, perm2, evenedges)) return 0;
	// canonised graph is now in gout. Determine the sign
	if (evenedges)
	{
	    sign_from_vertex_order = (j%2 == 0? 1 : -1);
	    if (i>j)
	    {
	      sign_from_edge_orientations = -sign_from_edge_orientations;
	    }
	    sign_from_canonisation = getgpermsign_even(temp, n-1, perm2); // gout <-> temp??
	    return sign_from_edge_orientations * sign_from_vertex_order * sign_from_canonisation;
	}
	else
	{
	    perm1[0] = inds[i][i];
	    ii =1;
	    for (k=0;k<n-2;k++) for (l=k+1;l<n-1;l++) if (TISELEMENT(gout[k],l)) {
	      perm1[ii++] = inds[perm2[k]][perm2[l]];
	    }
        return getpermsign(perm1, ii);  // TODO:use getgpermsign
    } 
}


/**************************************************************************/
void genwheelgraph(graph* gout, int n)
// generates wheel with 2n+1 spokes
{
	graph gg[MAXN];
	int i,j, N=2*n+2;
	for (i=0;i<2*n+2;i++) gg[i]=0;
	for (i=0;i<2*n+1;i++) { 
	  TADDELEMENT(gg[i],2*n+1); // to center
	  TADDELEMENT(gg[2*n+1],i);
	  j=(i+1) % (2*n+1);

	  TADDELEMENT(gg[i], j); // to next
	  TADDELEMENT(gg[j], i); // to next
 	  //TADDELEMENT(gg[i],(i-1)% (2*n+1) ); // to previous
 	}
	// canonise
	makecanon2(gg,gout,n*2+2,NULL, FALSE);
}
/**************************************************************************/


int readgraphmap(char* cFile, graphmap& m)
{
	ifstream f (cFile, ios::in);
	string s;
	int count=0;
	if (!f.is_open()) { cerr <<  cFile << ": Cannot open file\n"; return 0; }
	
	while (!f.eof()) {
	  getline(f, s);
	  if (s.length() > 1) {
	    m[s]=(++count);
	    //cout << s << endl;
	  }
	}
	f.close();
	return count;
}

int readgraphvect(char* cFile, graphvect& v)
/* File format: x (coeff) (graph in g6) */
{
	ifstream f (cFile, ios::in);
	string s, gra, stemp;
	double coeff;
	int count=0;
	if (!f.is_open()) { cerr <<  cFile << ": Cannot open file\n"; return 0; }
	
	while (!f.eof()) {
	  getline(f, s);
	  if (s.length() > 1) {
	    if (s[0] == 'x') {
	      istringstream iss(s);
	      iss >> stemp >> coeff >> gra;
	      v[gra] += coeff;
	      count++;
	    }
	  }
	}
	f.close();
	return count;
}

string gratostring(graph *g, int n)
{
    string res;
    res = ntog6(g,1,n);
    // purge the newline
    res.erase(res.length()-1,1);
    return res;
}

int stringtogra(string s, graph *gout)
{
  char c[100];
  strcpy(c, s.c_str());
  stringtograph(c, gout, 1);
  return graphsize(c);
}
void printset(FILE* f, setword s, int n)
{
	for (int i=0;i<n;i++)
		if (TISELEMENT(s, i)) 
			fprintf(f,"1");
		else	fprintf(f,"0");
	fprintf(f,"\n");
}
void edgereconnect(graph* g, int n, int n1, setword torec, double pref, graphvect& res)
{
	int i, v;
	setword temp[MAXN]={0};
	int perm[MAXN];
	string gra;

	// find next vertex to connect
	for (v=0;v<n1;v++) if (TISELEMENT(torec, v)) break;
	if (v==n1) {
//fprintf(stderr,"-%d-", nredges(g,n) );
//		fprintf(stderr,"in");
	  // no vertex left => output the (canonised) graph and return
	  if (!makecanon2(g, temp, n, perm, FALSE)) { return; } // TODO: even edge case
//	  fprintf(stderr,"in2");
//printperm(perm, n);
          invperm(perm,n);
	  int sign = getgpermsign(temp,n,perm);
//	  fprintf(stderr,"in3");
	  gra = gratostring(temp, n);
//cerr << gra <<  " " << sign*pref<<endl;
	  res[gra] += sign*pref;
	  //fprintf(f, "x%d ", sign*pref);
	  //writeg6(f, temp, 1, n);
	  //fprintf(stderr, ".");
//	  fprintf(stderr,"out");
	  return;
	}

	// v is the vertex this instance reconnects
	TDELELEMENT(torec, v);
//	printset(stderr, torec, n1);
	for (i=n1;i<n;i++) {
	  // add edge from v to i
	  TADDELEMENT(g[v],i);  TADDELEMENT(g[i],v); 
	  // recursively add the next edge
	  edgereconnect(g, n, n1, torec, pref, res);
	  // undo changes before adding a different edge
	  TDELELEMENT(g[v],i);  TDELELEMENT(g[i],v);
	}	
}

void insertatvert(graph *g1, int n1, graph *g2, int n2, int v, double pref, graphvect& res)
/* insert g2 into g1 at vertex v */
{
fprintf(stderr,".");
	// permute the vertices so that v becomes the last
	int perm[MAXN], i;
	setword temp[MAXN]={0};
	for (i=0;i<n1;i++) perm[i] = (i<v?i: (i==v?n1-1:i-1) );
	int sign = graphperm(g1, n1, perm, temp);
//fprintf(stderr,"%d", nredges(temp,n1) );
	// edges to be reconnected
	setword torec = temp[n1-1];
	for (i=0;i<n1-1;i++) TDELELEMENT(temp[i], n1-1);
	// merge graphs
	for (i=0;i<n2;i++) {
	  temp[n1-1+i] = (g2[i] >> (n1-1));
	}
//fprintf(stderr,"+%d", nredges(temp,n1+n2-1) );
//printset(stderr, torec, n1);
	// reconnect edges
	edgereconnect(temp, n1+n2-1, n1-1, torec, sign*pref, res);
}

void graphprebracket(const graphvect G1, const graphvect G2, graphvect& res)
{
    graph g1[MAXN], g2[MAXN];
    graphvect::const_iterator g1it, g2it;
    int n1, n2;
    for (g1it = G1.begin(); g1it != G1.end(); ++g1it) {
        n1 = stringtogra(g1it->first, g1);
        for (g2it = G2.begin(); g2it != G2.end(); ++g2it) {
            n2 = stringtogra(g2it->first, g2);
            for (int i=0;i<n1;i++) {
                insertatvert(g1,n1,g2,n2,i,g1it->second * g2it->second,res);
	    }
        }
    }
}

void scalarmult(graphvect& v, double lambda)
{
    graphvect::iterator vit;
    for (vit = v.begin(); vit != v.end(); vit++)
        vit->second *= lambda;
}
//TODO: fehlende initialisierung von doubles!!!
//TODO: extrasign for odd graphs
graphvect bracket(const graphvect G1, const graphvect G2)
{
    graphvect res;
    int sign = -1, deg1=0, deg2=0;

    // determine sign (assumes homogeneous vectors given
    graphvect::const_iterator it;
    graph temp[MAXN];
    it = G1.begin();
    if (it != G1.end()) {
	int n = stringtogra(it->first, temp);
	deg1 = 2*n - nredges(temp, n) -2;
    }
    it = G2.begin();
    if (it != G2.end()) {
	    int n = stringtogra(it->first, temp);
	    deg2 = 2*n - nredges(temp, n) -2;
    }
    if ( (deg1*deg2)%2 != 0) sign=1;

    graphprebracket(G2, G1, res);
    scalarmult(res, sign);
    graphprebracket(G1, G2, res);

    return res;
}

void printgraphvect(FILE* f, const graphvect v)
{
	graphvect::const_iterator vit;
	for (vit = v.begin(); vit != v.end(); vit++)
	    fprintf(f, "x %.20f %s\n", vit->second, vit->first.c_str());
}
void writegraphvect(char* cFile, const graphvect v)
{
	FILE* f = fopen(cFile, "w");
	printgraphvect(f, v);
	fclose(f);
}

#define ABS(x) (x<0?-x:x)
void cleanv(graphvect& v)
{
	graphvect::iterator vit;
	for (vit = v.begin(); vit != v.end(); vit++) {
		if (ABS(vit->second)<TOLER) 
			v.erase(vit);
	}
	
}

void addto(graphvect& v1, const graphvect v2)
{
	graphvect::const_iterator vit;
	for (vit = v2.begin(); vit != v2.end(); vit++) {
		v1[vit->first] += vit->second;
	}
}	

void dualbasis2coinvariants(graphvect& v)
{
	graphvect::iterator vit;
	graph g[MAXN];
	int n;
	for (vit = v.begin(); vit != v.end(); vit++) {
		n = stringtogra(vit->first, g);
		vit->second /= ordersymmgroup(g,n);
	}
}	
void coinvariants2dualbasis(graphvect& v)
{
	graphvect::iterator vit;
	graph g[MAXN];
	int n;
	for (vit = v.begin(); vit != v.end(); vit++) {
		n = stringtogra(vit->first, g);
		vit->second *= ordersymmgroup(g,n);
	}
}
void countsymm(permutation *p, int n)
/* Called by allgroup. Just counts the number of group elements. */
{
	ggrcount++;
}
int ordersymmgroup(graph* g, int n)
{
	static DEFAULTOPTIONS_GRAPH(options);
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
	setword workspace[150];
	statsblk stats;
	//	set *gv; /////
		grouprec *group;
/* The following cause nauty to call two procedures which
		store the group information as nauty runs. */
	options.userautomproc = groupautomproc;
	options.userlevelproc = grouplevelproc;

	nauty(g,lab,ptn,NULL,orbits,&options,&stats,
				      workspace,150,1,n,NULL);
/* Get a pointer to the structure in which the group information
				has been stored. If you use TRUE as an argument, the
				structure will be "cut loose" so that it won't be used
				again the next time nauty() is called. Otherwise, as
				here, the same structure is used repeatedly. */
	group = groupptr(FALSE);
/* Expand the group structure to include a full set of coset
				representatives at every level. This step is necessary
				if allgroup() is to be called. */
	makecosetreps(group);
/* Call the procedure writeautom() for every element of the group.
				The first call is always for the identity. */
	ggrcount = 0;
	allgroup(group,countsymm);
	if (ggrcount<1) fprintf(stderr, "Error in ordersymm: ggrcount=%d\n", ggrcount);
	return ggrcount;
}
