
/* 
Wrapper library for nauty, to be called by julia.
*/

//#define MAXN 1000    /* Define this before including nauty.h */
#include "nauty.h"   /* which includes <stdio.h> and other system files */
#include "gtools.h" 

typedef  unsigned long long juliagraph ;

typedef int bool;

int* gAutomInfo;
int automCount;


bool isjuliabitset(juliagraph* g, int i,int j,int n)
{
	int bitpos, wordpos, bitpos2,s;
	juliagraph w;
	bitpos2 = j * n +i;
	s = sizeof(juliagraph) * 8;
	wordpos = bitpos2 / s;
	bitpos = bitpos2 % s;
	//fprintf(stderr, "    _ %d %d %d %d\n", wordpos, bitpos, bitpos2, s);
	w = g[wordpos];
	//fprintf(stderr, "Hallo%llx", w);

	return (w & (1LL << bitpos)) != 0LL;
	//return g[wordpos] & (1 << (sizeof(juliagraph)-bitpos-1));
}

int
test_plus(int a, int b)
{
	graph g[MAXN*MAXM];
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
	static DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;
	int n,m,v;

	return a + b+3;

	/* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
	Here we change those options that we want to be different from the
	defaults.  writeautoms=TRUE causes automorphisms to be written.     */
//	options.writeautoms = TRUE;
	
	/* The nauty parameter m is a value such that an array of
	m setwords is sufficient to hold n bits.  The type setword
	is defined in nauty.h.  The number of bits in a setword is
	WORDSIZE, which is 16, 32 or 64.  Here we calculate
	m = ceiling(n/WORDSIZE).                                  */
	//m = SETWORDSNEEDED(n);
	//42
	/* The following optional call verifies that we are linking
	to compatible versions of the nauty routines.            */
	//nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
	/* Now we create the cycle.  First we zero the graph, than for
	each v, we add the edge (v,v+1), where values are mod n. */
	//EMPTYGRAPH(g,m,n);
	//for (v = 0; v < n; ++v)  ADDONEEDGE(g,v,(v+1)%n,m);
	//printf("Generators for Aut(C[%d]):\n",n);
	/* Since we are not requiring a canonical labelling, the last
	parameter to densenauty() is not required and can be NULL. */
	//densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
	/* The size of the group is returned in stats.grpsize1 and
	stats.grpsize2. */
	//printf("Automorphism group size = ");
	//writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
	//printf("\n");
	//}
	// xit(0); 
}

int
copyGraph( juliagraph* data, int n, graph* gout)
{
	int m,u,v;
	m = SETWORDSNEEDED(n);
	//nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
	EMPTYGRAPH(gout,m,n);
	//fprintf(stderr, "Hallo%llx", data[0]);
	//for (v = 0; v < n; ++v) 
	//	ADDONEEDGE(gout,v,(v+1)%n,m);
	for (u=1;u<n;u++)
	{
		for (v=0;v<u;v++)
		{
			if (isjuliabitset(data, v, u, n))
			{
				ADDONEEDGE(gout, u,v,m);
				//fprintf(stderr,"x");
				//fprintf(stderr, "%d %d \n",u,v);
			}
		}
	}
	return m;
}

int
copyGraphdi( juliagraph* data, int n, graph* gout)
{
	int m,u,v;
	m = SETWORDSNEEDED(n);
	//nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
	EMPTYGRAPH(gout,m,n);
	//fprintf(stderr, "Hallo%llx", data[0]);
	//for (v = 0; v < n; ++v) 
	//	ADDONEEDGE(gout,v,(v+1)%n,m);
	for (u=0;u<n;u++)
	{
		for (v=0;v<n;v++)
		{
			if (isjuliabitset(data, v, u, n))
			{
				ADDONEARC(gout, u,v,m);
				//fprintf(stderr,"x");
				//fprintf(stderr, "%d %d \n",u,v);
			}
		}
	}
	return m;
}

void
revCopyGraph( graph* g, int m, int n, juliagraph* dataout)
{
	int i,j,cnt;
	juliagraph x, maxx;
	set* gj;

	//fprintf(stderr, "Hallo%llx", data[0]);
	//for (v = 0; v < n; ++v) 
	//	ADDONEEDGE(gout,v,(v+1)%n,m);
	
	// clear out buffer
	cnt = n*n / (8* sizeof(juliagraph));
	if (cnt* 8 * sizeof(juliagraph) < n*n)
		cnt++;

	for (j=0;j<cnt;j++)
		dataout[j] = 0LL;

	cnt = 0;
	maxx = 1LL << (sizeof(juliagraph)*8-1);
	x = 1;
	for (j=0;j<n;j++)
	{
		gj = GRAPHROW(g,j,m);
		for (i = 0; i < n; ++i)
        {
            if (ISELEMENT(gj,i)) 
            {
            	dataout[cnt] |= x;
            	//fprintf(stderr,"*%llx",x);
            }
            if (x== maxx )
            {
            	x = 1;
            	cnt++;
            }
            else
            {
            	x <<=1;
            }
		}
	}
}

void
revCopyGraphdi( graph* g, int m, int n, juliagraph* dataout)
{
	int i,j,cnt;
	juliagraph x, maxx;
	set* gj;

	//fprintf(stderr, "Hallo%llx", data[0]);
	//for (v = 0; v < n; ++v) 
	//	ADDONEEDGE(gout,v,(v+1)%n,m);
	
	// clear out buffer
	cnt = n*n / (8* sizeof(juliagraph));
	if (cnt* 8 * sizeof(juliagraph) < n*n)
		cnt++;

	for (j=0;j<cnt;j++)
		dataout[j] = 0LL;

	cnt = 0;
	maxx = 1LL << (sizeof(juliagraph)*8-1);
	x = 1;
	for (j=0;j<n;j++)
	{
		gj = GRAPHROW(g,j,m);
		for (i = 0; i < n; ++i)
        {
            if (ISELEMENT(gj,i)) 
            {
            	dataout[cnt] |= x;
            	//fprintf(stderr,"*%llx",x);
            }
            if (x== maxx )
            {
            	x = 1;
            	cnt++;
            }
            else
            {
            	x <<=1;
            }
		}
	}
}

char*
generateg6( juliagraph* data, int n)
{
	// Leaks memory, so dont use... this is just for testing
	graph g[MAXN*MAXM];
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
	static DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;
	int m;
	char* c;

	m = copyGraph(data, n, g);
	c = ntog6(g, m, n);

	//fprintf(stderr, "%s",c);
	return c;

	/* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
	Here we change those options that we want to be different from the
	defaults.  writeautoms=TRUE causes automorphisms to be written.     */
//	options.writeautoms = TRUE;
	
	/* The nauty parameter m is a value such that an array of
	m setwords is sufficient to hold n bits.  The type setword
	is defined in nauty.h.  The number of bits in a setword is
	WORDSIZE, which is 16, 32 or 64.  Here we calculate
	m = ceiling(n/WORDSIZE).                                  */
	//m = SETWORDSNEEDED(n);
	//42
	/* The following optional call verifies that we are linking
	to compatible versions of the nauty routines.            */
	//nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
	/* Now we create the cycle.  First we zero the graph, than for
	each v, we add the edge (v,v+1), where values are mod n. */
	//EMPTYGRAPH(g,m,n);
	//for (v = 0; v < n; ++v)  ADDONEEDGE(g,v,(v+1)%n,m);
	//printf("Generators for Aut(C[%d]):\n",n);
	/* Since we are not requiring a canonical labelling, the last
	parameter to densenauty() is not required and can be NULL. */
	//densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
	/* The size of the group is returned in stats.grpsize1 and
	stats.grpsize2. */
	//printf("Automorphism group size = ");
	//writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
	//printf("\n");
	//}
	// xit(0); 
}


int testing(juliagraph *data, int n)
{
	int m;
	graph g[MAXN*MAXM];
	m = copyGraph(data, n, g);
	revCopyGraph(g, m, n, data);
	return 0;
}
/* Just reads the graph to internal format, and then writes it again inplace.
   If all goes well, the result is just unchanged.
   This is used solely as a test of the conversion routines from/to Julia format.  */
int testingdi(juliagraph *data, int n)
{
	int m;
	graph g[MAXN*MAXM];
	m = copyGraphdi(data, n, g);
	revCopyGraphdi(g, m, n, data);
	return 0;
}


void
userautomproct(int count, permutation *p, int *orbits,
              int numorbits, int stabvertex, int n)
/* called by nauty;  operates on glIsEven */
{
    int i;
    for (i=0;i<n;++i) {
        gAutomInfo[n*automCount+i] = p[i];
    }
    automCount++;
}


/*
	Canonize the input (in place) and output the permutation, and optionally generators of the automorphism group in permout
	We assume that permout is large enough.
	Indices of vertices are 0 based.
	The return value is the number of generators.
	All permutations (generators plus the canonizing one) are written consecutively.
	I.e., if r is returned, (r+1)n ints will be taken in permout
*/
int canonlabel(juliagraph *data, int n, bool* colors, int* permout, bool get_isoms, bool get_canonperm)
{
	int m, i;
	graph g[MAXN*MAXM];
	graph gcanon[MAXN*MAXM];
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
	DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;

	m = copyGraph(data, n, g);

	//fprintf(stderr, "params: %d %d\n" , get_isoms, get_canonperm);
	   	options.getcanon = TRUE;
	   	if (get_isoms)
	   	{
	   		//fprintf(stderr, "automorphisms collected");
        	options.userautomproc = userautomproct;
        }
        if (colors)
        {
        	for (i=0;i<n;i++)
        	{
        		lab[i]=i;
        		ptn[i]=colors[i];
        	}
        	options.defaultptn = FALSE;
        }
        //else
        //	fprintf(stderr, "automorphisms not collected");
        gAutomInfo = permout;
        automCount=0;

        //EMPTYSET(active,m);
        densenauty(g,lab,ptn,orbits,&options,&stats,m,n,gcanon);
        // display permutation
        //fprintf(stderr, ";");

        if (get_canonperm)
        {
        	for (i=0;i<n;++i) {
            	permout[n*automCount+i] = lab[i];
        	}
        }
        //sprintf(buf, "num%d",stats.numorbits);
        //strcat(autominfo, buf);
        //gt_numorbits = stats.numorbits;

        revCopyGraph(gcanon, m,n,data);
        return automCount;
}


int canonlabeldi(juliagraph *data, int n, bool* colors, int* permout, bool get_isoms, bool get_canonperm)
{
	int m, i;
	graph g[MAXN*MAXM];
	graph gcanon[MAXN*MAXM];
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
	//DEFAULTOPTIONS_DIGRAPH(options);
	DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;

	m = copyGraphdi(data, n, g);


	   	options.getcanon = TRUE;
	   	if (get_isoms)
	   	{
	   		//fprintf(stderr, "automorphisms collected");
        	options.userautomproc = userautomproct;
        }
        if (colors)
        {
        	for (i=0;i<n;i++)
        	{
        		lab[i]=i;
        		ptn[i]=colors[i];
        	}
        	options.defaultptn = FALSE;
        }
        //else
        //	fprintf(stderr, "automorphisms not collected");
        gAutomInfo = permout;
        automCount=0;
        options.digraph = TRUE;

        //EMPTYSET(active,m);
        densenauty(g,lab,ptn,orbits,&options,&stats,m,n,gcanon);
        // display permutation
        //fprintf(stderr, ";");

        if (get_canonperm)
        {
        	for (i=0;i<n;++i) {
            	permout[n*automCount+i] = lab[i];
        	}
        }
        //sprintf(buf, "num%d",stats.numorbits);
        //strcat(autominfo, buf);
        //gt_numorbits = stats.numorbits;

        revCopyGraphdi(gcanon, m,n,data);
        return automCount;
}



