// Header for thomas graph functions
#ifndef  _TGRA_H_    /* only process this file once */
#define  _TGRA_H_


#define MAX(i,j) ((i)>(j)?(i):(j))
#define MIN(i,j) ((i)>(j)?(j):(i))
#define TISELEMENT(setadd,pos)  (((setadd) & BITT[pos]) != 0)
#define TADDELEMENT(setadd,pos)  ((setadd) |= BITT[pos])
#define TDELELEMENT(setadd,pos)  ((setadd) &= ~BITT[pos])


//#ifdef __cplusplus
//extern "C" {
//#endif

#include <map>
#include <string>
#include "nauty/gtools.h"

typedef std::map<std::string, int> graphmap;             // A map mapping each graph to an index (... in a matrix)
typedef std::map<std::string, double> graphvect;         // represents a linear combination of graphs

int readgraphmap(char* cFile, graphmap& m);     // returns count
int readgraphvect(char* cFile, graphvect& v);
void writegraphvect(char* cFile, const graphvect v);
void printgraphvect(FILE* f, const graphvect v);
void addto(graphvect& v1, const graphvect v2);
#define TOLER 1e-5
void cleanv(graphvect& v);	// REMOVE all entries ob absolute value < TOLER
void scalarmult(graphvect& v, double lambda);
graphvect bracket(const graphvect G1, const graphvect G2);

void dualbasis2coinvariants(graphvect& v);
void coinvariants2dualbasis(graphvect& v);

int ordersymmgroup(graph* g, int n);

void invperm(permutation*, int);
int getpermsign(permutation* perm, int n);
int getgpermsign(graph*, int, permutation*);
int graphperm(graph*, int, permutation*, graph*);
// checks whether graph is valid (has no odd symmetries)
boolean isNonzeroGraph(graph *g, int n, boolean evenedges);

boolean isNonzeroGraphOppositeHairParity(graph *g, int n, boolean evenedges);

// takes graph g, contracts edge (ij) and puts the result in gout. The return value is the sign, or 0 in case of invalid graph (double edge etc)
int contract(graph*, int, int, int, graph*, boolean evenedges);

void genwheelgraph(graph* gout, int n);

std::string gratostring(graph *g, int n);
int stringtogra(std::string s, graph *gout);

boolean makecanon2(graph *g, graph *gcan, int n, int* perm, boolean evenedges);

//#ifdef __cplusplus
//}
//#endif


#endif
