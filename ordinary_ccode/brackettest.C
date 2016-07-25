#include <iostream>
#include <stdio.h>
#include <iostream>
#include "tgra.h"

using namespace std;

int main(int argc, char* argv[])
{
	// read in three files and check the Jacobi id. on them
	graphvect v1,v2,v3;
	if (argc<3) return -1;
	readgraphvect(argv[1], v1);
	readgraphvect(argv[2], v2);
	readgraphvect(argv[3], v3);
	
	printf("v1:\n");
	printgraphvect(stdout, v1);
	printf("\nv2:\n");
	printgraphvect(stdout, v2);
	printf("\nv3:\n");
	printgraphvect(stdout, v3);


	graphvect res1, res2, res3, restmp;
	
	restmp = bracket(v2,v3);
	cleanv(restmp);
/*	printgraphvect(stdout, restmp);
	FILE* f = fopen("testen.g6","w");
	printgraphvect(f, restmp);
	fclose(f);
	*/	
	//printgraphvect(stderr, restmp);	
	
	cerr << "X";
	res1 = bracket(v1,restmp);
	cerr << "XX";
	
	restmp=bracket(v3,v1);
	//printgraphvect(stderr, restmp);
	cleanv(restmp);
	cerr << "Y";
	res2 = bracket(v2,restmp);
	cerr << "YY";
	//res2 = bracket(v3,v1);
	
	restmp=bracket(v1,v2);
	cleanv(restmp);
	cerr << "Z";
	res3 = bracket(v3,restmp);
	cerr << "ZZ";
	cleanv(res1); cleanv(res2); cleanv(res3);
	//printgraphvect(stderr, res3);
  //  scalarmult(res3,-1);	
	addto(res1, res2);
	addto(res1, res3);
	cleanv(res1);
	
	printf("\nres1:\n");
	printgraphvect(stderr, res1);
	printf("\nres2:\n");
	//printgraphvect(stdout, res2);
	printf("\nres3:\n");
	//printgraphvect(stdout, res3);

	return 0;
}