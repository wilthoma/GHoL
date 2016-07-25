#include <iostream>
#include <stdio.h>
#include <iostream>
#include "tgra.h"
#include <iomanip>

using namespace std;

int main(int argc, char* argv[])
{
	setprecision(20);

	if (argc !=4) {
		cout << "Usage:\n compbracket v1.hg6 v2.hg6 outfile.hg6\n";
		return -1;
	}
	
	graphvect v1,v2,v3, v4;
	readgraphvect(argv[1], v1);
	readgraphvect(argv[2], v2);
	cleanv(v1);
	cleanv(v2);
	dualbasis2coinvariants(v1);
	dualbasis2coinvariants(v2);
	
	printgraphvect(stdout, v1);
	printgraphvect(stdout, v2);
	
	v4 = bracket(v1,v2);
	//v4 = bracket(v1,v3);
	cleanv(v4);
	coinvariants2dualbasis(v4);
	//printgraphvect(stdout, v3);
	cleanv(v4);
	writegraphvect(argv[3], v4);

	return 0;
}