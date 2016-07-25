#include <iostream>
#include <stdio.h>
#include <iostream>
#include "tgra.h"

using namespace std;

int main(int argc, char* argv[])
{
	// Simply print out the wheel graph
	if (argc != 2) {
		cout << "Usage: genwheelgra n \n Prints the wheel with 2n+1 spokes\n";
		return -1;
	}
	int n = atoi(argv[1]);
	graph g[MAXN];
	genwheelgraph(g, n);
	cout << gratostring(g, 2*n+2) << endl;
	return 0;
}