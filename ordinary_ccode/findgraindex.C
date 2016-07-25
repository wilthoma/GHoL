#include <iostream>
#include <stdio.h>
#include <iostream>
#include "tgra.h"

using namespace std;

int main(int argc, char* argv[])
{
	// Simply print out the wheel graph
	if (argc != 2) {
		cout << "Usage: findgraindex file.g6\n";
		return -1;
	}
	graphmap gm;
	readgraphmap(argv[1],gm);
	
	while (TRUE) {
		string grain;
		cin >> grain;
		if (gm.count(grain)>0)
			cout << gm[grain] << endl;
		else cout << "Not found\n";
	}
	return 0;
}