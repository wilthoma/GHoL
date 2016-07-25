#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>

#include "tgra.h"

using namespace std;

int main(int argc, char* argv[])
{
	// Simply print out the wheel graph
	if (argc != 3) {
		cout << "Usage: \n mattohg6 matfile.txt gra.g6 > outfile.hg6 \n Converts indices in matrix file to .hg6\n";
		return -1;
	}
	
	// read graph file
	vector<string> graDB;
	string line;
	graDB.clear();
	ifstream infile (argv[2], ios_base::in);
	while (getline(infile, line, '\n'))
	{
		graDB.push_back (line);
	}
	infile.close();

	ifstream matfile (argv[1], ios_base::in);
	int count =0;
	while (getline(matfile, line, '\n'))
	{
		cout << "x " << line << " " << graDB[count++] << endl;
		
	}
	matfile.close();
		
	
	return 0;
}