#define USAGE \
"Usage:\n\
   makemat list.g6 differential.hg6\n\
   makemat list.g6 differential.hg6 > sparsematrix.txt\n\
\n\
Here:\n\
   list.g6 is a file containing a list of graphs in degree n-1 (generated with gradiff -l ...)\n\
   differential.hg6 is a file containing the differential from degree n to n-1 (generated with gradiff ...)\n\
\n\
The output will be in the standard sparse matrix format 'rowindex colindex entry'\n"

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <iomanip>

using namespace std;

typedef map<string, int> matmap;

int readDB(char* cFile, matmap& m)
{
	ifstream f (cFile, ios::in);
	string s;
	int count=0;
	if (!f.is_open()) { fprintf(stderr, "%s: Cannot open file\n", cFile); return 0; }
	
	while (!f.eof()) {
	  getline(f, s);
	  if (s.length() > 1) {
	    if (s[0]!= 'x')
	      m[s]=(++count);
	    //cout << s << endl;
	  }
	}
	f.close();
	return count;
}

int
main(int argc, char *argv[])
{
	if (argc!=3) {
	  cerr << USAGE;
	  return 1;
	}
	matmap mimage;
	int imcount = readDB(argv[1], mimage);
	fprintf(stderr, "%d graphs read from %s\n", imcount, argv[1]);
	
	// read the input, format: either (graph in g6) or x(int i) (graph in g6)
	int count =0;
	double coeff;
	string s, gra, stemp;
	istringstream iss(s);
	ifstream f (argv[2]);
	if (!f.is_open()) { fprintf(stderr, "%s: Cannot open file\n", argv[2]); return 0; } 
	while (!f.eof()) {
	  getline(f, s);
	  if (s.length() > 1) {
	    if (s[0] == 'x') {
	      istringstream iss(s);
	      //iss.seekp(ios:beg);
	      iss >> stemp >> coeff >> gra;
	      cout << count << " " << mimage[gra] << " " << setprecision(20) << coeff << endl;
	    }
	    else count++; 
	  }
	}
	f.close();
	
	// print out the size
	cout << count << " " << imcount << " " << 0 << endl;
	
	cerr << count << " graphs read from " << argv[2] << endl; 
	return 0;
}
