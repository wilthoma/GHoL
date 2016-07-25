#include <iostream>
#include <stdio.h>
#include <iostream>
#include "tgra.h"
#include <iomanip>

using namespace std;
int markCC(graph* g, int n, int startv, bool* vistd)
{
	vistd[startv]=true;
	int nrmarked = 1;
	for (int j=0;j<n;j++) {
		if (startv !=j && TISELEMENT(g[startv],j)) {
			if (!vistd[j])
				nrmarked += markCC(g,n,j,vistd);
		}
	}
	return nrmarked;
}
int nrICC(graph* g, int n, int v1, int v2, int* ICCSizes=NULL)
{
	bool vistd[MAXN]={0};
	//int ICCsizes[MAXN]={0};
	vistd[v1]=true;
	vistd[v2]=true;
	int cccount = 0;
	for (int i=0;i<n;i++) {
		if (!vistd[i]) {
			int temp=markCC(g,n,i,vistd);
			if (ICCSizes)
				ICCSizes[cccount]=temp;
			cccount++;
		}
	}
	
	return cccount;
}

int getarmlength(graph* g, int n, int lowv, int hiv, int src, int tgt, int& finaltgt)
{
	if (!TISELEMENT(g[tgt], hiv))
	{
		finaltgt = tgt;
		return 0;
	}
	
	// find next target
	for (int i=0;i<n;i++) {
		if (i!=hiv && i!= src && i!=tgt && TISELEMENT(g[tgt],i))
		{
			return (getarmlength(g,n,lowv, hiv,tgt,i, finaltgt)+1);
		}
	}
	
	cerr << "Error in getarmlength()\n";
	return -10000;
}

bool getgrlengths(graph* g, int n, int* hv, int* hvc, int* lengthsout)
/* output l1 l2 r1 r2 middle */
{
	int lowval=6, lowv, hiv, lengths[MAXN], tgts[MAXN], cnt=0;
	if (hvc[0]==lowval)
		{ lowv=hv[0]; hiv=hv[1]; }
	else if (hvc[1]==lowval)
		{ lowv=hv[1]; hiv=hv[0]; }
	else return false;
	
	// loop through all verts adjacent to lowv, find the root index
	for (int i=0;i<n;i++) {
		if (i!=hiv && i!= lowv && TISELEMENT(g[lowv],i)) {
			lengths[cnt] = getarmlength(g,n,lowv,hiv,lowv,i,tgts[cnt]);
			cnt++;
		}
	}
	
	if (cnt!=5) { cerr << "Error: can handle only 5 arm-type\n"; return false; }
	
	// find the double final tgts
	int doubletgts[MAXN], dtcnt=0;
	for (int i=0;i<cnt-1;i++)
		for (int j=i+1;j<cnt;j++)
			if (tgts[i]==tgts[j]) {
				doubletgts[dtcnt]=tgts[i];
				// dangerous
				lengthsout[2*dtcnt] = lengths[i];
				lengthsout[2*dtcnt+1] = lengths[j];
				dtcnt++;
			}
	
	if (dtcnt!= 2) {
		 cerr << "Error: can handle only 5 arm-type (2dt) "<<dtcnt<<endl;
		 cerr << gratostring(g, n) << endl;
		 for (int i=0;i<5;i++) cerr << "--"<<tgts[i]<<" "<<lengths[i];
		 cerr << endl;
		 return false; }
	
	// find the single target
	for (int i=0;i<cnt;i++) {
		if (tgts[i]!=doubletgts[0] && tgts[i]!=doubletgts[1])
		{
			lengthsout[4] = lengths[i];
			break;
		} 
	}
	
	return true;
}
bool checkgra(graph* g, int n)
{
	// Find the two high valence vertices
	int hv[10], hcount =0, hvc[10];
	int ICCsizes[MAXN]={0}, armlengths[MAXN];
	
	for (int i=0;i<n;i++) {
		int ecount =0;
		for (int j=0;j<n;j++) 
			if (i!=j && TISELEMENT(g[i],j)) 
				ecount++;
		if (ecount >= 4) {
			hvc[hcount]=ecount;
			hv[hcount++] = i;
		}
	}
	
	if (hcount < 2) return true;
	
	if (!TISELEMENT(g[hv[0]],hv[1])) return false;
	
	int nrconnc = nrICC(g,n,hv[0],hv[1], ICCsizes);
	
	if (nrconnc>1) return false; //needs three conn. comp part
	//if (nrconnc==3 && (hvc[1]==n-4 || hvc[0]==n-4) ) return false; //no
//	if (nrconnc==3 && (hvc[1]==n-3 || hvc[0]==n-3) ) return false;
//	if (nrconnc==3 && (hvc[1]==n-2 || hvc[0]==n-2) ) return false;
	//if (nrconnc==1 && (hvc[1]==n-4 || hvc[0]==n-4) ) return false; // no
	
	if (nrconnc==1 
		   //&& (hvc[1]!=4 && hvc[0]!=4) 
		   //    && (hvc[1]!=5 && hvc[0]!=5)
		   && (hvc[1]!=6 && hvc[0]!=6)
	   ) return false; //no
	
/**/	if (nrconnc==2 
		   //&& (hvc[1]!=4 && hvc[0]!=4) 
		   //    && (hvc[1]!=5 && hvc[0]!=5)
		   && (hvc[1]!=6 && hvc[0]!=6)
	   ) return false; //no
	if (nrconnc==2 
		   && ICCsizes[1]!=5 && ICCsizes[0]!=5
		   && ICCsizes[1]!=4 && ICCsizes[0]!=4
		   //&& ICCsizes[1]!=3 && ICCsizes[0]!=3
		   //&& ICCsizes[1]!=2 && ICCsizes[0]!=2
		   //&& ICCsizes[1]!=7 && ICCsizes[0]!=7
		   //&& ICCsizes[1]!=8 && ICCsizes[0]!=8
	   ) return false;	
	
	//if (nrconnc==1 && !iscenterlengthzero(g,n, hv, hvc)) return false;
	
	if (nrconnc==1)
	{
		getgrlengths(g,n,hv,hvc,armlengths);
		//for (int i=0;i<5;i++) cerr << armlengths[i] << " ";
		//cerr << endl;
		if (armlengths[4] > 2)	return false;
		//if (armlengths[0]!=0 &&  armlengths[1]!=0) return false;
		//if (armlengths[2]!=0 &&  armlengths[3]!=0) return false;
		
		//if (armlengths[0]+armlengths[1]+armlengths[2]+armlengths[3]<8) return false;
	}
	
	//if (nrconnc==2 && (hvc[1]>n-5 || hvc[0]>n-5) ) return false;
//	if (nrconnc==2 && (hvc[1]==n-4 || hvc[0]==n-4) ) return false; // no
//	if (nrconnc==2 && (hvc[1]==n-3 || hvc[0]==n-3) ) return false;
//	if (nrconnc==2 && (hvc[1]==n-2 || hvc[0]==n-2) ) return false; //noch ok, bringt 10
	
	
	//if (nrconnc==2) return false; // no
	
	//if (nrconnc==1 && (hvc[1]==n-3 || hvc[0]==n-3) ) return false; // no
	//if (nrconnc==1 && (hvc[1]==n-3 || hvc[0]==n-3) ) return false; //no
	
	
	
	
	
//	if (nrconnc!=2) return false;
//	if (hvc[1]!=5 && hvc[0]!=5) return false;
//	if (ICCsizes[1]!=7 && ICCsizes[0]!=7) return false;
	
	
	return true;

}

int main(int argc, char* argv[])
{
	
	graphmap gm;
	gm.clear();
	readgraphmap("gra13_0p2.g6", gm);
	
	graphmap::const_iterator it;
	for (it=gm.begin();it!=gm.end();it++)
	{
		graph g[30];
		int n= stringtogra(it->first, g);
		
		//for (int i=0;i<n;i++) {
		//	int ecount =0;
		//	for (int j=0;j<n;j++) 
		//		if (i!=j && TISELEMENT(g[i],j)) 
		//			ecount++;
			//if (ecount >= n-4)
		//	{
				
				if (checkgra(g, n))
					cout << it->second  <<endl;
		//		break;
		//	}
		//}
	}
	
	return 0;
	// read in three files and check the Jacobi id. on them
	graphvect v1,v2,v3, v4;
	readgraphvect("wheel3.hg6", v1);
			//graph g[MAXN];
			//int n= stringtogra(v1.begin()->first, g);
			//cout << ordersymmgroup(g, n) << endl;
			//return 0;
	readgraphvect("wheel5.hg6", v2);
	cleanv(v1);
	cleanv(v2);
	dualbasis2coinvariants(v1);
	dualbasis2coinvariants(v2);
	
	printgraphvect(stdout, v1);
	printgraphvect(stdout, v2);
	
	v3 = bracket(v1,v2);
	v4 = bracket(v1,v3);
	cleanv(v4);
	coinvariants2dualbasis(v4);
	//printgraphvect(stdout, v3);
	cleanv(v4);
	writegraphvect("wheelbr3_3_7.hg6", v4);

	return 0;
}