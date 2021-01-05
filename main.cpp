#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
#include <NTL/LLL.h>

#include "lattice.h"
#include "DeepLLL.h"
#include "DeepBKZ.h"
#include "DualDeepLLL.h"
#include "DualDeepBKZ.h"
#include "SelfDualDeepBKZ.h"


int main(){
	int i, j, k, n=120, b, f;
	bool res; 
	lattice L; L.SetDims(n, n);
	
	/* input file */
	std::ifstream File("BKZ20_SVP120_seed0.txt");
	File >> L.basis;
	
	/* Compute Gram-Schmidt information */
	// L.SetGSO();

	
#if 1
	int ff = 26; 
	for (i=0; i<=0; ++i) {
		b = 40; 
		// res =L.DeepBKZ(1, n, b, 0.99, n, 5);
		// res = L.DualDeepBKZ(b, 0.99, 1, n, n, 5);
		// f = floor(1.0*ff*b/(n-ff+1));
		// f = round(1.0*ff*b/(n-ff+1));
		// cout << "f = " << f << endl; 
		// f = ceil(ff*b/(n-ff+1)); 
		// cout << "f = " << f << endl; 
		res = L.SelfDualDeepBKZ(b, 0.99, f, n, 20);
		// res = L.DeepBKZ(1, n, b, 0.99, n, 20);
		// res = L.DualDeepBKZ(b, 0.99, 1, n, n, 20);		
		if (res = false) {
			cout << "DeepBKZ error in main function" << endl; 
		}
	}
#endif 
	
	return 0;
}
