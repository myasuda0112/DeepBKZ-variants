#ifndef DeepBKZ_H
#define DeepBKZ_H

#include "lattice.h"
#include "DeepLLL.h"
#include "Enum.h"
#include <time.h>

/**********
DeepBKZ 
**********/
inline bool lattice::DeepBKZ(int start, int end, int b, FLOAT alpha, int gamma, int abort)
{
	bool res;
	int z, j, k, h, i, m, n = NumRows, l, tour = 0, N;
	double current; 
	long double rho, rho1;
	Vec<int> v; v.SetLength(n);
	Vec<double> R; R.SetLength(n);
	Vec<UINT> sv; sv.SetLength(n);
	mat_ZZ Basis;
	clock_t time_start = clock(), time_end;
	clock_t time_start1, time_end1;
	double total, time = 0.0; 
		
	cout << "\nDeepBKZ-" << b << endl;
	SetGSO();
	res = DeepLLL(start, end, start+1, alpha, gamma); /* Initial DeepLLL */
	if (res != true) {
		cout << "DeepLLL error1 in DeepBKZ" << endl;
		return false; 
	}
	current = pow(B(start), 0.5);
	cout << "Initial norm = " << current << endl;
	
	/* Main loop */
	GSA_slope(rho, start, end); 
	tour = 0; N = 0; 
	z = start-1; j = start-1;
	while( z < end-1 ){
		if (j == end-1) {
			j = start-1; ++tour;
			if (tour % 10 == 0) {
				cout << "tour = " << tour << ": GSA slope = " << rho << endl;
			}
			
			GSA_slope(rho1, start, end);
			/* Early termination */
			if (abort != 0) {
				if (rho > rho1) {
					++N;
					if (N >= abort) {
						cout << "Auto-abort termination" << endl;
						time_end = clock(); 
						total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
						cout << "Total time (seconds) = " << total << endl;
						cout << "DeepLLL time (seconds) = " << time << endl; 
						return true;
					}
				} else {
					N = 0; 
				}
				rho = rho1; 
			}
		}
		++j; k = min(j+b-1, end); h = min(k+1, end);
		
		/* Enumerate projected shortest vector */
		for (i=1; i<=n; ++i) { R(i) = (double) alpha*B(j); }
		res = ENUM(v, R, j, k);
		
		time_start1 = clock(); 
		if( res != true ) {
			++z;
			res = DeepLLL(start, h, h-1, alpha, gamma); /* DeepLLL */
			if (res != true) {
				cout << "DeepLLL error-1 in DeepBKZ" << endl; 
				return false; 
			}
		} else {
			z = start-1;
			
			/* Construction of a short lattice vector */
			for (l=1; l<=n; ++l) { sv(l) = 0; }
			for (i=j; i<=k; ++i) {
				for (l=1; l<=n; ++l) { sv(l) += v(i)*basis(i, l);}  
			}
			
			/* Insertion of a short lattice vector */
			Basis.SetDims(h+1, n);
			for (i=1; i<=j-1; ++i) {
				for (l=1; l<=n; ++l) {Basis(i, l) = basis(i, l);}
			}
			for (l=1; l<=n; ++l) { Basis(j, l) = sv(l);}
			for (i=j+1; i<=h+1; ++i) { 
				for (l=1; l<=n; ++l) {Basis(i, l) = basis(i-1, l);}
			}
			
			/* Modified LLL to remove the linear dependency  */
			NTL::LLL_FP(Basis, alpha, 0, 0, 0);
			
			/* Conversion from NTL::mat_ZZ to UINT */ 
			for (i=1; i<=h; ++i) {
				for (l=1; l<=n; ++l) {basis(i, l) = to_int(Basis(i+1, l));}
			}
			
			/* DeepLLL after insertion of short lattice vector */
			SetGSO(); 
			res = DeepLLL(start, h, j, alpha, gamma);
			if (res != true) {
				cout << "DeepLLL error-2 in DeepBKZ" << endl; 
				return false; 
			}
		}
		time_end1 = clock();
		time += static_cast<double>(time_end1 - time_start1)/CLOCKS_PER_SEC;
		
		/* for debug */
		if (B(start) < alpha*current*current) {
			current = pow(B(start), 0.5); 
			cout << "Norm = " << current << endl;
			cout << basis(start) << endl;
		}	
	}
	
	time_end = clock(); 
	total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
	cout << "Total time (seconds) = " << total << endl;
	cout << "DeepLLL time (seconds) = " << time << endl;
	return true;
}



#endif // DeepBKZ_H
