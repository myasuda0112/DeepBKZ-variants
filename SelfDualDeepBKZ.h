#include <sys/time.h>
#include "lattice.h"
#include "DeepLLL.h"
#include "Enum.h"
#include "DeepBKZ.h"
#include "DualEnum.h"
// #include "DualDeepBKZ.h"

bool DualENUM(Vec<int>& result, Mat<FLOAT> mu, Vec<double> InverseB, Vec<double> A); 

/***********************************
SelfDualDeepBKZ with free dimensions
***********************************/
inline bool lattice::SelfDualDeepBKZ(int b, long double alpha, int f, int gamma, int abort)
{
	bool res;
	long double current;
	long double slope, slope1; 
	long long N, tour;
	long double rho, rhof, rhob, rho1;
	long insertion = 0;
	int d, z = 0, j = 0, k = 0, h = 0, hh, i, m=NumCols, n = NumRows, l, ff; 
	int flag, r, fflag, num;
	Vec<int> v, w, tmp1, tmp2; v.SetLength(n); w.SetLength(n); tmp1.SetLength(n); tmp2.SetLength(n);   
	Vec<long int> sv, tmph;  sv.SetLength(n); tmph.SetLength(n);
	Vec<double> R; R.SetLength(n);
	Vec<double> local_B, A;
	Mat<long double> local_mu; 
	Vec<int> x, zz;
	Mat<UINT> dBasis, NewBasis;
	mat_ZZ Basis; 
	
	local_mu.SetDims(b, b); local_B.SetLength(b); 
	A.SetLength(b); x.SetLength(b);
	
	res = SetGSO();
	if (res != true) {
		cout << "Error: SetGSO in SelfDualDeepBKZ" << endl; 
		return false; 
	}

	/* Initial DeepLLL */
	res = DeepLLL(1, n, 2, alpha, gamma);
	if (res != true) {
		cout << "Error: DeepLLL in SelfDualDeepBKZ" << endl; 
		return false; 
	}

	current = pow(B(1), 0.5); 
	cout << "\nSelfDualDeepBKZ-" << b << endl;
	cout << "Norm = " << current << endl;	
	
	N = 0; GSA_slope(rho1, 1, n);
	tour = 0; 
	while(1) {
		Basis.SetDims(n, n);
		for (i=1; i<=n; ++i) {
			for (l=1; l<=n; ++l) {
				Basis(i, l) = basis(i, l); 
			}
		}
		NTL::LLL_FP(Basis, 0.999, 0, 0, 0);			
		/* Conversion from NTL::mat_ZZ to UINT */ 
		for (i=1; i<=n; ++i) {
			for (l=1; l<=n; ++l) {basis(i, l) = to_int(Basis(i, l));}
		}
		SetGSO();
		
		++tour; 
		/* Computation of GSA slopes */
		GSA_slope(rho, 1, n);
		GSA_slope(rhof, 1, n/2); 
		GSA_slope(rhob, n/2+1, n);
		
		if (tour % 10 == 0) {
			cout << "tour = " << tour << ", rho = " << rho << endl;
		// cout << "rhof = " << rhof << endl;
		// cout << "rhob = " << rhob << endl;
		}
		
		if (rho > rho1) {
			N = 0;
		} else {
			++N; 
			if (N >= abort) {
				cout << "Auto-abort" << endl; 
				return true; 
			}
		}
		rho1 = rho;
		
		if (rhof <= rhob) { /* forward tours */
			// cout << "Forward tours" << endl;
			flag = 0; 
			for (j=1; j<=n-f; ++j) {
				k = min(j+b-1, n-f+1); h = min(k+1, n-f+1);
				
				/* Enumerate projected shortest vector */
				for (i=1; i<=n; ++i) { R(i) = (double) alpha*B(j); }
				res = ENUM(v, R, j, k);
				
				if( res != true ) {
					res = DeepLLL(1, h, h-1, alpha, gamma); /* DeepLLL */
					if (res != true) {
						cout << "DeepLLL error-1 in DeepBKZ" << endl; 
						return false; 
					}
				} else {
					++flag;
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
					res = DeepLLL(1, h, j, alpha, gamma);
					if (res != true) {
						cout << "DeepLLL error-2 in DeepBKZ" << endl; 
						return false; 
					}
				}
			}
			
			if (flag == 0) {
			 	return true; 
			}
			/* for debug */
			if (B(1) < alpha*current*current) {
				current = pow(B(1), 0.5); 
				cout << "Norm = " << current << endl;
				cout << basis(1) << endl;
				++fflag; 
			}	
		} else { /* backward tours */
			// cout << "Backward tours" << endl;
			flag = 0; 
			for (j = n; j>= f+1; --j) {
				k = max(j-b+1, f);
				d = j-k+1; /* d = dimension of a local projected lattice */
				hh = k-1; 
				
				/* Construction of a local projected basis */
				local_mu.SetDims(d, d); local_B.SetLength(d); 
				A.SetLength(d); x.SetLength(d);
				for (i=1; i<=d; ++i) {
					local_B(i) = 1.0/B(k+i-1);
					/* Full Enumeration Setting */
					A(i) = alpha/B(j);
					x(i) = 0;
				}
				for (i=1; i<=d; ++i) {
					for (h=1; h<=d; ++h) {local_mu(i, h) = mu(i+k-1, h+k-1);}
				}
					
				/* Dual Enumeration */
				res = DualENUM(x, local_mu, local_B, A);
				
				if (res == false) {
					res = DualDeepLLL(alpha, hh, n, hh+1, gamma); 
					if (res != true) {
						cout << "DualDeepLLL error-2 in Self-dual DeepBKZ" << endl; 
						return false;
					}
				} else {
					++flag;
					ff = 0; 
					for (i=1; i<=d; ++i) {
						if (x(i) != 0) { ++ff; }
					}
					if (ff == 0) {
						cout << "DualEnum output error in Self_dual DeepBKZ" << endl; 
						return false; 
					}
					
					// cout << "j=" << j << ", x=" << x << endl; /* for debug */
					
					/* Insertion of a short dual vector */
					NewBasis.SetDims(d, m);
					dBasis.SetDims(d, m); 
					for (i=1; i<=d; ++i) {
						for (h=1; h<=m; ++h) {dBasis(i, h) = basis(k+i-1, h);}
					}
					res = Insert(NewBasis, dBasis, x);
					if (res != true) {
						cout << "Insert error in Dual_DeepBKZ" << endl; 
						return false;
					}
					for (i=1; i<=d; ++i) {
						for (h=1; h<=m; ++h) {basis(k+i-1, h) = NewBasis(i, h);}
					}
					
					/* Dual DeepLLL after insertion */
					// res = DualDeepLLL(alpha, end, start, start-1, gamma);
					res = DualDeepLLL(alpha, hh, n, j, gamma);
					if (res != true) {
						cout << "DualDeepLLL error-2 in Dual_DeepBKZ" << endl; 
						return false;
					}
				}
			}
			if (flag == 0) {
			 	return true; 
			}
		}
	}
	return true; 
}
