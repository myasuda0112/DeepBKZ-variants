#ifndef DualENUM_H
#define DualENUM_H

#include <float.h>

#include "lattice.h"

#if 1
/***********************************************************************************
Dual Enumeration_fast (faster version by using double type) 
-Input/Ouput
	Vec<int>& result		: the coordinates of a short dual lattice vector
	Mat<long double> mu		: GSO coefficients
	Vec<double> InverseB	: Inverse of squared GSO lengths
	Vec<double> A			: upper bound of enumeration area
-Requirement
	Appropriate squared GSO lengths are given
-Content
	Enumeration of the coordinates of short dual lattice vectors with 
	squared length < A
-Reference
	[1] Algorithm 2 in "Practical, Predictable Lattice Basis Reduction", 
	D. Micciancio and M. Walter (presented at EUROCRYPT2016)
	[2] Algorithm 2 in "Lattice Enumeration using Extreme Pruning", 
	N. Gama, P.Q. Nguyen and O. Regev (presented at EUROCRYPT2010)
************************************************************************************/
inline bool DualENUM(Vec<int>& result, Mat<FLOAT> mu, Vec<double> InverseB, Vec<double> A)
{
	int i, j, h, k, flag=0, last_nonzero=1, d=result.length(), tmp0;
	Vec<int> v, w; v.SetLength(d), w.SetLength(d); 
	Mat<FLOAT> hmu; hmu.SetDims(d, d);
	Mat<double> hhmu; hhmu.SetDims(d, d);
	Vec<double> c, rho; rho.SetLength(d), c.SetLength(d);
	double tmp;
	Mat<double> sigma; sigma.SetDims(d + 1, d);
	Vec<int> r; r.SetLength(d);
	
	/* Compute of hmu(k, j) */
	for (k=1; k <= d; ++k) {
		for (j=1; j<=k-1; ++j) {hmu(k, j) = 0.0;}
		hmu(k, k) = 1.0; 
		for (j=k+1; j<=d; ++j) {
			tmp = 0.0; for (h = k; h <= j-1; ++h) {tmp += mu(j, h)*hmu(k, h);}
			hmu(k, j) = -tmp;
		}
	}
	for (k=1; k<=d; ++k) {
		for (j=1; j<=d; ++j) {
			hhmu(k, j) = hmu(k, j);
		}
	}

	/* Step 1 */
	for (i = 1; i <= d+1; ++i) { for (j = 1; j <= d; ++j) { sigma(i, j) = 0.0; }}
	for (i = 1; i <= d+1; ++i) { r(i) = i-1; }
	
	/* Initialization */
	v(1) = 1;
	for (i=2; i<=d; ++i) {v(i) = 0;} /* v = (1,0,...,0) */
	for (i=1; i<=d; ++i) {
		rho(i) = 0.0, c(i) = 0.0, w(i) = 0; /* c(i): centers & w(i): jumps */
	}
	
	k=1; flag = 0; 
	while(true){
		/* compute squared norm rho(k) of current node */
		tmp = v(k)-c(k); 
		tmp = tmp*tmp; 
		if (k>=2) {rho(k) = rho(k-1) + tmp*InverseB(k);}
		if (k==1) {rho(1) = tmp*InverseB(1);}
		
		if (rho(k) <= A(k)) {
			if (k==d) {
				if (abs(rho(d)) < DBL_EPSILON) {
					if (flag == 0) {
						for (i=1; i<=d; ++i) { result(i) = 0; }
						return false;
					} else {
						// cout << "flag = " << flag << endl; 
						return true; 
					}
				} else {
					/* Solution is found */
					for (i=1; i<=d; ++i) { 
						result(i) = v(i);
						A(i) = min(A(i), 0.99*rho(d)); 
					}
					++flag;
				} 
			} else {
				++k; /* going down the tree */
				tmp0 = r(k); 
				r(k+1) = min(r(k+1), tmp0);
				for (i=tmp0; i<=k-1; ++i) {
					sigma(i+2, k) = sigma(i+1, k) + v(i)*hhmu(i, k); 
				}
				c(k) = -sigma(k+1, k);
				v(k) = (int) round(c(k));
				w(k) = 1;
			}
		} else {
			--k; /* going up the tree */
			if (k==0) {
				if (flag != 0) {
					// cout << "flag = " << flag << endl; 
					return true;
				} else {
					for (i=1; i<=d; ++i) { result(i) = 0; } /* No solution is found */
					return false;
				}
			}
			r(k+1) = k;
			if (k >= last_nonzero) {
				last_nonzero = k;
				++v(k);
			} else {
				if (v(k) > c(k)) {v(k) -= w(k);} 
				else {v(k) += w(k);}
				++w(k);
			}
		}
	}
}
#endif 


#if 0
/***********************************************************************************
Dual Enumeration_fast (faster version by using double type) 
-Input/Ouput
	Vec<int>& result		: the coordinates of a short dual lattice vector
	Mat<long double> mu		: GSO coefficients
	Vec<double> InverseB	: Inverse of squared GSO lengths
	Vec<double> A			: upper bound of enumeration area
-Requirement
	Appropriate squared GSO lengths are given
-Content
	Enumeration of the coordinates of short dual lattice vectors with 
	squared length < A
-Reference
	[1] Algorithm 2 in "Practical, Predictable Lattice Basis Reduction", 
	D. Micciancio and M. Walter (presented at EUROCRYPT2016)
	[2] Algorithm 2 in "Lattice Enumeration using Extreme Pruning", 
	N. Gama, P.Q. Nguyen and O. Regev (presented at EUROCRYPT2010)
************************************************************************************/
inline bool DualENUM(Vec<int>& result, Mat<long double> mu, Vec<double> InverseB, Vec<double> A)
{
	int i, j, h, k, flag, last_nonzero=1, d=result.length();
	Vec<int> v, w; v.SetLength(d), w.SetLength(d); 
	Mat<long double> hmu; hmu.SetDims(d, d);
	Mat<double> hhmu; hhmu.SetDims(d, d);
	Vec<double> c, rho; rho.SetLength(d), c.SetLength(d);
	double tmp; 
	
	/* Compute of hmu(k, j) */
	for (k=1; k <= d; ++k) {
		for (j=1; j<=k-1; ++j) {hmu(k, j) = 0.0;}
		hmu(k, k) = 1.0; 
		for (j=k+1; j<=d; ++j) {	
			tmp = 0.0; for (h = k; h <= j-1; ++h) {tmp += mu(j, h)*hmu(k, h);}
			hmu(k, j) = -tmp;
		}
	}
	for (k=1; k<=d; ++k) {
		for (j=1; j<=d; ++j) {hhmu(k, j) = hmu(k, j);}
	}

	/* Initialization */
	v(1) = 1; for (i=2; i<=d; ++i) {v(i) = 0;} /* v = (1,0,...,0) */
	for (i=1; i<=d; ++i) {
		rho(i) = 0.0, c(i) = 0.0, w(i) = 0; /* c(i): centers & w(i): jumps */
	}
	
	k=1; flag = 0; 
	while(true){
		/* compute squared norm rho(k) of current node */
		tmp = v(k)-c(k); tmp = tmp*tmp; 
		if (k>=2) {rho(k) = rho(k-1) + tmp*InverseB(k);}
		if (k==1) {rho(1) = tmp*InverseB(1);}
		
		if (rho(k) < A(k)) {
			if (k==d) {
				/* Solution is found */
				for (i=1; i<=d; ++i) {result(i) = v(i);}
				for (i=1; i<=d; ++i) {A(i) = min(0.99*rho(d), A(i));}
				flag = 1; 
				cout << "v = " << result << endl; 
				// return true; /* Once solution is found, program ends */ 
			} else {
				k = k+1; /* going down the tree */
				tmp = 0.0; 
				for (i=1; i<k; ++i) {tmp += v(i)*hhmu(i, k);}
				c(k) = -tmp;
				v(k) = round(c(k)); 
				w(k)=1;
			}
		} else {
			k = k-1; /* going up the tree */
			if (k==0) {
				if (flag == 1) {
					return true;
				} else {
					/* There is no solution */
					for (i=1; i<=d; ++i) {result(i) = 0;}
					return false; /* No solution is found */
				}
			}
			if (k >= last_nonzero) {
				last_nonzero = k;
				v(k) = v(k)+1;
			} else {
				if (v(k) > c(k)) {v(k)=v(k)-w(k);} 
				else {v(k)=v(k)+w(k);}
				w(k) = w(k)+1;
			}
		}
	}
}
#endif




#endif //DualENUM_H
