#include "lattice.h"
#include "DualEnum.h"


/************************
Matrix Multiplication
************************/
inline void MatMul(Mat<UINT>& C, Mat<UINT> A, Mat<UINT> B)
{
	int i, j, k, nn=A.NumRows(), mm=A.NumCols(), ll=B.NumCols();	
	for (i=1; i<=nn; ++i) {
		for (j=1; j<=ll; ++j) {
			C(i, j) = 0; 
			for (k=1; k<=mm; ++k) {C(i, j) += A(i, k)*B(k, j);}
		}
	}
}


/****************************************************************************************
Insertion of a new vector to obtain a new (primal) basis
- Input
	Mat<long long>& NewBasis	: New basis
	Mat<long long> Basis		: Current Basis
	Vec<int> x					: coefficient vectors of a short dual lattice vector
- Output
	true(1)			: successful
	false(0)		: program error occurs
- Content
	x is inserted into the dual basis to obtain a new basis
- Reference	
	[1] Algorithm 2 in "Practical, Predictable Lattice Basis Reduction", 
	D. Micciancio and M. Walter (presented at EUROCRYPT2016)
*****************************************************************************************/
inline bool Insert(Mat<UINT>& NewBasis, Mat<UINT> Basis, Vec<int> x)
{
	int i, j, k, mm = Basis.NumCols(), nn = Basis.NumRows();
	double alpha = 0.99;
	double beta = 4.0/(4*alpha - 1); 
	double tmp, gamma;
	bool res;
	Mat<UINT> U; U.SetDims(nn, nn);
	mat_ZZ tmp_base; tmp_base.SetDims(nn, nn+1); 
	 
	/* Construction of gamma */
	tmp = 0.0; for (i=1; i<= nn; ++i) {tmp += x(i)*x(i);}
	tmp = pow(tmp, 0.5); /* Euclid norm of x */
	tmp = tmp*pow(beta, (nn-2)/2.0);
	gamma = roundl(2*tmp);
	
	/* Construction of matrix */
	for (i=1; i<=nn; ++i) {for (j=1; j<=nn; ++j) {tmp_base(i, j) = 0;}}
	for (i=1; i<=nn; ++i) {tmp_base(i, i) = 1;}
	for (i=1; i<=nn; ++i) {tmp_base(i, nn+1) = gamma*x(i);} 
	
	/* LLL-reduction (NTL-function) */
	NTL::LLL_FP(tmp_base, alpha, 0, 0, 0); 
	
	/* Extraction of a unimodular matrix */
	for (i=1; i<=nn; ++i) {
		for (j=1; j<=nn; ++j) {U(i, j) = to_int(tmp_base(i, j));}
	}

	/* Construction of a new basis */
	MatMul(NewBasis, U, Basis);
	
	return true; 
}

/*******************************************************************
Dual DeepBKZ (by using fast dual DeepLLL and Enumeration)
- Input
	lattice class	: lattice basis information
	int beta		: block size
	double			: reduction parameter
- Output
	true(1)			: successful
	false(0)		: program error occurs
- Content
	The dual of the output basis is alpha-DeepBKZ-reduced 
*******************************************************************/
inline bool lattice::DualDeepBKZ(int beta, FLOAT alpha, int end, int start, int gamma, int abort)
{
	int i, j, k, d, h, f, flag, m = NumCols, nn = NumRows, count, zz, hh, N, tour; 
	bool res;
	Vec<double> local_B, A;
	Mat<FLOAT> local_mu;
	Vec<int> x;
	Mat<UINT> Basis, NewBasis;
	mat_ZZ tmp_base; tmp_base.SetDims(nn, m);
	long double rho, rho1; 
		
	cout << "DualDeepBKZ-" << beta << endl;	
	/* Initial Dual DeepLLL */
	res = DualDeepLLL(alpha, end, start, start-1, gamma);
	if (res != true) {
		cout << "DualDeepLLL Error-1 in Dual_DeepBKZ" << endl; 
		return false;
	}
	
	N = 0;
	GSA_slope(rho, start, end);
	zz = start+1; j = start+1;
	tour = 0; 
	while (zz > end-1) {
		if (j==end+1) {
			++tour; 
			j = start+1;
			
			GSA_slope(rho1, end, start);
			if (tour % 10 == 0) {
				cout << "tour = " << tour << ", rho = " << rho1 << endl;
			}
			/* Early termination */
			if (abort != 0) {
				if (rho > rho1) {
					++N;
					if (N >= abort) {
						cout << "Auto-abort termination" << endl; 
						return true; 
					}
				} else {
					N = 0; 
				}
				rho = rho1; 
			}
		}
		--j;
		k = max(j-beta+1, end);
		hh = max(k-1, end); 
		d = j-k+1; /* d = dimension of a local projected lattice */
		
		/* Construction of a local projected basis */
		local_mu.SetDims(d, d); local_B.SetLength(d); 
		A.SetLength(d); x.SetLength(d);
		for (i=1; i<=d; ++i) {
			local_B(i) = static_cast<double>(1.0/B(k+i-1));
			/* Full Enumeration Setting */
			A(i) = static_cast<double>(alpha/B(j));
			x(i) = 0;
		}
		for (i=1; i<=d; ++i) {
			for (h=1; h<=d; ++h) {local_mu(i, h) = mu(i+k-1, h+k-1);}
		}
			
		/* Dual Enumeration */ 
		DualENUM(x, local_mu, local_B, A);
		f = 0; 
		for (i=1; i<=d; ++i) {
			if (x(i) != 0) {++f;}
		}
		
		if (f == 0) {
			--zz; 
			res = DualDeepLLL(alpha, hh, start, hh+1, gamma); 
			if (res != true) {
				cout << "DualDeepLLL error-2 in Dual_DeepBKZ" << endl; 
				return false;
			}
		} else {
			zz = start-1; 
			// cout << "j=" << j << ", x=" << x << endl; /* for debug */
			
			/* Insertion of a short dual vector */
			NewBasis.SetDims(d, m);
			Basis.SetDims(d, m); 
			for (i=1; i<=d; ++i) {
				for (h=1; h<=m; ++h) {Basis(i, h) = basis(k+i-1, h);}
			}
			res = Insert(NewBasis, Basis, x);
			if (res != true) {
				cout << "Insert error in Dual_DeepBKZ" << endl; 
				return false;
			}
			for (i=1; i<=d; ++i) {
				for (h=1; h<=m; ++h) {basis(k+i-1, h) = NewBasis(i, h);}
			}
			
			/* Dual DeepLLL after insertion */
			// res = DualDeepLLL(alpha, end, start, start-1, gamma);
			res = DualDeepLLL(alpha, hh, start, j, gamma); 
			if (res != true) {
				cout << "DualDeepLLL error-2 in Dual_DeepBKZ" << endl; 
				return false;
			}
		}
	}
	return true;
}

