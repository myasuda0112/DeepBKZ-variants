#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <NTL/LLL.h>
#include <vector>
#include <math.h> 
#include <algorithm>
// #include "Arithmetic.h"

using namespace std;
using namespace NTL;

#if 0
/* fast version */
typedef int UINT;
typedef double FLOAT;
#endif 

#if 1
/* stable version */
typedef long int UINT;
typedef long double FLOAT;
#endif 

static const FLOAT ETA=0.501;

/**************************
Class for a lattice basis
**************************/
class lattice {
public:

	Mat<UINT> basis;	/* lattice basis */
	Mat<FLOAT> mu;		/* GSO coefficients */
	Vec<FLOAT> B;		/* squared lengths of GSO vectors */
	
	FLOAT tmp, tmp1, tmp2, tmp3;
	Vec<UINT> tmpp;
	Vec<FLOAT> D, DD, BB, mut, tmpld;

	int NumRows; /* Number of row vectors */
	int NumCols; /* Number of column vectors */
	
	/* lattice.h */
	void SetDims(int rows, int cols); /* Initialization */
	bool SetGSO(); 	  /* Gram-Schmidt orthogonalization */
	void GSA_slope(long double &slope, int start, int end); 
	
	/* DeepLLL.h */
	void size_reduce(int k);
	bool DeepLLL(int start, int end, int h, FLOAT alpha, int gamma);

	/* Enum.h */
	bool ENUM(Vec<int> &result, Vec<double> SR, int g, int h);
	
	/* DeepBKZ.h */
	bool DeepBKZ(int start, int end, int b, FLOAT alpha, int gamma, int abort);
	
	/* DualDeepLLL.h */
	bool DualDeepLLLGSOUpdate(Vec<FLOAT>& inverseB, int k, int l, Vec<FLOAT> hmu, Vec<FLOAT> HD); 
	bool DualDeepLLL(FLOAT alpha, int end, int start, int hh, int gamma);
	
	/* DualDeepBKZ.h */
	bool DualDeepBKZ(int beta, FLOAT alpha, int end, int start, int gamma, int abort); 
	
	/* Self-dual DeepBKZ.h */
	bool SelfDualDeepBKZ(int b, long double alpha, int f, int gamma, int abort); 
};

/**********************************
Initialization of a lattice basis
**********************************/
inline void lattice::SetDims(int rows, int cols)
{
	int n = rows; 
	NumRows = rows;
	NumCols = cols;
	
	basis.SetDims(rows, cols);
	mu.SetDims(rows, rows);
	B.SetLength(rows);
	
	tmpp.SetLength(NumCols);
	DD.SetLength(n);
	D.SetLength(n);
	BB.SetLength(n);
	mut.SetLength(n);
	tmpld.SetLength(n);
}

/***********************************
Inner product
***********************************/
inline bool InnerProduct(FLOAT &r, Vec<UINT> &a, Vec<UINT> &b)
{
	int i, n = a.length();
	if (n != b.length()) {
		cout << "Length error in InnerProduct" << endl; 
		return false; 
	}
	
	r = 0.0; 
	for (i=1; i<=n; ++i) { r += a(i)*b(i); }
	
	return true; 
}
/*****************************************
Computation of Gram-Schmidt information,  
based on Cholesky decomposition
******************************************/
inline bool lattice::SetGSO()
{
	int i, j, k, n = NumCols;
	bool res;
	Mat<FLOAT> r, s; r.SetDims(n, n), s.SetDims(n, n); 
	
	for (i=1; i <= n; ++i) {
		for (j=1; j < i; ++j) {
			res = InnerProduct(r(i, j), basis(i), basis(j));
			if (res == false) {
				cout << "InnerProduct error-1 in SetGSO" << endl; 
				return false; 
			}
			for (k=1; k< j; ++k) { r(i, j) -= mu(j, k)*r(i, k); }
			mu(i, j) = r(i, j)/r(j, j); 
		}
		mu(i, i) = 1.0; 
		
		res = InnerProduct(s(i, 1), basis(i), basis(i));
		if (res == false) {
			cout << "InnerProduct error-2 in SetGSO" << endl; 
			return false; 
		}
		for (j=2; j<=i; ++j) { s(i, j) = s(i, j-1) - mu(i, j-1)*r(i, j-1); }
		B(i) = s(i, i);
		r(i, i) = B(i);
	}
	return true; 
}

/*****************************
Computation of the GSA slope
******************************/
inline void lattice::GSA_slope(long double &slope, int start, int end)
{
	int i, mm = end-start+1;
	long double tmp1, tmp2, tmp3, tmp4; 
	
	slope = 0.0;
	tmp1 = tmp2 = 0.0; 
	for (i=start; i<=end; ++i) { 
		tmp1 += i*log(B(i));
		tmp2 += log(B(i));
	}
	tmp3 = tmp4 = 0.0; 
	for (i=start; i<=end; ++i) {
		tmp3 += i; 
		tmp4 += i*i;
	}
	slope = (mm*tmp1 - tmp3*tmp2)/(mm*tmp4 - tmp3*tmp3);
}

#endif // LATTICE_H
