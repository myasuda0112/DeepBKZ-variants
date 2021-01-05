#include "lattice.h"

/*************************************
Efficient GSO update in Dual DeepLLLL
*************************************/
inline bool lattice::DualDeepLLLGSOUpdate(Vec<FLOAT>& inverseB, int k, int l, Vec<FLOAT> hmu, Vec<FLOAT> HD)
{
    int i, j, h, n = NumRows;
    Mat<FLOAT> xi; xi.SetDims(n, n);
    FLOAT tmp, tmp1, tmp2, tmp3;
    
    for (i=1; i<=n; ++i) {
        for (j=1; j<=n; ++j) {xi(i, j) = mu(i, j);}
    }

    for (i=l+1; i<=n; ++i) {
        tmp = 0.0;
        for (h=k; h<=l; ++h) {tmp += hmu(h)*mu(i, h);}
        xi(i, l) = tmp;
    }
    for (j=k; j<=l-1; ++j) {
        tmp2 = 1.0/(HD(j+1)*B(j+1));
        tmp3 = HD(j)*B(j+1);
        
        for (i=j+1; i<=l-1; ++i) {
            tmp = mu(i+1, j+1)*tmp3;
            tmp1 = 0.0;
            for (h=k; h<=j; ++h) {tmp1 += hmu(h)*mu(i+1, h);}
            tmp1 = tmp - (tmp1*hmu(j+1));
            xi(i, j) = tmp1*tmp2;
        }
        
        xi(l, j) = -hmu(j+1)*tmp2;
        
        for (i=l+1; i<=n; ++i) {
            tmp = mu(i, j+1)*tmp3;
            tmp1 = 0.0;
            for (h=k; h<=j; ++h) {tmp1 += hmu(h)*mu(i, h);}
            tmp1 = tmp - (tmp1*hmu(j+1));
            xi(i, j) = tmp1*tmp2;
        }
    }
    for (j=1; j<=k-1; ++j) {
        for (i=k; i<=l-1; ++i) {xi(i, j) = mu(i+1, j);}
        xi(l, j) = mu(k, j);
    }
    
    for (i=1; i<=n; ++i) {
        for (j=1; j<=n; ++j) {mu(i, j) = xi(i, j);}
    }

    for (j=k; j<=l-1; ++j) {B(j) = (HD(j+1)*B(j+1))/HD(j);}
    B(l) = 1.0/(HD(l));
    for (j=k; j<=l; ++j) {inverseB(j) = 1.0/B(j);}
    
    return true; 
}

/*****************************************************
DualDeepLLL (with deep insertion restriction) 
- Input
	FLOAT alpha		: reduction parameter
	int end			: end index
	int start		: start index
	int hh			: deep start index
	int gamma		: insertion restriction parameter
- Output 
Dual DeepLLL-reduced basis
******************************************************/
inline bool lattice::DualDeepLLL(FLOAT alpha, int end, int start, int hh, int gamma)
{
    int i, j, h, k, l, flag, n = NumRows;
	UINT q; 
	FLOAT tmp, D;
    bool res;
    Vec<UINT> vec_tmp; vec_tmp.SetLength(n);
    Vec<FLOAT> hmu_tmp, HD; hmu_tmp.SetLength(n), HD.SetLength(n);
    Vec<FLOAT> inverseB; inverseB.SetLength(n);
	Mat<FLOAT> hmu; hmu.SetDims(n, n);
    
    res = SetGSO();
    if (res != true) {
        cout << "SetGSO error in DualDeepLLL" << endl;
        return false;
    }
    for (i=1; i<= n; ++i) {inverseB(i) = 1.0/B(i);}
    
    for (i=1; i<=n; ++i) { hmu(i, i) = 1.0; }
	k = hh; 
    while(k >= end) {
		// cout << "k = " << k << endl;
        hmu(k, k) = 1.0;
        for (j=k+1; j<=n; ++j) {
            /* Compute of hmu(k, j) */
            tmp = 0.0;
            for (h = k; h <= j-1; ++h) {tmp += mu(j, h)*hmu(k, h);}
            hmu(k, j) = -tmp;
            
			/* Size reduction at k */
			if (fabs(hmu(k, j)) > ETA) {
				q = (UINT) roundl(hmu(k, j));
				
				/* Update of basis vector */
				for (l=1; l<=n; ++l) { basis(j, l) += q*basis(k, l);} 
				
				/* Update of hmu and mu */
                for (l=j; l<=n; ++l) {hmu(k, l) -= q*hmu(j, l);}
                for (h=1; h<=k; ++h) {mu(j, h) += q*mu(k, h);}
            }
        }
        
        flag = 0;
        for (i=k; i<=n; ++i) {HD(i)=0.0;}
        D = inverseB(k);
        HD(k) = D;
        for (h=k+1; h<=n; ++h) {
            D +=hmu(k, h)*hmu(k, h)*inverseB(h);
            HD(h) = D;
            if (B(h)*D < alpha && n-gamma<= h <= k+gamma) {
                l=h; ++flag;
            }
        }
        
        if (flag != 0) {
			/* Dual deep insertion at (k, l) */
			vec_tmp = basis(k);
			for (h = k; h != l; ++h) { basis(h) = basis(h+1); }
			basis(l) = vec_tmp;
			
			/* GSO update */
            for (h=1; h<=n; ++h) {hmu_tmp(h) = hmu(k, h);};
            res = DualDeepLLLGSOUpdate(inverseB, k, l, hmu_tmp, HD);
            k = min(l, n-1);
        } else {
            --k;
        }
    }
	
	/* Size reduction for stability */
    for(k=end;k>=1;--k){
        hmu(k, k) = 1.0;
        for (j=k+1; j<=n; ++j) {
            
            tmp = 0.0;
            for (h = k; h <= j-1; ++h) {tmp += mu(j, h)*hmu(k, h);}
            hmu(k, j) = -tmp;
			
            /* Size reduction at k */
			if (fabs(hmu(k, j)) > ETA) {
				q = (UINT) roundl(hmu(k, j));
				
				/* Update of basis vector */
				for (l=1; l<=n; ++l) { basis(j, l) += q*basis(k, l);} 
				
				/* Update of hmu and mu */
                for (l=j; l<=n; ++l) {hmu(k, l) -= q*hmu(j, l);}
                for (h=1; h<=k; ++h) {mu(j, h) += q*mu(k, h);}
            }
        }   
    }
    
    return true;
}










