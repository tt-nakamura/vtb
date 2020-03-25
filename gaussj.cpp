// uses NTL
//   http://www.shoup.net/ntl

template<class T>
long gaussj(Mat<T>& A)
// transform A into reduced row echelon form
// returns rank of A
{
    long i,j,k,l,m(A.NumRows()), n(A.NumCols());
    T t;
    for(k=l=0; k<m && l<n; l++) {
        for(i=k; i<m && IsZero(A[i][l]); i++);
        if(i==m) continue;
        if(i>k) swap(A[i], A[k]);
        inv(t, A[k][l]);
        for(j=l; j<n; j++) A[k][j] *= t;
        for(i=0; i<m; i++) {
            if(IsZero(A[i][l]) || i==k) continue;
            for(j=n-1; j>=l; j--) {
                mul(t, A[i][l], A[k][j]);
                A[i][j] -= t;
            }
        }
        k++;
    }
    return k;
}

template<class T>
long gaussj(Mat<T>& A, Vec<long>& col)
// transform left part of A into identity matrix
// using row and column operations and assuming m<=n
// col = indices of column permutation
// returns rank of A
{
    long i,j,k,m(A.NumRows()), n(A.NumCols());
    T t;
    col.SetLength(n);
    for(i=0; i<n; i++) col[i] = i;
    for(k=0; k<m; k++) {
        for(i=k; i<m && IsZero(A[i][k]); i++);
        if(i==m) {
            for(j=k+1; j<n && IsZero(A[k][j]); j++);
            if(j==n) break;
            for(i=0; i<m; i++) swap(A[i][j], A[i][k]);
            swap(col[j], col[k]);
        }
        else if(i>k) swap(A[i], A[k]);
        inv(t, A[k][k]);
        for(j=k; j<n; j++) A[k][j] *= t;
        for(i=0; i<m; i++) {
            if(IsZero(A[i][k]) || i==k) continue;
            for(j=n-1; j>=k; j--) {
                mul(t, A[i][k], A[k][j]);
                A[i][j] -= t;
            }
        }
    }
    return k;
}
