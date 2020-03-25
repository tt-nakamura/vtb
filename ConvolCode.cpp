// uses NTL
//   http://www.shoup.net/ntl

#include "ConvolCode.h"

#define CONVOL_DEGMAX 8

const long ConvolCode::coeff[2] = {0171, 0133};

ConvolCode::ConvolCode(const long *c, long N)
// c = array of coefficients of generator polynomials
// N = number of generator polynomials
{
    long i,j;
    g.SetLength(N);
    for(i=m=0; i<N; i++) {
        GF2XFromLong(g[i],c[i]);
        if(deg(g[i]) > m) m = deg(g[i]);
    }
    if(m > CONVOL_DEGMAX) Error("degree too large");
    C.SetDims(N,1<<(m+1));
    for(i=0; i<N; i++)
        for(j=1; j<C.NumCols(); j++)
            if(weight(j&c[i])&1) set(C[i][j]);
}

void ConvolCode::encode(Vec<GF2>& t, const GF2X& s)
// t = encoding of s
{
    long i,j,k;
    GF2X h;
    t.SetLength((deg(s)+m+1)*g.length());
    clear(t);
    for(i=0; i<g.length(); i++) {
        mul(h,s,g[i]);
        for(j=0, k=i; j<=deg(h); j++, k+=g.length())
            t[k] = h[j];// interleave
    }
}

void ConvolCode::encode(Vec<GF2>& t, const std::string& s) {
    GF2X f;
    GF2XFromBytes(f, (const unsigned char*)s.c_str(), s.length());
    encode(t,f);
}

void ConvolCode::decode(GF2X& s, const Vec<GF2>& t)
// s = decoding of t by Viterbi method
{
    long i,j,k,l,h,a,b,L(t.length()/g.length()),M(1<<m);
    Mat<long> A,B;
    A.SetDims(L+1,M);
    B.SetDims(L,M);
    A[0][0] = 0;
    for(j=1; j<M; j++) A[0][j] = -t.length();
    for(i=k=0; i<L; i++, k+=g.length()) {
        for(j=0; j<M; j++) {
            h = j|M;
            a = A[i][j>>1];
            b = A[i][h>>1];
            for(l=0; l<g.length(); l++) {
                if(C[l][j] == t[k+l]) a++;
                if(C[l][h] == t[k+l]) b++;
            }
            if(a>b) { A[i+1][j] = a; B[i][j] = j>>1; }
            else    { A[i+1][j] = b; B[i][j] = h>>1; }
        }
    }
    clear(s);
    for(i=L-1, j=0; i>=0; j=B[i--][j])
        if(j&1) SetCoeff(s,i);
}

void ConvolCode::decode(std::string& s, const Vec<GF2>& t) {
    GF2X f;
    decode(f,t);
    s.resize(NumBytes(f));
    BytesFromGF2X((unsigned char*)s.data(), f, s.length());
}