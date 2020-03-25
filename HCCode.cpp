// uses NTL
//   http://www.shoup.net/ntl

#include "HCCode.h"
#include "gaussj.cpp"

#define HC_GF2E_DEGMAX 8

HCCode::HCCode(long n, long k)
// n = data length after encoding
// k = data length before encoding
// assume GF2E has been initialized
{
    long i,j,l,q2,g,s,m(GF2E::degree());
    if(m > HC_GF2E_DEGMAX) Error("GF2E degree too large");
    if(m&1) Error("m must be even");
    if(NumBits(n-1) > 3*(m>>1)) Error("n is too large");
    GF2EX h;
    GF2E x1,y1;
    GF2X f;
    Vec<GF2E> x,y,z;
    Vec<long> col;
    q = (1<<(m>>1));
    q2 = (1<<m);
    g = (q2-q)>>1;
    s = n-k;
    u = s>>1;
    n_err_max = u-g;
    t = n_err_max+1;
    if(n<=3*g+1) Error("n is too small");
    if(k<g)  Error("k is too small");
    if(t<=1) Error("k is too large");
    H.SetDims(n,s);
    a.SetLength(s);
    b.SetLength(s);
    x.SetLength(n);
    SetCoeff(h,q);
    SetCoeff(h,1);
    FindRoots(y,h);
    y.SetLength(n);
    SetCoeff(h,0);
    FindRoots(z,h);
    for(i=q, j=1; i<n; j++) {
        GF2XFromLong(f,j);
        conv(x1,f);
        power(y1, x1, q+1);
        for(l=0; l<q && i<n; l++, i++) {
            x[i] = x1;
            mul(y[i], y1, z[l]);
        }
    }
    K.SetDims(q+1, (s+g-1)/(q+1));
    for(i=j=l=0; l<s; l++) {
        a[l] = i;
        b[l] = j;
        K[i][j] = l;
        if(i>0) { i--; j++; }
        else if(j>=q) { i=q; j-=q-1; }
        else { i=j+1; j=0; }
    }
    for(i=0; i<n; i++) {
        for(j=0; j<s; j++) {
            power(x1, x[i], a[j]);
            power(y1, y[i], b[j]);
            mul(H[i][j], x1, y1);
        }
    }
    kernel(G,H);
    gaussj(G,col);
    Mat<GF2E> H1(H);
    for(i=0; i<n; i++)
        if(i!=col[i]) H[i] = H1[col[i]];
}

void HCCode::encode(Vec<GF2E>& c, const Vec<GF2E>& b)
// c = encoding of b
{
    long i,j;
    GF2E t;
    VectorCopy(c,b,length());
    for(j=dim(); j<length(); j++) {
        for(i=0; i<dim(); i++) {
            mul(t, b[i], G[i][j]);
            c[j] += t;
        }
    }
}

void HCCode::encode(Vec<GF2>& t, const GF2X& s)
// t = encoding of s (divided into data of length k)
// length of t is multiple of n
{
    long i,j,m(GF2E::degree());
    GF2X f(s),g;
    Vec<GF2> v;
    Vec<GF2E> b,c;
    t.SetLength(0);
    b.SetLength(dim());
    for(i=0; !IsZero(f);) {
        trunc(g,f,m);
        conv(b[i++], g);
        f >>= m;
        if(!IsZero(f) && i < b.length()) continue;
        while(i < b.length()) clear(b[i++]);
        encode(c,b);
        for(i=j=0; j < c.length(); j++) {
            conv(g,c[j]);
            VectorCopy(v,g,m);
            t.append(v);
        }
    }
}

void HCCode::encode(Vec<GF2>& t, const std::string& s) {
    GF2X f;
    GF2XFromBytes(f, (const unsigned char*)s.c_str(), s.length());
    encode(t,f);
}

long HCCode::decode(Vec<GF2E>& c, const Vec<GF2E>& d)
// c = decoding of d by Skorobogatov-Vladut method
// returns the number of corrected characters
// returns -1 if unable to correct
{
    long i,j,k,m,n;
    GF2E z,w;
    Vec<long> l;
    Vec<GF2E> s;
    Mat<GF2E> S,T;
    VectorCopy(c,d,dim());
    mul(s,d,H);
    if(IsZero(s)) return 0;
    S.SetDims(t,u);
    for(i=0; i<t; i++) {
        for(j=i; j<u; j++) {
            m = a[i] + a[j];
            n = b[i] + b[j];
            if(m<=q) S[i][j] = s[K[m][n]];
            else {
                m -= q+1;
                S[i][j]  = s[K[m][n+q]];                
                S[i][j] -= s[K[m][n+1]];
            }
        }
    }
    for(i=1; i<t; i++)
        for(j=0; j<i; j++) S[i][j] = S[j][i];
    kernel(T,S);
    if(T.NumRows()==0) return -1;
    for(i=0; i<length(); i++) {
        clear(z);
        for(j=0; j<t; j++) {
            mul(w, H[i][j], T[0][j]);
            z += w;
        }
        if(IsZero(z)) l.append(i);
    }
    k = l.length();
    S.SetDims(k+1, H.NumCols());
    for(i=0; i<k; i++) S[i] = H[l[i]];
    S[k] = s;
    kernel(T,S);
    if(T.NumRows()!=1) return -1;
    inv(z, T[0][k]);
    for(i=j=0; i<k; i++) {
        if(IsZero(T[0][i])) continue;
        j++;
        if(l[i] >= dim()) continue;
        mul(w, z, T[0][i]);
        c[l[i]] += w;
    }
    return j;
}

long HCCode::decode(GF2X& s, const Vec<GF2>& t)
// s = decoding of t (divided into data of length n)
// length of s is multiple of k
// returns the number of corrected chars
// negative value is returned if uncorrected chars exist
{
    long i,j,k,l,m(GF2E::degree());
    GF2X f,g;
    Vec<GF2> u,v;
    Vec<GF2E> c,d;
    d.SetLength(length());
    conv(f,t);
    for(i=k=l=0; !IsZero(f);) {
        trunc(g,f,m);
        conv(d[i++], g);
        f >>= m;
        if(!IsZero(f) && i < d.length()) continue;
        while(i < d.length()) clear(d[i++]);
        j = decode(c,d);
        if(j<0) l=1; else k+=j;
        for(i=j=0; j < c.length(); j++) {
            conv(g,c[j]);
            VectorCopy(u,g,m);
            v.append(u);
        }
    }
    conv(s,v);
    return (l ? -k:k);
}

long HCCode::decode(std::string& s, const Vec<GF2>& t) {
    long k;
    GF2X f;
    k = decode(f,t);
    s.resize(NumBytes(f));
    BytesFromGF2X((unsigned char*)s.data(), f, s.length());
    return k;
}