// uses NTL
//   http://www.shoup.net/ntl

#include "RSCode.h"
#include "gaussj.cpp"

#define RS_GF2E_DEGMAX 8

RSCode::RSCode(long n, long k)
// n = data length after encoding
// k = data length before encoding
// assume GF2E has been initialized
{
    long i,j,m(GF2E::degree());
    GF2EX g,h,h1;
    GF2X f;
    if(m > RS_GF2E_DEGMAX) Error("GF2E degree too large");
    if(NumBits(n-1) >= m) Error("n is too large");
    H.SetDims(n,n-k);
    e.SetLength(1<<m);
    SetCoeff(g, H.NumCols());
    SetX(h);
    set(f);
    for(i=0; i<n; i++) {
        conv(h[0],f);
        InvMod(h1,h,g);
        VectorCopy(H[i], h1, H.NumCols());
        conv(f,h[0]);
        e[LongFromGF2X(f)] = i;
        f <<= 1;
    }
    kernel(G,H);
    gaussj(G);
}

void RSCode::encode(Vec<GF2E>& c, const Vec<GF2E>& b)
// c = encoding of b
{
    long i,j;
    GF2E t;
    VectorCopy(c, b, length());
    for(j=dim(); j<length(); j++) {
        for(i=0; i<dim(); i++) {
            mul(t, b[i], G[i][j]);
            c[j] += t;
        }
    }
}

void RSCode::encode(Vec<GF2>& t, const GF2X& s)
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

void RSCode::encode(Vec<GF2>& t, const std::string& s) {
    GF2X f;
    GF2XFromBytes(f, (const unsigned char*)s.c_str(), s.length());
    encode(t,f);
}

long RSCode::decode(Vec<GF2E>& c, const Vec<GF2E>& d)
// c = decoding of d by euclidean algorithm
// returns the number of corrected characters
// returns -1 if unable to correct
{
    long i,j,k(H.NumCols()>>1);
    GF2EX s,t,v(1),u,q,r;
    GF2E z,w;
    Vec<Pair<GF2EX, long> > f;
    VectorCopy(c, d, dim());
    mul(s.rep, d, H);
    s.normalize();
    if(IsZero(s)) return 0;
    SetCoeff(t, H.NumCols());
    while(deg(s) >= k) {// euclid
        DivRem(q,r,t,s);
        mul(t,q,v);
        sub(t,u,t);
        u = v;
        v = t;
        t = s;
        s = r;
    }
    if(IsZero(s) || IsOne(v)) return -1;
    diff(t,v);
    MakeMonic(v);
    CanZass(f,v);
    if(f.length() < deg(v)) return -1;
    for(i=0; i<f.length(); i++) {
        j = e[LongFromGF2X(rep(ConstTerm(f[i].a)))];
        if(j>=dim()) continue;
        eval(z, s, ConstTerm(f[i].a));
        eval(w, t, ConstTerm(f[i].a));
        z /= w;
        c[j] -= z;
    }
    return f.length();
}

long RSCode::decode(GF2X& s, const Vec<GF2>& t)
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

long RSCode::decode(std::string& s, const Vec<GF2>& t) {
    long k;
    GF2X f;
    k = decode(f,t);
    s.resize(NumBytes(f));
    BytesFromGF2X((unsigned char*)s.data(), f, s.length());
    return k;
}