// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/GF2X.h>
using namespace NTL;

void GF2XFromLong(GF2X& x, long a) {
    clear(x);
    x.SetLength(NumBits(a));
    for(long i=0; a; a>>=1, i++) if(a&1) set(x[i]);
}

long LongFromGF2X(const GF2X& a) {
    if(IsZero(a)) return 0;
    long i,x(1);
    for(i=deg(a)-1; i>=0; i--) {
        x<<=1;
        if(IsOne(a[i])) x|=1;
    }
    return x;
}

void GF2XFromString(GF2X& x, const std::string& a) {
    GF2XFromBytes(x, (const unsigned char*)a.c_str(), a.length());
}

void StringFromGF2X(std::string& x, const GF2X& a) {
    x.resize(NumBytes(a));
    BytesFromGF2X((unsigned char*)x.data(), a, x.length());
}