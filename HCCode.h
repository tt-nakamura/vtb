// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/GF2XFactoring.h>
#include<NTL/GF2EXFactoring.h>
#include<NTL/mat_GF2E.h>
using namespace NTL;

struct HCCode {// Hermite curve code
    Vec<long> a,b;// exponents of (x,y) in H
    Mat<GF2E> G;// generator matrix
    Mat<GF2E> H;// check matrix
    Mat<long> K;// indices for syndromes
    long q;// sqrt(order of finite field)
    long u,t,n_err_max;
    HCCode(long, long);
    void encode(Vec<GF2E>&, const Vec<GF2E>&);
    void encode(Vec<GF2>&, const GF2X&);
    void encode(Vec<GF2>&, const std::string&);
    long decode(Vec<GF2E>&, const Vec<GF2E>&);
    long decode(GF2X&, const Vec<GF2>&);
    long decode(std::string&, const Vec<GF2>&);
    inline long dim() { return G.NumRows(); }// k
    inline long length() { return H.NumRows(); }// n
};

void GF2XFromLong(GF2X&, long);
