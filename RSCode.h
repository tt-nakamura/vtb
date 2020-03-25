// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/GF2XFactoring.h>
#include<NTL/GF2EXFactoring.h>
#include<NTL/mat_GF2E.h>
using namespace NTL;

struct RSCode {// Reed-Solomon code
    Mat<GF2E> G;// generator matrix
    Mat<GF2E> H;// check matrix
    Vec<long> e;// table of exponents
    RSCode(long, long);
    void encode(Vec<GF2E>&, const Vec<GF2E>&);
    void encode(Vec<GF2>&, const GF2X&);
    void encode(Vec<GF2>&, const std::string&);
    long decode(Vec<GF2E>&, const Vec<GF2E>&);
    long decode(GF2X&, const Vec<GF2>&);
    long decode(std::string&, const Vec<GF2>&);
    inline long length() { return H.NumRows(); }// n
    inline long dim() { return G.NumRows(); }// k
};

long LongFromGF2X(const GF2X&);