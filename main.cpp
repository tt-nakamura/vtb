#include "ConvolCode.h"
#include "RSCode.h"
#include "HCCode.h"

main() {
    long i,N(20);
    std::string s,t;
    GF2X f;
    Vec<GF2> u;

    s = "Hello world of error correcting codes.";

    ConvolCode C;
    C.encode(u,s);
    for(i=0; i<N; i++) u[RandomBnd(u.length())]++;
    C.decode(t,u);
    std::cout << t << std::endl;

    BuildSparseIrred(f,8);
    GF2E::init(f);
    RSCode R(32,16);
    R.encode(u,s);
    for(i=0; i<N; i++) u[RandomBnd(u.length())]++;
    std::cout << R.decode(t,u) << std::endl;
    std::cout << t << std::endl;

    BuildSparseIrred(f,4);
    GF2E::init(f);
    HCCode H(64,32);
    H.encode(u,s);
    for(i=0; i<N; i++) u[RandomBnd(u.length())]++;
    std::cout << H.decode(t,u) << std::endl;
    std::cout << t << std::endl;
}
