#include "ConvolCode.h"
#include "RSCode.h"
#include "HCCode.h"
#include<fstream>

main() {
    long i,j,k,n,L(1<<12),M(100),N(50);
    double e,v,x,x1(0),x2(0.25),dx((x2-x1)/N);
    GF2X f,g;
    Vec<GF2> c;

    ConvolCode C;
    std::ofstream f1("fig4a.txt");
    for(i=0; i<=N; i++) {
        x = x1 + i*dx;
        n = long(L*x);
        e = v = 0;
        for(j=0; j<M; j++) {
            random(f,L);
            C.encode(c,f);
            for(k=0; k<n; k++) c[RandomBnd(c.length())]++;
            C.decode(g,c);
            g -= f;
            k = weight(g);
            e += k;
            v += k*k;
        }
        e /= M;
        v /= M;
        v -= e*e;
        v = sqrt(v);
        e /= L;
        v /= L;
        f1 << x << '\t' << e << '\t' << v << '\n';
        std::cout << x << '\t' << e << '\t' << v << '\n';
    }

    BuildSparseIrred(f,8);
    GF2E::init(f);
    RSCode R(32,16);
    std::ofstream f2("fig4b.txt");
    for(i=0; i<=N; i++) {
        x = x1 + dx*i;
        n = long(L*x);
        e = v = 0;
        for(j=0; j<M; j++) {
            random(f,L);
            R.encode(c,f);
            for(k=0; k<n; k++) c[RandomBnd(c.length())]++;
            R.decode(g,c);
            g -= f;
            k = weight(g);
            e += k;
            v += k*k;
        }
        e /= M;
        v /= M;
        v -= e*e;
        v = sqrt(v);
        e /= L;
        v /= L;
        f2 << x << '\t' << e << '\t' << v << '\n';
        std::cout << x << '\t' << e << '\t' << v << '\n';
    }

    BuildSparseIrred(f,4);
    GF2E::init(f);
    HCCode H(64,32);
    std::ofstream f3("fig4c.txt");
    for(i=0; i<=N; i++) {
        x = x1 + dx*i;
        n = long(L*x);
        e = v = 0;
        for(j=0; j<M; j++) {
            random(f,L);
            H.encode(c,f);
            for(k=0; k<n; k++) c[RandomBnd(c.length())]++;
            H.decode(g,c);
            g -= f;
            k = weight(g);
            e += k;
            v += k*k;
        }
        e /= M;
        v /= M;
        v -= e*e;
        v = sqrt(v);
        e /= L;
        v /= L;
        f3 << x << '\t' << e << '\t' << v << '\n';
        std::cout << x << '\t' << e << '\t' << v << '\n';
    }
}