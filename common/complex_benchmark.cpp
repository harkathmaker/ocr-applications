#include <complex>
#include <iostream>
#include <sys/time.h>
#include <unistd.h>

#define START_SEED 10
#define BENCH_LOOP 100000

class mycomplex {
public:
    mycomplex(double _re, double _img): re(_re), img(_img) {}

    inline mycomplex operator *(const mycomplex& other) {
        // real = (ac - bd)
        // imag = (bc + ad)
        return mycomplex(re*other.re - img*other.img, img*other.re + re*other.img);
    }

    inline mycomplex operator+(const mycomplex& other) {
        return mycomplex(re+other.re,img+other.img);
    }

    inline void operator+=(const mycomplex& other) {
        re += other.re;
        img += other.img;
    }

    inline double real() const {
        return re;
    }
    inline double imag() const {
        return img;
    }
private:
    double re;
    double img;
};

double mysecond() {
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void bench_mycomplex(unsigned int n) {
    double *d = new double[n*2];
    double *r = new double[n*2];

    mycomplex *c = (mycomplex*)d;
    mycomplex *rc = (mycomplex*)r;

    double t = mysecond();

    srand(START_SEED);
    for(unsigned int i=0;i<n;i++) {
        d[i] = 1 + rand() % 100;
        r[i] = 0;
    }

    printf("Setup time (mycomplex): %f\n",mysecond()-t);
    t = mysecond();

    for(unsigned int k=0;k<BENCH_LOOP;k++) {
        for(unsigned int i=0;i<n-1;i++) {
            rc[i] = c[i]*c[i+1];
        }
    }

    printf("Multiply time (mycomplex): %f\n",mysecond()-t);
    t = mysecond();

    for(unsigned int k=0;k<BENCH_LOOP;k++) {
        for(unsigned int i=0;i<n-1;i++) {
            rc[i] += c[i]+c[i+1];
        }
    }

    printf("Add time (mycomplex): %f\n",mysecond()-t);

    printf("Final result: d[n/2] = %f + %fi\n\n",rc[n/2].real(),rc[n/2].imag());

    delete d,r;
}

void bench_complex(unsigned int n) {
    double *d = new double[n*2];
    double *r = new double[n*2];
    
    std::complex<double> *c = (std::complex<double>*)d;
    std::complex<double> *rc = (std::complex<double>*)r;

    double t = mysecond();

    srand(START_SEED);
    for(unsigned int i=0;i<n*2;i++) {
        d[i] = 1 + rand() % 100;
        r[i] = 0;
    }

    printf("Setup time (complex): %f\n",mysecond()-t);
    t = mysecond();

    for(unsigned int k=0;k<BENCH_LOOP;k++) {
        for(unsigned int i=0;i<n-1;i++) {
            rc[i] = c[i]*c[i+1];
        }
    }

    printf("Multiply time (complex): %f\n", mysecond()-t);
    t = mysecond();


    for(unsigned int k=0;k<BENCH_LOOP;k++) {
        for(unsigned int i=0;i<n-1;i++) {
            rc[i] += c[i]+c[i+1];
        }
    }

    printf("Add time (complex): %f\n", mysecond()-t);

    printf("Final result: c[n/2] = %f + %fi\n\n",rc[n/2].real(),rc[n/2].imag());

    delete d,r;
}

void bench_doubles(unsigned int n) {
    double *d = new double[n*2];
    double *r = new double[n*2];

    double t = mysecond();

    srand(START_SEED);
    for(unsigned int i=0;i<n*2;i++) {
        d[i] = 1 + rand() % 100;
        r[i] = 0;
    }

    printf("Setup time (doubles): %f\n",mysecond()-t);
    t = mysecond();

    for(unsigned int k=0;k<BENCH_LOOP;k++) {
        for(unsigned int i=0;i<n*2-1;i+=2) {
            // real = (ac - bd)
            r[i] = (d[i]*d[i+2] - d[i+1]*d[i+3]);
            // imag = (bc + ad)
            r[i+1] = (d[i+1]*d[i+2] + d[i]*d[i+3]);
        }
    }

    printf("Multiply time (doubles): %f\n",mysecond()-t);
    t = mysecond();

    for(unsigned int k=0;k<BENCH_LOOP;k++) {
        for(unsigned int i=0;i<2*n-1;i+=2) {
            r[i] += d[i]+d[i+2];
            r[i+1] += d[i+1]+d[i+3];
        }
    }

    printf("Add time (doubles): %f\n",mysecond()-t);

    printf("Final result: d[n/2] = %f + %fi\n\n",r[(n/2)*2],r[(n/2)*2+1]);

    delete d,r;
}

int main() {
    double testComplex[2] = { 1.05, -2 };
    std::complex<double> *c = (std::complex<double>*)testComplex;
    printf("c = %f + %fi\n",c->real(),c->imag());
    printf("sizeof(complex<double>):%lu\n",sizeof(std::complex<double>));
    printf("sizeof(complex<float>):%lu\n",sizeof(std::complex<float>));
    printf("sizeof(mycomplex):%lu\n",sizeof(mycomplex));

    bench_complex(100000);
    bench_doubles(100000);
    bench_mycomplex(100000);

    return 0;
}
