#pragma once
#include <stdint.h>
#include <complex>
void CooleyTukey_outofplace(uint16_t N, uint16_t stride, std::complex<double>* input, std::complex<double>* output);
void CooleyTukey_inplace(uint16_t N, uint16_t stride, std::complex<double>* data);

void GentlemanSande_outofplace(uint16_t N);
void GentlemanSande_inplace(uint16_t N);


typedef struct Comp {
    /* comp of the form: a + bi */
    double a, b;
} Comp;

Comp comp_create(double a, double b) {
    Comp res;
    res.a = a;
    res.b = b;
    return res;
}

Comp comp_add(Comp c1, Comp c2) {
    Comp res = c1;
    res.a += c2.a;
    res.b += c2.b;
    return res;
}

Comp comp_sub(Comp c1, Comp c2) {
    Comp res = c1;
    res.a -= c2.a;
    res.b -= c2.b;
    return res;
}

Comp comp_mul(Comp c1, Comp c2) {
    Comp res;
    res.a = c1.a * c2.a - c1.b * c2.b;
    res.b = c1.b * c2.a + c1.a * c2.b;
    return res;
}


void comp_print(Comp comp) {
    printf("%.6f + %.6f i\n", comp.a, comp.b);
}

/* const double PI = acos(-1); */
#define PI 3.141592653589793
#define SQR(x) ((x) * (x))

/* Calculate e^(ix) */
Comp comp_euler(double x) {
    Comp res;
    res.a = cos(x);
    res.b = sin(x);
    return res;
}

#define comp_mul_self(c, c2) \
do { \
    double ca = c->a; \
    c->a = ca * c2->a - c->b * c2->b; \
    c->b = c->b * c2->a + ca * c2->b; \
} while (0)

void dft(const Comp* sig, Comp* f, int n, int inv) {
    Comp ep = comp_euler(2 * (inv ? -PI : PI) / (double)n);
    Comp ei, ej, * pi = &ei, * pj = &ej, * pp = &ep;
    int i, j;
    pi->a = pj->a = 1;
    pi->b = pj->b = 0;
    for (i = 0; i < n; i++)
    {
        f[i].a = f[i].b = 0;
        for (j = 0; j < n; j++)
        {
            f[i] = comp_add(f[i], comp_mul(sig[j], *pj));
            comp_mul_self(pj, pi);
        }
        comp_mul_self(pi, pp);
    }
}

void fft(const Comp* sig, Comp* f, int s, int n, int inv) {
    int i, hn = n >> 1;
    Comp ep = comp_euler((inv ? PI : -PI) / (double)hn), ei;
    Comp* pi = &ei, * pp = &ep;
    if (!hn) *f = *sig;
    else
    {
        fft(sig, f, s << 1, hn, inv);
        fft(sig + s, f + hn, s << 1, hn, inv);
        pi->a = 1;
        pi->b = 0;
        for (i = 0; i < hn; i++)
        {
            Comp even = f[i], * pe = f + i, * po = pe + hn;
            comp_mul_self(po, pi);
            pe->a += po->a;
            pe->b += po->b;
            po->a = even.a - po->a;
            po->b = even.b - po->b;
            comp_mul_self(pi, pp);
        }
    }
}

void print_result(const Comp* sig, const Comp* sig0, int n) {
    int i;
    double err = 0;
    for (i = 0; i < n; i++)
    {
        Comp t = sig0[i];
        t.a /= n;
        t.b /= n;
        comp_print(t);
        t = comp_sub(t, sig[i]);
        err += t.a * t.a + t.b * t.b;
    }
    printf("Error Squared = %.6f\n", err);
}

void test_dft(const Comp* sig, Comp* f, Comp* sig0, int n) {
    int i;
    puts("## Direct DFT ##");
    dft(sig, f, n, 0);
    for (i = 0; i < n; i++)
        comp_print(f[i]);
    puts("----------------");
    dft(f, sig0, n, 1);
    print_result(sig, sig0, n);
    puts("################");
}

void test_fft(const Comp* sig, Comp* f, Comp* sig0, int n) {
    int i;
    puts("## Cooley–Tukey FFT ##");
    fft(sig, f, 1, n, 0);
    for (i = 0; i < n; i++)
        comp_print(f[i]);
    puts("----------------------");
    fft(f, sig0, 1, n, 1);
    print_result(sig, sig0, n);
    puts("######################");
}