#pragma once
#include "FFTs.h"
#include <vector>
#include <iostream>
#include <math.h>

constexpr double _PI  = 3.1415926535897931;
constexpr double _2PI = 6.2831853071795864;
//Naive FFT implementations

//During FFT demonstration, indices layout an in-place algorithm, suffering from time/decimation frequency that reverse the natural order(bit reversal permutations)
//An out-of-place implementation can be donde manipulating indices in order to retrieve natural order. This is why in-place CT-FFT has "wrong" access indices to the output array compared to theoretical demonstration

//To do it inplace, you could split the series just in half, but it requires bit permutation. You can do it with an input in natural order, but it requires to play with stride and offsets (doing it out-of-place).
//In fact, if you halve the input stream in even/odd series recursively, when the recursion ends, the resulting series are bit reversed, so, obviously, output is in natural order.

//Consider a 4 - point FFT
//                                                                                                 and now we start computing the FFT backwards. The "input" is reversed.
//x[0] -> even   | now we have two DFT  even  { x[0], x[2] } = {x_e[0], x_e[1]}     | and another four series   even-even {x[0]} = {x_ee'[0]}           |                              X[0]
//x[1]           |                                                                  |                           even-odd  {x[2]} = {x_eo'[0]}           |                              X[1]
                                                                                                                                                                                       
//x[2] -> even   |                      odd   { x[1], x[3] } = {x_o[0], x_o[1]}     |                           odd-even  {x[1]} = {x_oe'[0]}           |                              X[2]
//x[3]           |                                                                  |                           odd-odd	  {x[3]} = {x_oo'[0]}           |                              X[3]


//As I said, out-of-place algorithms consume O(N) memory, but don't need bit reversal.
//In-place algorithms do not consume extra memory, but need tricks to cope with the butterfly bit reversal

void CooleyTukey_outofplace(uint16_t N, uint16_t stride, std::complex<double>* input, std::complex<double>* output)
{

	if (N == 1)
	{
		*output = *input;
		return;
	}

	const double theta = -_2PI / static_cast<double>(N);
	std::complex<double> omega = { cos(theta), sin(theta) };

	const uint16_t half_N{ N / 2u };

	uint16_t new_s = stride << 1;
	CooleyTukey_outofplace(half_N, new_s, input, output); // even
	CooleyTukey_outofplace(half_N, new_s, input + stride, output + half_N); //odd

	for (uint16_t i = 0; i < half_N; ++i)
	{
		std::complex<double> aux = output[i];
		std::complex<double> wi = pow(omega, static_cast<float>(i));
		output[i] = aux + wi * output[i + half_N];
		output[i + half_N] = aux - wi * output[i + half_N];
	}
		
}

void CooleyTukey_inplace(uint16_t N, uint16_t stride, std::complex<double>* data)
{
    if (N <= 1)
        return;

    const uint16_t half_N{ N / 2u };
    const double theta{ -_2PI / static_cast<double>(N) };

    CooleyTukey_inplace(half_N, stride, data);//even
    CooleyTukey_inplace(half_N, stride + half_N, data);//odd

    std::complex<double> omega = { cos(theta), sin(theta) };

    for (uint16_t i = 0; i < half_N; ++i)
    {
        const std::complex<double> wi = pow(omega, static_cast<float>(i));
        const std::complex<double> Y_i{ data[i + stride]};
        const std::complex<double> Z_i{ wi * data[i + stride + half_N]};

        data[stride + i] = Y_i + Z_i;
        data[stride + i + half_N] = Y_i - Z_i;
    }
}

void GentlemanSande_outofplace(uint16_t N)
{

}
void GentlemanSande_inplace(uint16_t N)
{

}

Comp comp_create(double a, double b)
{
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
    puts("## Cooley�Tukey FFT ##");
    fft(sig, f, 1, n, 0);
    for (i = 0; i < n; i++)
        comp_print(f[i]);
    puts("----------------------");
    fft(f, sig0, 1, n, 1);
    print_result(sig, sig0, n);
    puts("######################");
}