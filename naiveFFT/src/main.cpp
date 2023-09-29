#pragma once
#include <complex>
#include <cmath>
#include <vector>
#include <iostream>
typedef std::complex<double> complex_t;

#include "FFTs.h"

void F(int N, int q, complex_t* x)
// N : sequence length
// q : block start point (initial value is 0)
// x : input/output sequence
{
    const int m = N / 2;
    const double theta0 = 2 * 3.141516 / N;

    if (N > 1) {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p * theta0), -sin(p * theta0));
            const complex_t a = x[q + p + 0];
            const complex_t b = x[q + p + m];
            x[q + p + 0] = a + b;
            x[q + p + m] = (a - b) * wp;
        }
        F(N / 2, q + 0, x); // even position components
        F(N / 2, q + m, x); // odd position components
    }
}

void bit_reverse(int N, complex_t* x) // bit reversal sorting
// N : sequence length
// x : input/output sequence
{
    int n_half = N >> 1;

    for (int i = 0, j = 1; j < N - 1; j++) 
    {
        for (int k = n_half; k > (i ^= k); k >>= 1)
        { }


        if (i < j) std::swap(x[i], x[j]); // swap x[i] and x[j]
    }
}






int main()
{
    std::vector<std::complex<double>> x{ {5.0, 0.0}, {3.0, 0.0}, {2.0, 0.0}, {1.0, 0.0}/*, {10.0, 0.0}, {7.0, 0.0}, {9.0, 0.0}, {4.0, 0.0}*/ };
    int N = 4;
    int q = 0;
    //F(N, q, x.data());
    //bit_reverse(N, x.data());
    //cono(N);

    int n, i, k;
    std::vector<Comp> sig{ {5.0, 0.0}, {3.0, 0.0}, {2.0, 0.0}, {1.0, 0.0}/*, {10.0, 0.0}, {7.0, 0.0}, {9.0, 0.0}, {4.0, 0.0}*/ };
    std::vector<Comp> f, sig0;

    n = sig.size();
    f.resize(n);
    sig0.resize(n);


    //std::vector<std::complex<double>> x{ {5.0, 0.0}, {3.0, 0.0}, {2.0, 0.0}, {1.0, 0.0}, {10.0, 0.0}, {7.0, 0.0}, {9.0, 0.0}, {4.0, 0.0} };
    std::vector<std::complex<double>> y;
    y.resize(n);
    test_dft(sig.data(), f.data(), sig0.data(), n);
    test_fft(sig.data(), f.data(), sig0.data(), n);

    CooleyTukey_outofplace(n, 1, x.data(), y.data());
    bit_reverse(n, x.data());
    CooleyTukey_inplace(n, 0, x.data());
    return 0;
}