#pragma once
#include <complex>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include "FFTs.h"

void buildData(std::uint32_t N, std::complex<double>* data)
{
    for (uint32_t i = 0; i < N; ++i)
    {
        data[i] = i;
    }
}

typedef std::complex<double> complex_t;

void bit_reverse(int N, complex_t* x) // bit reversal sorting
// N : sequence length
// x : input/output sequence
{
    int n_half = N >> 1;

    for (int i = 0, j = 1; j < N - 1; j++)
    {
        for (int k = n_half; k > (i ^= k); k >>= 1)
        {
        }


        if (i < j) std::swap(x[i], x[j]); // swap x[i] and x[j]
    }
}


typedef std::complex<double> complex_t;

void fft0(int n, int s, bool eo, complex_t* x, complex_t* y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    const int m = n / 2;
    const double theta0 = 2 * _PI / n;

    if (n == 2) {
        complex_t* z = eo ? y : x;
        for (int q = 0; q < s; q++) {
            const complex_t a = x[q + 0];
            const complex_t b = x[q + s];
            z[q + 0] = a + b;
            z[q + s] = a - b;
        }
    }
    else if (n >= 4) {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p * theta0), -sin(p * theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s * (p + 0)];
                const complex_t b = x[q + s * (p + m)];
                y[q + s * (2 * p + 0)] = a + b;
                y[q + s * (2 * p + 1)] = (a - b) * wp;
            }
        }
        fft0(n / 2, 2 * s, !eo, y, x);
    }
}

void fft(int n, complex_t* x) // Fourier transform
// n : sequence length
// x : input/output sequence
{
    complex_t* y = new complex_t[n];
    fft0(n, 1, 0, x, y);
    delete[] y;
    //for (int k = 0; k < n; k++) x[k] /= n;
}

int main()
{
    std::uint32_t N{ static_cast<uint32_t>(pow(2.0, 15.0))};
    std::vector<std::complex<double>> x;
    std::vector<std::complex<double>> y;
    std::vector<std::complex<double>> d;
    std::vector<std::complex<double>> g;

    std::vector<std::complex<double>> z;
    std::vector<std::complex<double>> w;
    
    x.resize(static_cast<size_t>(N));
    y.resize(static_cast<size_t>(N));
    d.resize(static_cast<size_t>(N));
    g.resize(static_cast<size_t>(N));
    z.resize(static_cast<size_t>(N));
    w.resize(static_cast<size_t>(N));
    
    buildData(N, x.data());
    buildData(N, y.data());
    buildData(N, d.data());
    buildData(N, g.data());

    std::chrono::steady_clock::time_point start;
    std::chrono::steady_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    ///////////COOLEY TUKEY/////////
    start = std::chrono::steady_clock::now();
    bit_reverse(N, x.data());
    cooleyTukey_inplace(N, 0, x.data());

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "(DIT) Cooley-Tukey in-place with previous bit-reverse permutation    " << elapsed_seconds.count() << "seconds" << '\n';

    /////////GENTLEMAN SANDE/////////
    start = std::chrono::steady_clock::now();
    gentlemanSande_inplace(N, 0, y.data());
    bit_reverse(N, y.data());

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "(DIF) Gentleman-Sande in-place with post bit-reverse permutation    " << elapsed_seconds.count() << "seconds" << '\n';


    start = std::chrono::steady_clock::now();
    DIF_slow_recursive_Stockham0(N, 1, d.data(), z.data());

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "(DIF) Stockham (self-sorting) slow      " << elapsed_seconds.count() << "seconds" << '\n';


    start = std::chrono::steady_clock::now();
    DIF_Stockham(N, 1, 0, g.data(), w.data());
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "(DIF) Stockham optimized memory access pattern:     " << elapsed_seconds.count() << "seconds" << '\n';

    return 0;
}