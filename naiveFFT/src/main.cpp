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