#pragma once
#include <cstdint>
#include <complex>
constexpr double _PI = 3.1415926535897931;
constexpr double _2PI = 6.2831853071795864;
//TODO: Stockham

void cooleyTukey_outofplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);
void cooleyTukey_inplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* data);

void gentlemanSande_inplace(std::uint32_t N, std::uint32_t offset, std::complex<double>* data);
void gentlemanSande_outofplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);


void naiveDIFStockham0(std::uint32_t N, std::uint32_t stride, std::uint32_t eo, std::complex<double>* inout, std::complex<double>* work);
void naiveDIFStockham1(std::uint32_t N, std::uint32_t stride, std::uint32_t eo, std::complex<double>* input, std::complex<double>* output);

void DIF_slow_recursive_Stockham0(std::uint32_t N, std::uint32_t stride, std::complex<double>* inout, std::complex<double>* work);
void DIF_slow_recursive_Stockham1(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);

void DIF_fast_recursive_Stockham0(std::uint32_t N, std::uint32_t stride, std::complex<double>* inout, std::complex<double>* support);
void DIF_fast_recursive_Stockham1(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);


//slot: [0, 1]
//During the first call, X is our input array, and Y is our work area. However, this algorithm will consume the input data, and final result will be written in X.
//bool slot will keep track of which pointer slot is our actual output data pointer. If slot = 0, x is pointing to our inout array, if slot = 1, y is pointing to our inout array.
void DIF_Stockham(std::uint32_t N, std::uint32_t stride, bool slot, std::complex<double>* x, std::complex<double>* y);