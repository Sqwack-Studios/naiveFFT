#pragma once
#include <cstdint>
#include <complex>

//TODO: Stockham

void cooleyTukey_outofplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);
void cooleyTukey_inplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* data);

void gentlemanSande_inplace(std::uint32_t N, std::uint32_t offset, std::complex<double>* data);
void gentlemanSande_outofplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);


void naiveDIFStockham0(std::uint32_t N, std::uint32_t stride, std::uint32_t eo, std::complex<double>* inout, std::complex<double>* work);
void naiveDIFStockham1(std::uint32_t N, std::uint32_t stride, std::uint32_t eo, std::complex<double>* input, std::complex<double>* output);

void DIF_recursive_Stockham0(std::uint32_t N, std::uint32_t stride, std::complex<double>* inout, std::complex<double>* work);
void DIF_recursive_Stockham1(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output);