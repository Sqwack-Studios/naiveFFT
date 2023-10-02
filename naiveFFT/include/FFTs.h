#pragma once
#include <cstdint>
#include <complex>

//TODO: Stockham

void cooleyTukey_outofplace(std::uint16_t N, std::uint16_t stride, std::complex<double>* input, std::complex<double>* output);
void cooleyTukey_inplace(std::uint16_t N, std::uint16_t stride, std::complex<double>* data);

void gentlemanSande_inplace(std::uint16_t N, std::uint16_t offset, std::complex<double>* data);
void gentlemanSande_outofplace(std::uint16_t N, std::uint16_t stride, std::complex<double>* input, std::complex<double>* output);


