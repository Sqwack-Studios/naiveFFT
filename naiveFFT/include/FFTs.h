#pragma once
#include <cstdint>
#include <complex>

//TODO: Stockham

void CooleyTukey_outofplace(std::uint16_t N, std::uint16_t stride, std::complex<double>* input, std::complex<double>* output);
void CooleyTukey_inplace(std::uint16_t N, std::uint16_t stride, std::complex<double>* data);

void GentlemanSande_inplace(std::uint16_t N, std::uint16_t offset, std::complex<double>* data);
void GentlemanSande_outofplace(uint16_t N);


typedef struct Comp {
    /* comp of the form: a + bi */
    double a, b;
} Comp;
Comp comp_create(double a, double b);
Comp comp_add(Comp c1, Comp c2);
Comp comp_sub(Comp c1, Comp c2);
Comp comp_mul(Comp c1, Comp c2);


void comp_print(Comp comp);


/* Calculate e^(ix) */
Comp comp_euler(double x);
void dft(const Comp* sig, Comp* f, int n, int inv);
void fft(const Comp* sig, Comp* f, int s, int n, int inv);

void print_result(const Comp* sig, const Comp* sig0, int n);

void test_dft(const Comp* sig, Comp* f, Comp* sig0, int n);
void test_fft(const Comp* sig, Comp* f, Comp* sig0, int n);