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

void cooleyTukey_outofplace(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output)
{

	if (N == 1)
	{
		*output = *input;
		return;
	}

	const double theta = -_2PI / static_cast<double>(N);
	std::complex<double> omega = { cos(theta), sin(theta) };

	const std::uint32_t half_N{ N / 2u };

    std::uint32_t new_s = stride << 1;
	cooleyTukey_outofplace(half_N, new_s, input, output); // even
	cooleyTukey_outofplace(half_N, new_s, input + stride, output + half_N); //odd

	for (std::uint32_t i = 0; i < half_N; ++i)
	{
		const std::complex<double> aux = output[i];
		const std::complex<double> wi = pow(omega, static_cast<float>(i));
        const std::complex<double> b = output[i + half_N] * wi;

		output[i] = aux + b;
		output[i + half_N] = aux - b;
	}
		
}

void cooleyTukey_inplace(std::uint32_t N, std::uint32_t offset, std::complex<double>* data)
{
    if (N <= 1)
        return;

    const std::uint32_t half_N{ N / 2u };
    const double theta{ -_2PI / static_cast<double>(N) };

    cooleyTukey_inplace(half_N, offset, data);//even
    cooleyTukey_inplace(half_N, offset + half_N, data);//odd

    std::complex<double> omega = { cos(theta), sin(theta) };

    for (std::uint32_t i = 0; i < half_N; ++i)
    {
        const std::complex<double> wi = pow(omega, static_cast<float>(i));
        const std::complex<double> Y_i{ data[i + offset]};
        const std::complex<double> Z_i{ wi * data[i + offset + half_N]};

        data[offset + i] = Y_i + Z_i;
        data[offset + i + half_N] = Y_i - Z_i;
    }
}


//This is super slow
void gentlemanSande_outofplace(std::uint32_t N, std::uint32_t q, std::complex<double>* x, std::complex<double>* y)
{
    if (N <= 1)
        return;

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    for (std::uint32_t i = 0u; i < half_N; ++i)
    {
        const std::uint32_t aI{ static_cast<std::uint32_t>(q + i) };
        const std::uint32_t bI{ static_cast<std::uint32_t>(q + i + half_N) };
        const std::complex<double> w_i = pow(omega, static_cast<float>(i));
        const std::complex<double> a = x[aI];
        const std::complex<double> b = x[bI];

        y[aI] = a + b;
        y[bI] = (a - b) * w_i;
    }

    gentlemanSande_outofplace(half_N, q, x, y);
    gentlemanSande_outofplace(half_N, q + half_N, x, y);

    for (std::uint32_t i = 0u; i < half_N; ++i)
    {
        x[q + 2u * i] = y[q + i];
        x[q + 2u * i + 1u] = y[q + i + half_N];
    }
}

void gentlemanSande_inplace(std::uint32_t N, std::uint32_t offset, std::complex<double>* data)
{
    if (N <= 1)
        return;

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };

 
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    for (std::uint32_t i = 0u; i < half_N; ++i)
    {
        const std::uint32_t target{ static_cast<std::uint32_t>(i + offset) };
        const std::uint32_t targetDisplaced{ static_cast<std::uint32_t>(target + half_N) };

        const std::complex<double> w_i = pow(omega, static_cast<float>(i));

        const std::complex<double> a{ data[target]};
        const std::complex<double> b{ data[targetDisplaced] };

        data[target] = a + b;
        data[targetDisplaced] = (a - b) * w_i;
    }

    gentlemanSande_inplace(half_N, offset, data);         //split1
    gentlemanSande_inplace(half_N, offset + half_N, data);//split2

}

void stockhamButterflyPass(std::uint32_t half_N, std::uint32_t stride, std::uint32_t eo, std::complex<double> w, std::complex<double>* in, std::complex<double>* out)
{
    for (std::uint32_t p = 0u; p < half_N; ++p)
    {
        const std::complex<double> w_p{ pow(w, static_cast<float>(p)) };
        const std::complex<double> a{ in[eo + stride * (p)] };
        const std::complex<double> b{ in[eo + stride * (p + half_N)] };

        out[eo + stride * (2u * p)] = a + b;
        out[eo + stride * (2u * p + 1u)] = (a - b) * w_p;
    }
}

void naiveDIFStockham0(std::uint32_t N, std::uint32_t stride, std::uint32_t eo, std::complex<double>* inout, std::complex<double>* work)
{
    if (N <= 1u)
        return;

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    stockhamButterflyPass(half_N, stride, eo, omega, inout, work);

    const std::uint32_t nStride{ stride << 1u };
    naiveDIFStockham1(half_N, nStride, eo         , work, inout);
    naiveDIFStockham1(half_N, nStride, eo + stride, work, inout);

}
void naiveDIFStockham1(std::uint32_t N, std::uint32_t stride, std::uint32_t eo, std::complex<double>* input, std::complex<double>* output)
{
    if (N <= 1u)
    {
        output[eo] = input[eo];
    }

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    stockhamButterflyPass(half_N, stride, eo, omega, input, output);

    const std::uint32_t nStride{ stride << 1u };

    naiveDIFStockham0(half_N, nStride, eo         , output, input);
    naiveDIFStockham0(half_N, nStride, eo + stride, output, input);
}