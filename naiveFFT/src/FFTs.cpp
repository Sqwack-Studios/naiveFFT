#pragma once
#include "FFTs.h"
#include <vector>
#include <iostream>
#include <math.h>

constexpr double PI = 3.141592653589793;
//Naive FFT implementations

//During FFT demonstration, indices layout an in-place algorithm, suffering from time/decimation frequency that reverse the natural order(bit reversal permutations)
//An out-of-place implementation can be donde manipulating indices in order to retrieve natural order. This is why in-place CT-FFT has "wrong" access indices to the output array compared to theoretical demonstration

//To do it inplace, you could split the series just in half, but it requires bit permutation. You can do it with an input in natural order, but it requires to play with stride and offsets (doing it out-of-place).
//In fact, if you halve the input stream in even/odd series recursively, when the recursion ends, the resulting series are bit reversed, so, obviously, output is in natural order.

//Consider a 4 - point FFT
//																											and now we start computing the FFT backwards. The "input" is reversed.
//x[0] -> even   | now we have two DFT  even  { x[0], x[2] } = {x_e[0], x_e[1]}  	| and another four series   even-even {x[0]} = {x_ee'[0]}   		|																X[0]
//x[1]           |                                         							|							even-odd  {x[2]} = {x_eo'[0]}			|																X[1]

//x[2] -> even	 |						odd   { x[1], x[3] } = {x_o[0], x_o[1]}		|							odd-even  {x[1]} = {x_oe'[0]}			|																X[2]
//x[3]			 |																	|                           odd-odd	  {x[3]} = {x_oo'[0]}			|																X[3]


//As I said, out-of-place algorithms consume O(N) memory, but don't need bit reversal.
//In-place algorithms do not consume extra memory, but need tricks to cope with the butterfly bit reversal

void CooleyTukey_outofplace(uint16_t N, uint16_t stride, std::complex<double>* input, std::complex<double>* output)
{

	if (N == 1)
	{
		*output = *input;
		return;
	}

	const double theta = -2.0 * PI / N;
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

}

void GentlemanSande_outofplace(uint16_t N)
{

}
void GentlemanSande_inplace(uint16_t N)
{

}