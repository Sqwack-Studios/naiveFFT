#pragma once
#include "FFTs.h"
#include <vector>
#include <iostream>
#include <math.h>


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



//So, the thing is that using Gentleman-Sade formulation, we have that for an input sequence x[p], an output sequence X_k is generated so that:
//X_k = X_2k + X_2k+1, where:
//X_2k = F_N/2 (x[p] + x[p + n/2])     are the half-sequence lenght fourier transforms
//X_2k+1 = F_N/2 ((x[p] - x[p + n/2]))Wp_N
//This means that we have a relationship between the input/output sequences that we can exploit if we want to build an out-of-place algorithm.
//Moreover, we can use that to create an algorithm that is self-sorting
// 
//For a N length sequence where N = 2^L, a 2-radix implementation will have up to L steps.
//Each step we have an input array x_h, and we output into array X that will get self-sorted based in previous relationship.
//The output array for each step will be input for the next step, so:
//x_h+1(q, p) = x_h(q,p) + x_h(q, p + n/2)
//x_h+1(q, p) = ( x_h(q,p) - x_h(q, p + n/2) )Wp_N
//
//Where, foreach H step, q = 0, 1,..., stride - 1
//                     p = 0, 1,..., n/2 - 1
//Where p is the butterfly pass index, and q represents an offset for even/odd splits of the sequence.
//
//Let's see what happens for an 8-point FFT: we have two memory arrays of size N; v0 and v1 
//     ______________________________________________________________
//    |      ||      ||      ||      ||      ||      ||      ||      |
//v0  | v0[0]|| v0[1]|| v0[2]|| v0[3]|| v0[4]|| v0[5]|| v0[6]|| v0[7]| -> inout array
//    |______||______||______||______||______||______||______||______|
//     ______________________________________________________________
//    |      ||      ||      ||      ||      ||      ||      ||      |
//v1  | v1[0]|| v1[1]|| v1[2]|| v1[3]|| v1[4]|| v1[5]|| v1[6]|| v1[7]| -> support array
//    |______||______||______||______||______||______||______||______|
//
//   [h = 0, stage-0; stride = 1; N = 8; half_N = 4 ]
//fft0(8, 1, 0, v0, v1)
//we access whole array: v0[p]
//we write whole array: v1[p]
//
//   [h = 1, stage-1; stride = 2; N = 4; half_N = 2 ]
//fft1(4, 2, 0, v1, v0)                     fft1(4, 2, 1, v1, v0)
//
// for each substage, p = 0, 1; q depends of even/odd split
// for the even part, we compute butterfly components accessing v1 addresses:
// v1[0], v1[4]; v1[2], v2[6]  -> write into v0[0], v0[2]; v0[4], v0[6]; we only read/write into even addresses
// 
// odd part, q = 1:
// v1[1], v1[5]; v1[3], v2[7]  -> write into v0[1], v0[3]; v0[5], v0[7]; we only read/write into odd addresses
//
//   [h = 2, stage-2; stride = 4; N = 2; half_N = 1 ]
//fft0(2, 4, 0, v0, v1)    fft0(2, 4, 2, v0, v1)    fft0(4, 2, 1, v0, v1)    fft0(4, 2, 3, v0, v1)    
//
//There is, however, a problem as we will see now. The deeper we get into each stage, read/writes distance is large.
// 
// stage-2, even-even path  fft0(2, 4, 0, v0, v1) 
// p = 0, q = 0, we access v0 addresses:
// v0[0], v0[4]; v1[0], v1[4] ... etc
//
// stage-2, even-odd path  fft0(2, 4, 2, v0, v1) 
// p = 0, q = 2, we access v0 addresses:
// v0[2], v0[6]; v1[2], v1[6]
//We can make coalesced memory accesses if we make a single function that computes both, even/odd terms on the fly without making a function call for each specific even/odd split.
//If we look closely, for the stage-1
//even/odd access is interleaved. We could compute even->odd->even->odd->even->odd to help with memory localty.
//To do this, we loop over each q value, which is our even/odd index selector when we are computing our butterfly passes (during p looping)

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

void stockhamIterativeButterflyPass(std::uint32_t half_N, std::uint32_t stride, std::complex<double> w, std::complex<double>* in, std::complex<double>* out)
{
    for (std::uint32_t q = 0; q < stride; ++q)
    {
        for (std::uint32_t p = 0; p < half_N; ++p)
        {
            const std::complex<double> w_p{ pow(w, static_cast<float>(p)) };
            const std::complex<double> a{ in[q + stride * p] };
            const std::complex<double> b{ in[q + stride * (p + half_N)] };

            out[q + stride * (2u * p)] = a + b;
            out[q + stride * (2u * p + 1u)] = (a - b) * w_p;
        }
    }
}

void stockhamFastIterativeButterflyPass(std::uint32_t half_N, std::uint32_t stride, std::complex<double> w, std::complex<double>* in, std::complex<double>* out)
{
    for (std::uint32_t p = 0; p < half_N; ++p)
    {
           const std::complex<double> w_p {pow(w, static_cast<float>(p)) };
        for (std::uint32_t q = 0; q < stride; ++q)
        {
            const std::complex<double> a{ in[q + stride *  p          ] };
            const std::complex<double> b{ in[q + stride * (p + half_N)] };
            
            out[q + stride * (2u * p)]      = a + b;
            out[q + stride * (2u * p + 1u)] = (a - b) * w_p;
        }
    }
}


//Naive implementation: memory access pattern is horrible
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


//This implementation is aware of memory access patterns, but we made a critical mistake, as I will explain later
void DIF_slow_recursive_Stockham0(std::uint32_t N, std::uint32_t stride, std::complex<double>* inout, std::complex<double>* work)
{
    if (N <= 1u)
        return;

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    stockhamIterativeButterflyPass(half_N, stride, omega, inout, work);

    DIF_slow_recursive_Stockham1(half_N, stride << 1u, work, inout);
}

void DIF_slow_recursive_Stockham1(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output)
{
    if (N <= 1)
    {
        for (std::uint32_t q = 0; q < stride; ++q)
        {
            output[q] = input[q];
        }
        return;
    }

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    stockhamIterativeButterflyPass(half_N, stride, omega, input, output);

    DIF_slow_recursive_Stockham0(half_N, stride << 1u, output, input);
}


//We can do better. If we recall what we saw later, the deeper we get in the FFT stage, larger is the distance between read/write indices. Cache misses -> slow program.
//We tried to solve this by removing independent even/odd recursive calls, grouping them in the same call. This means that we process each even/odd split in the same stage.
//Let's see a 16-point radix-2 FFT, lets see what happens at stage = 2:
//   [h = 2, stage-2; stride = 4; N = 4; half_N = 2 ]
//fft0(2, 4, 0, v0, v1)    fft0(2, 4, 2, v0, v1)    fft0(4, 2, 1, v0, v1)    fft0(4, 2, 3, v0, v1)    
//
//We have up to 4 calls, 2 even and 2 odd. Now our stride is 4, so distance between read/writes is, at least, 4. And, as we can see, q index moves from [0, 3], that is: q = 0, ..., stride - 1
//So instead of making 4 calls, we just make a single call, computing each split foreach possible value of q. That is what is going on in DIF_slow_recursive implementation
//But it's slow.
//
//Following last example, let's see the memory access pattern for stage-2. Looping foreach p in each q value. q = 0, 1, 2, 3; p = 0, 1
//q = 0; p = 0
//     Access: v0[0], v0[8];  write: v1[0], v1[4]
//q = 0; p = 1
//     Access: v0[4], v0[12]; write: v1[8], v1[12] 
//q = 1; p = 0
//     Access: v0[1], v0[9];  write: v1[1], v1[5] 
//q = 1; p = 1
//     Access: v0[5], v0[13]; write: v1[9], v1[13] 
//q = 2; p = 0
//     Access: v0[2], v0[10]; write  v1[2], v1[6]
//q = 2; p = 1
//     Access: v0[6], v0[14]; write: v1[10], v1[14] 
//q = 3; p = 0
//     Access: v0[3], v0[11]; write: v1[3], v1[7] 
//q = 3; p = 1
//     Access: v0[7], v0[15]; write: v1[11], v1[15] 
//
//As we can see, we jump through memory a lot. Read/Write access distance is the stride, and the problem is increased the deeper, and larger the data set gets.
//But we can see a pattern. If we take a look at p = 0 items, we see that access for consecutive values of q has stride = 1!
//
//q = 0; p = 0
//     Access: v0[0], v0[8];  write: v1[0], v1[4]
//q = 1; p = 0
//     Access: v0[1], v0[9];  write: v1[1], v1[5] 
//q = 2; p = 0
//     Access: v0[2], v0[10]; write  v1[2], v1[6]
//q = 3; p = 0
//     Access: v0[3], v0[11]; write: v1[3], v1[7] 
/////////////////////
//q = 0; p = 1
//     Access: v0[4], v0[12]; write: v1[8], v1[12] 
//q = 1; p = 1
//     Access: v0[5], v0[13]; write: v1[9], v1[13] 
//q = 2; p = 1
//     Access: v0[6], v0[14]; write: v1[10], v1[14] 
//q = 3; p = 1
//     Access: v0[7], v0[15]; write: v1[11], v1[15] 
/////////////////////
//By swapping the loops, we have a better memory access pattern.
//


void DIF_fast_recursive_Stockham0(std::uint32_t N, std::uint32_t stride, std::complex<double>* inout, std::complex<double>* support)
{
    if (N <= 1u)
        return;

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };
    stockhamFastIterativeButterflyPass(half_N, stride, omega, inout, support);

    DIF_fast_recursive_Stockham1(half_N, stride << 1u, support, inout);
}

void DIF_fast_recursive_Stockham1(std::uint32_t N, std::uint32_t stride, std::complex<double>* input, std::complex<double>* output)
{
    if (N <= 1u)
    {
        for (std::uint32_t q = 0; q < stride; ++q)
        {
            output[q] = input[q];
        }
        return;
    }

    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    const std::complex<double> omega{ cos(theta0), sin(theta0) };

    stockhamFastIterativeButterflyPass(half_N, stride, omega, input, output);

    DIF_fast_recursive_Stockham0(half_N, stride << 1u, output, input);
}

//And finally, there's still something we can improve. Because we are swapping input/output in each call, based on the length of the sequence, we might have to make extra copies from X->Y, because Y will be our output array.
//We can have a boolean that tracks which pointer is the actual output, and once we get to the last stage of the FFT, we can simplify the butterfly pass and select the correct output array, avoiding extra memory accesses.

void DIF_Stockham(std::uint32_t N, std::uint32_t stride, bool slot, std::complex<double>* x, std::complex<double>* y)
{
    const std::uint32_t half_N{ N / 2u };
    const double theta0{ -_2PI / static_cast<double>(N) };
    if (N == 2)
    {
        //Simplified butterfly step. If N = 2, omega = -1, so our twiddle factor is always 1. 
        // We avoid a multiplication per butterfly, we avoid computing sincos and a pow() call.

       std::complex<double> * const Z{ slot ? y : x };
       for (std::uint32_t q = 0; q < stride; ++q)
       {
           const std::complex<double> a{ x[q] };
           const std::complex<double> b{ x[q + stride] };

           Z[q] = a + b;
           Z[q + stride] = a - b;

       }
        return;
    }
    else
    {
        const std::complex<double> omega{ cos(theta0), sin(theta0) };

        stockhamFastIterativeButterflyPass(half_N, stride, omega, x, y);

        DIF_Stockham(half_N, stride << 1u, !slot, y, x);
    }

}
