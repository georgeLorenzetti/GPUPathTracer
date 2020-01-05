#pragma once

#ifndef __CUDACC__
#define __launch_bounds__(x, y)


template <typename T, int TT>
class surface{};
template <typename T, int TT>
class texture{};
#endif

template <typename T>
__device__ __host__ inline static void swap(T& a, T& b)
{
	T c(a); a = b; b = c;
}