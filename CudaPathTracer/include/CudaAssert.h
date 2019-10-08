#pragma once
#include <precomp.h>

#define cudaAssert(...) cuda_assert((cuda##__VA_ARGS__), __FILE__, __LINE__, true);

inline cudaError cuda_assert(const cudaError code, const char* file, int line, bool abort){
	if (code != cudaSuccess){
		fprintf(stderr, "cuda_assert: %s %s %d\n", cudaGetErrorString(code), file, line);

		if (abort){
			cudaDeviceReset();
			exit(code);
		}
	}

	return code;
}