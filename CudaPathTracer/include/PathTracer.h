#pragma once
#include <precomp.h>

//put the parameters in a struct so launch_kernels doesnt need to take a million parameters
class PathTracer{
public:
	PathTracer(int cores);
	void CalcImageParameters();
	void Trace(Renderer* cuda_interop, glm::vec4* frame_buffer);

private:
	KernelParams kernel_params;
	Camera camera;
	int frame;
};