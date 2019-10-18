#pragma once
#include <precomp.h>

//put the parameters in a struct so launch_kernels doesnt need to take a million parameters
class PathTracer{
public:
	PathTracer(int cores);
	void CalcImageParameters();
	void Trace(Renderer* cuda_interop, glm::vec4* frame_buffer);

	void TranslateCamera(glm::vec3 direction, float delta_time);
	void RotateCamera(glm::vec3 direction, float delta_time);
private:
	KernelParams kernel_params;
	Camera camera;
	int frame;
};