#pragma once
#include <precomp.h>

using namespace glm;

PathTracer::PathTracer(int cores){
	//set up variables
	this->camera = Camera();
	this->frame = 0;

	//set up kernel params
	this->kernel_params.c_position = this->camera.c_position;
	this->kernel_params.sm_cores = cores;
	this->kernel_params.scene.Init();


	//allocate memory for kernel parameters
	cudaAssert(DeviceSynchronize());

	cudaAssert(Malloc(&kernel_params.shadow_ray_buffer, BUFFERSIZE * MAXBOUNCE * sizeof(Ray)));
	cudaAssert(Malloc(&kernel_params.ray_buffer, BUFFERSIZE * sizeof(Ray)));
	cudaAssert(Malloc(&kernel_params.ray_buffer_next, BUFFERSIZE * sizeof(Ray)));
	cudaAssert(Malloc(&kernel_params.ray_count, sizeof(int)));
	cudaAssert(Malloc(&(kernel_params.malloc_active_paths), sizeof(int)))
		;

}

void inline PathTracer::CalcImageParameters(){

	/*Calculate screen points and step sizes. Doing this in here to accomodate for resizing later*/
	//Get the screen's aspect ratio
	float aspectRatio = (float)SCRWIDTH / (float)SCRHEIGHT;

	//calculate homogenous coords of this->camera point and 3 corner points of the screen in world co-ordinate space
	float topLeftx = ((2 * (0 / SCRWIDTH)) - 1) * tan((this->camera.fov / 2) * (PI / 180)) * aspectRatio;
	float topLefty = (1 - (2 * (0 / SCRHEIGHT))) * tan((this->camera.fov / 2) * (PI / 180));
	vec4 topLeft = vec4(topLeftx, topLefty, -1, 1);
	topLeft = topLeft * this->camera.c_matrix;

	float topRightx = ((2 * (SCRWIDTH / SCRWIDTH)) - 1) * tan((this->camera.fov / 2) * (PI / 180)) * aspectRatio;
	float topRighty = (1 - (2 * (0 / SCRHEIGHT))) * tan((this->camera.fov / 2) * (PI / 180));
	vec4 topRight = vec4(topRightx, topRighty, -1, 1);
	topRight = topRight * this->camera.c_matrix;

	float bottomLeftx = ((2 * (0 / SCRWIDTH)) - 1) * tan((this->camera.fov / 2) * (PI / 180)) * aspectRatio;
	float bottomLefty = (1 - (2 * (SCRHEIGHT / SCRHEIGHT))) * tan((this->camera.fov / 2) * (PI / 180));
	vec4 bottomLeft = vec4(bottomLeftx, bottomLefty, -1, 1);
	bottomLeft = bottomLeft * this->camera.c_matrix;

	//calculate x,y,x step size for pixel interpolation
	vec3 stepH = vec3(topRight.x - topLeft.x, topRight.y - topLeft.y, topRight.z - topLeft.z);
	stepH = stepH * ((float)1 / SCRWIDTH);

	vec3 stepV = vec3(bottomLeft.x - topLeft.x, bottomLeft.y - topLeft.y, bottomLeft.z - topLeft.z);
	stepV = stepV * ((float)1 / SCRHEIGHT);

	this->kernel_params.top_left = vec3(topLeft.x, topLeft.y, topLeft.z);
	this->kernel_params.step_h = stepH;
	this->kernel_params.step_v = stepV;
	this->kernel_params.c_position = vec3(this->camera.c_position.x, this->camera.c_position.y, this->camera.c_position.z);
}

void PathTracer::Trace(Renderer* cuda_interop, vec4* frame_buffer){
	CalcImageParameters();
	launch_kernels(cuda_interop->cuda_array, frame_buffer, kernel_params, BUFFERSIZE, frame);

	frame++;
}

void PathTracer::TranslateCamera(vec3 direction, float delta_time){
	mat4 translate_matrix = glm::translate(mat4(1.0f), direction * this->camera.translation_speed * delta_time);
	this->camera.c_matrix = this->camera.c_matrix * translate_matrix;
	this->camera.c_position = translate_matrix * this->camera.c_position;
}

void PathTracer::RotateCamera(vec3 direction, float delta_time){
	mat4 rotate_matrix = glm::rotate(mat4(1.0f), this->camera.rotation_theta * delta_time, direction);
	this->camera.c_matrix =  this->camera.c_matrix * rotate_matrix;
}