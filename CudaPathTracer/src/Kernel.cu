#pragma once
#include <precomp.h>
#include <surface_functions.h>
using namespace glm;

//kernel variables
__device__ int counter_primary = 0;
__device__ int counter_extend = 0;
__device__ int counter_shade = 0;
__device__ int counter_connect = 0;
__device__ int start_position = 0;
__device__ int connect_ray_index = 0;

__device__ int debug_count = 0;

//draw variables
surface<void, cudaSurfaceType2D> screen;


//Device functions
__device__ unsigned int random_int(unsigned int& seed){
	seed ^= seed << 13;
	seed ^= seed << 17;
	seed ^= seed << 5;
	return seed;
}

__device__ float random_float(unsigned int& seed){
	return float(random_int(seed)) * 2.3283064365387e-10f;
}

__device__ vec3 sampleHemisphere(vec3 normal, unsigned int& seed){
	float rand1 = random_float(seed);
	float rand2 = random_float(seed);

	float sinTheta = sqrtf(1 - rand1 * rand1);
	float phi = DOUBLEPI * rand2;
	float x = sinTheta * cosf(phi);
	float z = sinTheta * sinf(phi);
	vec3 sampleDirection = vec3(x, rand1, z);
	vec3 Nb;
	vec3 Nt;

	//creating coordinate systen around the normal
	if (fabs(normal.x) > fabs(normal.y)){
		float denom = 1 / sqrtf(normal.x * normal.x + normal.z * normal.z);
		Nt = vec3(normal.z, 0, -1.0f * normal.x) * denom;
	} else{
		float denom = 1 / sqrtf(normal.y * normal.y + normal.z * normal.z);
		Nt = vec3(0, -1.0f * normal.z, normal.y) * denom;
	}
	Nb = cross(normal, Nt);

	Nt = vec3(normal.y * Nb.z - normal.z * Nb.y, normal.z * Nb.x - normal.x * Nb.z, normal.x * Nb.y - normal.y * Nb.x);

	vec3 result = vec3(sampleDirection.x * Nb.x + sampleDirection.y * normal.x + sampleDirection.z * Nt.x, sampleDirection.x * Nb.y + sampleDirection.y * normal.y + sampleDirection.z * Nt.y, sampleDirection.x * Nb.z + sampleDirection.y * normal.z + sampleDirection.z * Nt.z);
	result = normalize(result);
	return result;
}

__device__ bool intersect_triangle(Ray ray, vec3 v0, vec3 v1, vec3 v2, vec3& intersection_point, float& t){
	const float epsilon = 0.0000001;
	vec3 edge1, edge2, h, s, q;
	float a, f, u, v;
	edge1 = v1 - v0;
	edge2 = v2 - v0;
	h = cross(ray.direction, edge2);
	a = dot(edge1, h);

	if (a > -epsilon && a < epsilon)
		return false;

	f = 1.0 / a;
	s = ray.origin - v0;
	u = f * (dot(s, h));
	if (u < 0.0 || u > 1.0)
		return false;

	q = cross(s, edge1);
	v = f * dot(ray.direction, q);
	if (v < 0.0 || u + v > 1.0)
		return false;

	t = f * dot(edge2, q);
	intersection_point = ray.origin + ray.direction * t;
	if (t > epsilon){
		return true;
	} else
		return false;
}

//Debug Kernel
__global__ void print_frame(vec4* framebf, Ray* ray_buffer_next, Scene scene, int frame){
	//if(frame >=2 ){
	//for (int j = 0; j < connect_ray_index; j++){
	//	for (int i = 0; i < 4; i++){
	//		float current_t;
	//		vec3 intersection_point;

	//		bool intersected_something = intersect_triangle(ray_buffer_next[j],
	//			scene.t_vertices_gpu[scene.t_indices_gpu[i * 3]],
	//			scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 1]],
	//			scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 2]],
	//			intersection_point,
	//			current_t);

	//		if (i > 1 && intersected_something){
	//			printf("Yes \n");
	//		} else{
	//			//printf("no \n");
	//		}

	//	}
	//}
	//}

}

//Main Kernels
__global__ void set_global_variables(int ray_buffer_size){
	const unsigned int last_frame_stop = ray_buffer_size - counter_primary;
	start_position += last_frame_stop;
	start_position = start_position % (SCRWIDTH * SCRHEIGHT);

	counter_primary = 0;
	counter_extend = 0;
	counter_shade = 0;
	counter_connect = 0;
	connect_ray_index = 0;
	debug_count = 0;
}

__device__ void draw(vec4& colour, int x, int y){
	surf2Dwrite(colour, screen, x * sizeof(vec4), y);
}

__global__ void draw_frame(vec4* frame_buffer){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;
	if (x >= SCRWIDTH || y >= SCRHEIGHT) return;

	const int index = x + (y * SCRWIDTH);
	vec3 temp_colour = vec3();
	temp_colour.r = frame_buffer[index].r / frame_buffer[index].a;
	temp_colour.g = frame_buffer[index].g / frame_buffer[index].a;
	temp_colour.b = frame_buffer[index].b / frame_buffer[index].a;

	vec4 colour = vec4(temp_colour, 1.0f);
	draw(colour, y, x);
}

__global__ void colour(vec4* frame_buffer, int frame){

	int x = (blockIdx.x * blockDim.x) + threadIdx.x;
	int y = (blockIdx.y * blockDim.y) + threadIdx.y;
	int index = x + (y * SCRWIDTH);

	if (x >= SCRWIDTH || y >= SCRHEIGHT){
		return;
	}
	unsigned int seed = (index + frame * 147565741) * 720898027 * index;
	frame_buffer[index] = vec4(random_float(seed), random_float(seed), random_float(seed), 1.0f);
}

__global__ void generatePrimaryRays(Scene scene, vec3 topLeft, vec3 stepH, vec3 stepV, vec3 c_position, Ray* ray_buffer, int ray_buffer_size, int frame){

	while (true){
		int index = atomicAdd(&counter_primary, 1);
		int buffer_index = index + connect_ray_index;

		if (buffer_index >= ray_buffer_size){
			return;
		}

		const int x = (start_position + index) % SCRWIDTH;
		const int y = ((start_position + index) / SCRWIDTH) % SCRHEIGHT;

		vec3 pixelPoint = vec3(topLeft.x, topLeft.y, topLeft.z) + (stepH * (y + 0.5f)) + (stepV * (x + 0.5f));
		vec3 rayDirection = vec3(pixelPoint.x - c_position.x, pixelPoint.y - c_position.y, pixelPoint.z - c_position.z);
		vec3 rayOrigin = vec3(c_position.x, c_position.y, c_position.z);
		rayDirection = normalize(rayDirection);

		
		ray_buffer[buffer_index] = Ray(rayOrigin, rayDirection, x + (y * SCRWIDTH));

	}
}

__global__ void extend(Scene scene, Ray* ray_buffer, int ray_buffer_size, int triangle_count, int frame){
	while (true){
		int index = atomicAdd(&counter_extend, 1);
		unsigned int seed = (index + frame * 147565741) * 720898027 * index;
		if (index >= ray_buffer_size){
			return;
		}

		//no acceleration structure yet
		for (int i = 0; i < triangle_count; i++){
			float current_t;
			vec3 intersection_point;

			bool intersected_something = intersect_triangle(ray_buffer[index],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3]],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 1]],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 2]],
				intersection_point,
				current_t);

			if (!intersected_something || current_t < 0 || current_t >= ray_buffer[index].t){
				continue;
			}

			ray_buffer[index].t = current_t;
			ray_buffer[index].intersection_point = intersection_point;

			vec3 normal = scene.t_normals_gpu[i];
			float dot_product = dot(ray_buffer[index].direction, normal);
			if (dot_product > 0.0f){
				normal *= -1.0f;
			}
			ray_buffer[index].reflected_direction = normalize(sampleHemisphere(normal, seed));
			ray_buffer[index].intersection_index = i;
		}
	}
}

__global__ void shade(Scene scene, Ray* ray_buffer, int ray_buffer_size){
	while (true){
		int index = atomicAdd(&counter_shade, 1);

		if (index >= ray_buffer_size){
			return;
		}
		Ray* current_ray = &ray_buffer[index];
		float max_distance = MAXDISTANCE;

		if (ray_buffer[index].t < max_distance && ray_buffer[index].bounce <= MAXBOUNCE){
			current_ray->intersected_material = scene.t_mats_gpu[current_ray->intersection_index];
		}

		switch (current_ray->intersected_material.type){
			//Background
		case 0:
			current_ray->cumulative_colour = vec3(0.0f, 0.0f, 0.0f);
			current_ray->terminate_flag = true;
			break;
			//Light
		case 1:
			current_ray->cumulative_colour *= scene.emission;
			current_ray->terminate_flag = true;
			break;
			//Labertian
		case 2:
			vec3 normal = scene.t_normals_gpu[current_ray->intersection_index];
			float dot_product = dot(current_ray->direction, normal);
			if (dot_product > 0.0f){
				normal *= -1.0f;
			}
			vec3 BRDF = vec3(current_ray->intersected_material.colour.r * INVPI, current_ray->intersected_material.colour.g * INVPI, current_ray->intersected_material.colour.b * INVPI);
			vec3 inv_PDF = vec3(DOUBLEPI, DOUBLEPI, DOUBLEPI);
			current_ray->cumulative_colour = inv_PDF * BRDF * current_ray->cumulative_colour * (dot(current_ray->reflected_direction, normal));
			break;
		default:
			break;
		}
	}
}

__global__ void connect(Scene scene, Ray* ray_buffer, Ray* ray_buffer_next, int ray_buffer_size, int triangle_count, vec4* frame_buffer){
	while (true){
		int index = atomicAdd(&counter_connect, 1);

		if (index >= ray_buffer_size){
			return;
		}

		Ray* current_ray = &ray_buffer[index];
		switch (current_ray->terminate_flag){
		case true:
			//clamp to 255
			if (current_ray->cumulative_colour.r > 255.0f)
				current_ray->cumulative_colour.r = 255.0f;
			if (current_ray->cumulative_colour.g > 255.0f)
				current_ray->cumulative_colour.g = 255.0f;
			if (current_ray->cumulative_colour.b > 255.0f)
				current_ray->cumulative_colour.b = 255.0f;

			atomicAdd(&frame_buffer[current_ray->pixel_index].r, current_ray->cumulative_colour.r);
			atomicAdd(&frame_buffer[current_ray->pixel_index].g, current_ray->cumulative_colour.g);
			atomicAdd(&frame_buffer[current_ray->pixel_index].b, current_ray->cumulative_colour.b);
			atomicAdd(&frame_buffer[current_ray->pixel_index].a, 1.0f);
			atomicAdd(&debug_count, 1);
			break;
		case false:
			vec3 ray_origin = current_ray->intersection_point + (current_ray->reflected_direction * 0.00001f);
			int e_index = atomicAdd(&connect_ray_index, 1);
			ray_buffer_next[e_index] = Ray(ray_origin, current_ray->reflected_direction, current_ray->pixel_index, current_ray->cumulative_colour);
			ray_buffer_next[e_index].bounce = current_ray->bounce + 1;
			atomicAdd(&frame_buffer[current_ray->pixel_index].a, 1.0f);
			break;
		default:
			break;
		}
	}
}

cudaError print(cudaArray_const_t array, vec4* frame_buffer, KernelParams & kernel_params, int ray_buffer_size, int frame){
	print_frame << <1, 1 >> > (frame_buffer, kernel_params.ray_buffer, kernel_params.scene, frame);
	cudaAssert(DeviceSynchronize());
	cudaError c = cudaError();
	return c;
}

cudaError launch_kernels(cudaArray_const_t array, vec4* frame_buffer,  KernelParams & kernel_params, int ray_buffer_size, int frame){
	
	cudaError err = cudaAssert(BindSurfaceToArray(screen, array));
	if (err){
		return err;
	}

	generatePrimaryRays << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.top_left, kernel_params.step_h, kernel_params.step_v, kernel_params.c_position, kernel_params.ray_buffer, ray_buffer_size, frame);
	set_global_variables << <1, 1 >> > (ray_buffer_size);
	extend << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.ray_buffer, ray_buffer_size, kernel_params.scene.tri_count, frame);
	shade << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.ray_buffer, ray_buffer_size);
	connect << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.ray_buffer, kernel_params.ray_buffer_next, ray_buffer_size, kernel_params.scene.tri_count, frame_buffer);
	cudaAssert(DeviceSynchronize());

	std::swap(kernel_params.ray_buffer, kernel_params.ray_buffer_next);

	//const dim3 blockSize2d(8, 8);
	//const dim3 blocksPerGrid2d(
	//	(SCRWIDTH + blockSize2d.x - 1) / blockSize2d.x,
	//	(SCRHEIGHT + blockSize2d.y - 1) / blockSize2d.y);
	//colour << <blocksPerGrid2d, blockSize2d >> > (frame_buffer, frame);

	dim3 threads = dim3(16, 16);
 	dim3 blocks = dim3((SCRWIDTH + threads.x - 1) / threads.x, (SCRHEIGHT + threads.y - 1) / threads.y);
	draw_frame << <blocks, threads >> > (frame_buffer);
	cudaAssert(DeviceSynchronize());

	cudaError c = cudaError();
	return c;
}