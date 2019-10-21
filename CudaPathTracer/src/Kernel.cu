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
__device__ int count_shadow_ray = 0;
__device__ int connect_ray_index = 0;
__device__ int helper_count = 0;

__device__ int active_paths_gpu = 0;
//debug variables
__device__ int debug_count = 0;

//draw variables
surface<void, cudaSurfaceType2D> screen;


/**DEVICE FUNCTIONS**/

//get random int
__device__ unsigned int random_int(unsigned int& seed){
	seed ^= seed << 13;
	seed ^= seed << 17;
	seed ^= seed << 5;
	return seed;
}

//get random float
__device__ float random_float(unsigned int& seed){
	return float(random_int(seed)) * 2.3283064365387e-10f;
}

__device__ vec3 cosineWeightedSample(vec3 normal, unsigned int& seed) {
	float rand1 = random_float(seed);
	float rand2 = random_float(seed);

	const float r = sqrtf(rand1);
	const float theta = DOUBLEPI * rand2;

	const float x = r * cosf(theta);
	const float z = r * sinf(theta);

	vec3 sampleDirection = vec3(x, sqrtf(max(0.0f, 1 - rand1)), z);
	sampleDirection = normalize(sampleDirection);
	vec3 Nb;
	vec3 Nt;

	//createing coordinate systen around the normal
	if (fabs(normal.x) > fabs(normal.y)) {
		float denom = 1 / sqrtf(normal.x * normal.x + normal.z * normal.z);
		Nt = vec3(normal.z, 0, -1.0f * normal.x) * denom;
	}
	else {
		float denom = 1 / sqrtf(normal.y * normal.y + normal.z * normal.z);
		Nt = vec3(0, -1.0f * normal.z, normal.y) * denom;
	}
	Nb = cross(normal, Nt);

	Nt = vec3(normal.y * Nb.z - normal.z * Nb.y, normal.z * Nb.x - normal.x * Nb.z, normal.x * Nb.y - normal.y * Nb.x);

	vec3 result = vec3(sampleDirection.x * Nb.x + sampleDirection.y * normal.x + sampleDirection.z * Nt.x, sampleDirection.x * Nb.y + sampleDirection.y * normal.y + sampleDirection.z * Nt.y, sampleDirection.x * Nb.z + sampleDirection.y * normal.z + sampleDirection.z * Nt.z);
	result = normalize(result);
	return result;
}

//take uniform random sample from hemisphere
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

//triangle intersection
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

//draw to cuda surface
__device__ void Draw(vec4& colour, int x, int y){
	surf2Dwrite(colour, screen, x * sizeof(vec4), y);
}

/**DEBUG KERNELS/FUNCTIONS**/

//generic setup for printing things
__global__ void print_helper(Ray* ray_buffer){
	printf("%i %i \n", count_shadow_ray, connect_ray_index);
}

/**MAIN KERNELS**/

//reset kernel variables
__global__ void SetGlobalVariables(int ray_buffer_size){
	counter_primary = 0;
	counter_extend = 0;
	counter_shade = 0;
	counter_connect = 0;
	connect_ray_index = 0;
	debug_count = 0;
	count_shadow_ray = 0;
}

__global__ void SetupNextIteration(int* m_a_p){
	active_paths_gpu = connect_ray_index;
	helper_count = 0;
	*m_a_p = active_paths_gpu;
}

//process and draw each pixel colour
__global__ void draw_frame(vec4* frame_buffer){

	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;
	if (x >= SCRWIDTH || y >= SCRHEIGHT) return;

	const int index = x + (y * SCRWIDTH);
	vec3 temp_colour = vec3();

	atomicAdd(&(frame_buffer[index].a), 1.0f);

	//sample counter is stored in .a 
	temp_colour.r = frame_buffer[index].r / frame_buffer[index].a;
	temp_colour.g = frame_buffer[index].g / frame_buffer[index].a;
	temp_colour.b = frame_buffer[index].b / frame_buffer[index].a;

	vec3 exponent = vec3(1.0f / 2.2f);
	vec4 colour = vec4(pow(temp_colour, exponent), 1.0f);

	Draw(colour, y, x);
}

//genereate kernel
__global__ void GeneratePrimaryRays(Scene scene, vec3 topLeft, vec3 stepH, vec3 stepV, vec3 c_position, Ray* ray_buffer, int ray_buffer_size, int frame){

	while (true){
		int index = atomicAdd(&counter_primary, 1);
		int buffer_index = index + connect_ray_index;

		if (index >= ray_buffer_size){
			return;
		}

		const int x = (index) % SCRWIDTH;
		const int y = (index / SCRWIDTH) % SCRHEIGHT;

		vec3 pixelPoint = vec3(topLeft.x, topLeft.y, topLeft.z) + (stepH * (y + 0.5f)) + (stepV * (x + 0.5f));
		vec3 rayDirection = vec3(pixelPoint.x - c_position.x, pixelPoint.y - c_position.y, pixelPoint.z - c_position.z);
		vec3 rayOrigin = vec3(c_position.x, c_position.y, c_position.z);
		rayDirection = normalize(rayDirection);

		ray_buffer[index] = Ray(rayOrigin, rayDirection, x + (y * SCRWIDTH));
		atomicAdd(&active_paths_gpu, 1);
	}
}

//extend kernel
__global__ void Extend(Scene scene, Ray* ray_buffer, int ray_buffer_size, int triangle_count, int frame){
	while (true){
		int index = atomicAdd(&counter_extend, 1);
		unsigned int seed = (index + frame * 147565741) * 720898027 * index;
		if (index >= active_paths_gpu){
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

			float max_distance = MAXDISTANCE;

			if (ray_buffer[index].t < max_distance && ray_buffer[index].bounce <= MAXBOUNCE){
				ray_buffer[index].intersected_material = scene.t_mats_gpu[ray_buffer[index].intersection_index];
			}
		}
	}
}

//shade kernel
__global__ void Shade(Scene scene, Ray* shadow_ray_buffer, Ray* ray_buffer, Ray* ray_buffer_next, int ray_buffer_size, int triangle_count, int frame, vec4* frame_buffer){
	while (true){
		int index = atomicAdd(&counter_shade, 1);
		unsigned int seed = (index + frame * 147565741) * 720898027 * index;

		if (index >= active_paths_gpu){
			return;
		}
		Ray* current_ray = &ray_buffer[index];

		switch (current_ray->intersected_material.type){
			//Background
		case 0:
			current_ray->cumulative_colour = vec3(0);
			if (current_ray->bounce == 0){
				current_ray->cumulative_colour = scene.bg_colour;
			}

			current_ray->terminate_flag = true;
			break;
			//Light
		case 1:
			current_ray->cumulative_colour = vec3(0);
			if(current_ray->bounce == 0){
				current_ray->cumulative_colour = vec3(1.0f, 1.0f, 1.0f);
			}
			current_ray->terminate_flag = true;
			break;
			//Labertian
		case 2:
			float rand = random_float(seed);
			float split = 0;

			int counter = 1;
			vec3 random_point = vec3(0.0f);
			while (counter <= scene.light_tri_count){
				split += scene.light_areas_gpu[counter - 1];
				float proportion = split / scene.total_light_area;

				if (proportion > rand){

					//get random point on the light
					vec3 va = scene.t_vertices_gpu[scene.t_indices_gpu[scene.tri_count * 3 - (counter * 3)]];
					vec3 vb = scene.t_vertices_gpu[scene.t_indices_gpu[scene.tri_count * 3 - (counter * 3) + 1]];
					vec3 vc = scene.t_vertices_gpu[scene.t_indices_gpu[scene.tri_count * 3 - (counter * 3) + 2]];
					vec3 ab = vb - va;
					vec3 ac = vc - va;

					float w1 = random_float(seed);
					float w2 = random_float(seed);

					random_point = va + (w1 * ab) + (w2 * ac);
					break;
				}
				counter++;
			}

			vec3 BRDF = vec3(current_ray->intersected_material.colour.r * INVPI, current_ray->intersected_material.colour.g * INVPI, current_ray->intersected_material.colour.b * INVPI);
			float inv_pdf_hemisphere_sample = DOUBLEPI;


			vec3 shadow_ray_direction = random_point - current_ray->intersection_point;
			float distance_sqared = dot(shadow_ray_direction, shadow_ray_direction);
			shadow_ray_direction = normalize(shadow_ray_direction);
			vec3 shadow_ray_origin = current_ray->intersection_point + (0.0001f * shadow_ray_direction);

			vec3 normal = scene.t_normals_gpu[current_ray->intersection_index];
			if (dot(current_ray->reflected_direction, normal) < 0.0f){
				normal *= -1.0f;
			}

			vec3 light_normal = scene.t_normals_gpu[scene.tri_count - counter];
			if (dot(-1.0f * shadow_ray_direction, light_normal) < 0.0f){
				light_normal *= -1.0f;
			}

			float n_dot_l = dot(normal, shadow_ray_direction);
			float ln_dot_l = dot(light_normal, -1.0f * shadow_ray_direction);
			if (ln_dot_l > 0 && n_dot_l > 0){
			
				float area = scene.light_areas_gpu[scene.light_tri_count - counter];
				float inverse_area_pdf = scene.total_light_area / area;
				float solid_angle = (area * (ln_dot_l)) / distance_sqared;

				float pdf1 = 1 / solid_angle;
				float pdf2 = 1 / inv_pdf_hemisphere_sample;

				float combined_pdf = ((pdf1 / (pdf1 + pdf2))*pdf1) + ((pdf2 / (pdf1 + pdf2)) * pdf2);


				vec3 shadow_colour = BRDF * scene.emission * inverse_area_pdf * solid_angle * n_dot_l * current_ray->cumulative_colour;
				int shadow_index = atomicAdd(&count_shadow_ray, 1);
				shadow_ray_buffer[shadow_index] = Ray(shadow_ray_origin, shadow_ray_direction, current_ray->pixel_index, shadow_colour);
				shadow_ray_buffer[shadow_index].isShadow = true;
			}
			vec3 addition = inv_pdf_hemisphere_sample * BRDF * (dot(current_ray->reflected_direction, normal));
			current_ray->cumulative_colour *= addition;
			break;
		default:
			break;
		}

		switch(current_ray->terminate_flag){
		case true:
			atomicAdd(&(frame_buffer[current_ray->pixel_index].r), current_ray->cumulative_colour.r);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].g), current_ray->cumulative_colour.g);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].b), current_ray->cumulative_colour.b);
			break;
		case false:
			vec3 ray_origin = current_ray->intersection_point + (current_ray->reflected_direction * 0.00001f);
			int e_index = atomicAdd(&connect_ray_index, 1);
			ray_buffer_next[e_index] = Ray(ray_origin, current_ray->reflected_direction, current_ray->pixel_index, current_ray->cumulative_colour);
			ray_buffer_next[e_index].bounce = current_ray->bounce + 1;
			break;
		default:
			break;
		}
	}
}

__global__ void ShadeReference(Scene scene, Ray* ray_buffer, Ray* ray_buffer_next, int ray_buffer_size, vec4* frame_buffer){
	while (true){
		int index = atomicAdd(&counter_shade, 1);

		if (index >= active_paths_gpu){
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
			if (current_ray->bounce == 0){
				current_ray->cumulative_colour = scene.bg_colour;
			} else{
				current_ray->cumulative_colour = vec3(0.0f, 0.0f, 0.0f);
			}

			current_ray->terminate_flag = true;
			break;
		//Light
		case 1:
			current_ray->cumulative_colour = current_ray->cumulative_colour * scene.emission;
			current_ray->terminate_flag = true;
			break;
		//Labertian
		case 2:
			vec3 normal = scene.t_normals_gpu[current_ray->intersection_index];
			float dot_product = dot(current_ray->reflected_direction, normal);
			if (dot_product < 0.0f){
				normal *= -1.0f;
			}
			vec3 BRDF = vec3(current_ray->intersected_material.colour.r * INVPI, current_ray->intersected_material.colour.g * INVPI, current_ray->intersected_material.colour.b * INVPI);
			float inv_PDF = DOUBLEPI;

			vec3 addition = inv_PDF * BRDF * (dot(current_ray->reflected_direction, normal));
			current_ray->cumulative_colour *= addition;
			break;
		default:
			break;
		}

		switch (current_ray->terminate_flag){
		case true:
			atomicAdd(&(frame_buffer[current_ray->pixel_index].r), current_ray->cumulative_colour.r);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].g), current_ray->cumulative_colour.g);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].b), current_ray->cumulative_colour.b);

			atomicAdd(&debug_count, 1);

			break;
		case false:
			vec3 ray_origin = current_ray->intersection_point + (current_ray->reflected_direction * 0.00001f);
			int e_index = atomicAdd(&connect_ray_index, 1);
			ray_buffer_next[e_index] = Ray(ray_origin, current_ray->reflected_direction, current_ray->pixel_index, current_ray->cumulative_colour);
			ray_buffer_next[e_index].bounce = current_ray->bounce + 1;
			break;
		default:
			break;
		}
	}
}

//connect kernel
__global__ void Connect(Scene scene, Ray* shadow_ray_buffer, int triangle_count, vec4* frame_buffer){
	while (true){
		int index = atomicAdd(&counter_connect, 1);

		if (index >= count_shadow_ray){
			return;
		}

		Ray* current_ray = &shadow_ray_buffer[index];
		for (int i = 0; i < triangle_count; i++){
			float current_t;
			vec3 intersection_point;

			bool intersected_something = intersect_triangle(shadow_ray_buffer[index],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3]],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 1]],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 2]],
				intersection_point,
				current_t);
			
			if (!intersected_something || current_t < 0 || current_t >= current_ray->t){
				continue;
			}

			current_ray->t = current_t;
			current_ray->intersected_material = scene.t_mats_gpu[i];
		}

		if (current_ray->intersected_material.type == 1){
			atomicAdd(&(frame_buffer[current_ray->pixel_index].r), current_ray->cumulative_colour.r);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].g), current_ray->cumulative_colour.g);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].b), current_ray->cumulative_colour.b);
		}
	}
}

__global__ void g_singleAnswer(int* answer){ *answer = 2; }

//Launcher function
cudaError launch_kernels(cudaArray_const_t array, vec4* frame_buffer,  KernelParams & kernel_params, int ray_buffer_size, int frame){
	
	cudaError err = cudaAssert(BindSurfaceToArray(screen, array));
	if (err){
		return err;
	}



	int active_paths = BUFFERSIZE;
	cudaAssert(Memset(kernel_params.malloc_active_paths, 0, sizeof(int)));
	
	GeneratePrimaryRays << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.top_left, kernel_params.step_h, kernel_params.step_v, kernel_params.c_position, kernel_params.ray_buffer, ray_buffer_size, frame);
	
	while(active_paths > 0){
		SetGlobalVariables << <1, 1 >> > (ray_buffer_size);	
		Extend << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.ray_buffer, ray_buffer_size, kernel_params.scene.tri_count, frame);
#if 1		
		Shade << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.shadow_ray_buffer, kernel_params.ray_buffer, kernel_params.ray_buffer_next, ray_buffer_size, kernel_params.scene.tri_count, frame, frame_buffer);
		//print_helper << <1, 1 >> > (kernel_params.ray_buffer);
		Connect << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.shadow_ray_buffer, kernel_params.scene.tri_count, frame_buffer);
#else
		ShadeReference << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.ray_buffer, kernel_params.ray_buffer_next, ray_buffer_size, frame_buffer);
#endif
		SetupNextIteration << <1, 1 >> > (kernel_params.malloc_active_paths);
		cudaMemcpy(&active_paths, kernel_params.malloc_active_paths, sizeof(int), cudaMemcpyDeviceToHost);

		cudaAssert(DeviceSynchronize());
		std::swap(kernel_params.ray_buffer, kernel_params.ray_buffer_next);
	}

	dim3 threads = dim3(16, 16);
 	dim3 blocks = dim3((SCRWIDTH + threads.x - 1) / threads.x, (SCRHEIGHT + threads.y - 1) / threads.y);
	draw_frame << <blocks, threads >> > (frame_buffer);
	cudaAssert(DeviceSynchronize());

	cudaError c = cudaError();
	return c;
}