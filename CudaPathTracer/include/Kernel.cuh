#pragma once
#include <precomp.h>

struct Intersection{
	int material;
	glm::vec4 colour;
	glm::vec3 point;
	float t;

	__device__ inline static Intersection generate(int _material, glm::vec4 _colour, glm::vec3 _point, float _t){
		return{_material, _colour, _point, _t};
	}
};

struct Ray{
public:
	__device__ Ray(){};

	__device__ Ray(glm::vec3 _origin, glm::vec3 _direction, int _index){
		this->origin = _origin;
		this->direction = _direction;
		this->t = MAXDISTANCE;
		this->intersection_index = -1;
		this->cumulative_colour = glm::vec3(1.0f);
		this->pixel_index = _index;
		this->terminate_flag = false;
		this->bounce = 0;

		//initial intersection material is the background
		intersected_material = Material(0, glm::vec4(66.0f, 134.0f, 244.0f, 1.0f));
	}

	__device__ Ray(glm::vec3 _origin, glm::vec3 _direction, int _index, glm::vec3 _colour){
		this->origin = _origin;
		this->direction = _direction;
		this->t = MAXDISTANCE;
		this->intersection_index = -1;
		this->cumulative_colour = _colour;
		this->pixel_index = _index;
		this->terminate_flag = false;
		this->bounce = 0;

		//initial intersection material is the background
		this->intersected_material = Material(0, glm::vec4(66.0f, 134.0f, 244.0f, 1.0f));
	}

	glm::vec3 origin;
	glm::vec3 direction;
	float t;
	glm::vec3 intersection_point;
	int intersection_index;
	glm::vec3 cumulative_colour;
	int pixel_index;
	bool terminate_flag;
	int bounce;
	Material intersected_material;
	glm::vec3 reflected_direction;
	bool last_specular = false;
};

struct KernelParams{
	Scene scene;

	Ray* ray_buffer;
	Ray* ray_buffer_next;
	Ray* shadow_ray_buffer;
	int* ray_count;
	int* malloc_active_paths;

	int sm_cores;

	glm::vec3 top_left;
	glm::vec3 step_h;
	glm::vec3 step_v;
	glm::vec3 c_position;
};

cudaError print(cudaArray_const_t array, glm::vec4* frame_buffer, KernelParams & kernel_params, int ray_buffer_size, int frame);
cudaError launch_kernels(cudaArray_const_t array, glm::vec4* frame_buffer, KernelParams & kernel_params, BVH* bvh, int ray_buffer_size, int frame, bool n);