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

//debug variables
__device__ int debug_count = 0;
__device__ int y_count = 0;
__device__ bool set = false;

//draw variables
surface<void, cudaSurfaceType2D> screen;


/**DEVICE FUNCTIONS**/
//get random int
__device__ unsigned int random_int(unsigned int& seed) {
	seed ^= seed << 13;
	seed ^= seed << 17;
	seed ^= seed << 5;
	return seed;
}

//get random float
__device__ float random_float(unsigned int& seed) {
	return float(random_int(seed)) * 2.3283064365387e-10f;
}

__device__ vec3 cosineWeightedSample(const vec3& normal, unsigned int& seed) {
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
__device__ vec3 sampleHemisphere(const vec3& normal, unsigned int& seed) {
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

//triangle intersection
__device__ bool intersect_triangle(const Ray& ray, const vec3& v0, const vec3& v1, const vec3& v2, vec3& intersection_point, float& t) {
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
	if (t > epsilon) {
		return true;
	}
	else
		return false;
}

__device__ bool intersect_AABB(const vec3& origin, const vec3& direction_inverse, const vec3& aabb_min, const vec3& aabb_max) {
	vec3 t1 = (aabb_min - origin) * direction_inverse;
	vec3 t2 = (aabb_max - origin) * direction_inverse;

	vec3 min = glm::min(t1, t2);
	vec3 max = glm::max(t1, t2);

	float tmin = glm::max(min.x, glm::max(min.y, min.z));
	float tmax = glm::min(max.x, glm::min(max.y, max.z));

	return (tmax >= 0.0f && tmin < tmax);
}

__device__ int* intersect_MBVH_node(const vec3& origin, const vec3& direction, const vec3& direction_inverse, const float& t,
	const vec4& bounds_min_x, const vec4& bounds_max_x,
	const vec4& bounds_min_y, const vec4& bounds_max_y,
	const vec4& bounds_min_z, const vec4& bounds_max_z,
	bvec4& result) {

	union { vec4 t_min; int t_min_i[4]; };
	vec4 t1 = (bounds_min_x - origin.x) * direction_inverse.x;
	vec4 t2 = (bounds_max_x - origin.x) * direction_inverse.x;

	t_min = glm::min(t1, t2);
	vec4 t_max = glm::max(t1, t2);

	t1 = (bounds_min_y - origin.y) * direction_inverse.y;
	t2 = (bounds_max_y - origin.y) * direction_inverse.y;

	t_min = glm::max(t_min, glm::min(t1, t2));
	t_max = glm::min(t_max, glm::max(t1, t2));

	t1 = (bounds_min_z - origin.z) * direction_inverse.z;
	t2 = (bounds_max_z - origin.z) * direction_inverse.z;

	t_min = glm::max(t_min, glm::min(t1, t2));
	t_max = glm::min(t_max, glm::max(t1, t2));

	t_min_i[0] = ((t_min_i[0] & 0xFFFFFFFC) | 0b00);
	t_min_i[1] = ((t_min_i[1] & 0xFFFFFFFC) | 0b01);
	t_min_i[2] = ((t_min_i[2] & 0xFFFFFFFC) | 0b10);
	t_min_i[3] = ((t_min_i[3] & 0xFFFFFFFC) | 0b11);

	result[0] = (t_max[0] > 0.0f) && (t_min[0] <= t_max[0]) && (t_min[0] < t);
	result[1] = (t_max[1] > 0.0f) && (t_min[1] <= t_max[1]) && (t_min[1] < t);
	result[2] = (t_max[2] > 0.0f) && (t_min[2] <= t_max[2]) && (t_min[2] < t);
	result[3] = (t_max[3] > 0.0f) && (t_min[3] <= t_max[3]) && (t_min[3] < t);

	if (t_min[0] > t_min[1])
		swap(t_min[0], t_min[1]);
	if (t_min[2] > t_min[3])
		swap(t_min[2], t_min[3]);
	if (t_min[0] > t_min[2])
		swap(t_min[0], t_min[2]);
	if (t_min[1] > t_min[3])
		swap(t_min[1], t_min[3]);
	if (t_min[2] > t_min[3])
		swap(t_min[2], t_min[3]);

	return t_min_i;
}

__device__ bool shadow_traverse_BVH(const vec3* t_vertices_gpu, const vec3* t_normals_gpu, const Material* t_mats_gpu, int* t_indices_gpu, Ray* ray, BVHNode_CacheFriendly* bvh, int* bvh_tri_list) {

	int best_index = -1;
	float shortest_distance;

	// start from infinity
	shortest_distance = FLT_MAX;

	// create a stack for each ray
	// the stack is just a fixed size array of indices to BVH nodes
	int access_stack[64];

	int stack_index = 0;
	access_stack[0] = 0;
	stack_index++;

	vec3 intersection_point;

	// while the stack is not empty
	while (stack_index) {

		//pop from stack
		BVHNode_CacheFriendly current_node = bvh[access_stack[stack_index - 1]];
		stack_index--;

		if (!(current_node.u.leaf.count & 0x80000000)) { // if inner node
			if (intersect_AABB(ray->origin, (1.0f / ray->direction), current_node.bounds[0], current_node.bounds[1])) {
				access_stack[stack_index] = current_node.u.inner.index_right;
				stack_index++;
				access_stack[stack_index] = current_node.u.inner.index_left;
				stack_index++;
				if (stack_index > 64) {
					return false;
				}
			}
		}
		else { //else if leaf node
			for (int i = current_node.u.leaf.index_first_tri; i < current_node.u.leaf.index_first_tri + (current_node.u.leaf.count & 0x7fffffff); i++) {
				float current_t;
				vec3 intersection_point;

				bool intersected_something = intersect_triangle(*ray,
					t_vertices_gpu[t_indices_gpu[bvh_tri_list[i] * 3]],
					t_vertices_gpu[t_indices_gpu[bvh_tri_list[i] * 3 + 1]],
					t_vertices_gpu[t_indices_gpu[bvh_tri_list[i] * 3 + 2]],
					intersection_point,
					current_t);

				if (!intersected_something || current_t < 0 || current_t >= ray->t) {
					continue;
				}

				ray->t = current_t;
				ray->intersection_point = intersection_point;

				vec3 normal = t_normals_gpu[bvh_tri_list[i]];
				float dot_product = dot(ray->direction, normal);
				if (dot_product > 0.0f) {
					normal *= -1.0f;
				}

				ray->t = current_t;
				ray->intersected_material = t_mats_gpu[bvh_tri_list[i]];
			}
		}

	}

	return (best_index != -1);
}

__device__ bool shadow_traverse_MBVH(const vec3* t_vertices_gpu, const vec3* t_normals_gpu, const Material* t_mats_gpu, int* t_indices_gpu, Ray* ray, MBVHNode_CacheFriendly* mbvh, int* mbvh_tri_list) {

	int best_index = -1;
	float shortest_distance;

	// start from infinity
	shortest_distance = FLT_MAX;

	// create a stack for each ray
	// the stack is just a fixed size array of indices to BVH nodes
	int access_stack[64];

	int stack_index = 0;
	access_stack[0] = 0;
	stack_index++;

	vec3 intersection_point;

	// while the stack is not empty
	while (stack_index) {

		//pop from stack
		MBVHNode_CacheFriendly current_node = mbvh[access_stack[stack_index - 1]];
		stack_index--;

		if (!(current_node.u.leaf.count & 0x80000000)) { // if inner node

			//check all children for AABB collisions
			vec4 bounds_min_x = vec4(0.0f);
			vec4 bounds_max_x = vec4(0.0f);
			vec4 bounds_min_y = vec4(0.0f);
			vec4 bounds_max_y = vec4(0.0f);
			vec4 bounds_min_z = vec4(0.0f);
			vec4 bounds_max_z = vec4(0.0f);

			for (int i = 0; i < current_node.u.inner.child_count; i++) {
				bounds_min_x[i] = mbvh[current_node.u.inner.child_index + i].bounds[0].x;
				bounds_max_x[i] = mbvh[current_node.u.inner.child_index + i].bounds[1].x;

				bounds_min_y[i] = mbvh[current_node.u.inner.child_index + i].bounds[0].y;
				bounds_max_y[i] = mbvh[current_node.u.inner.child_index + i].bounds[1].y;

				bounds_min_z[i] = mbvh[current_node.u.inner.child_index + i].bounds[0].z;
				bounds_max_z[i] = mbvh[current_node.u.inner.child_index + i].bounds[1].z;
			}
			bvec4 result = bvec4();
			int* t_min_i = intersect_MBVH_node(ray->origin, ray->direction, 1.0f / ray->direction, ray->t, bounds_min_x, bounds_max_x, bounds_min_y, bounds_max_y, bounds_min_z, bounds_max_z, result);

			//add the hits to the stack
			if (any(result)) {
				for (int i = current_node.u.inner.child_count - 1; i >= 0; i--) {
					const int idx = (t_min_i[i] & 0b11);
					if (result[i] == 1) {
						access_stack[stack_index] = current_node.u.inner.child_index + i;

						stack_index++;
						if (stack_index > 64) {
							return false;
						}
					}
				}
			}
		}
		else { //else if leaf node
			for (int i = current_node.u.leaf.index_first_tri; i < current_node.u.leaf.index_first_tri + (current_node.u.leaf.count & 0x7fffffff); i++) {
				float current_t;
				vec3 intersection_point;

				bool intersected_something = intersect_triangle(*ray,
					t_vertices_gpu[t_indices_gpu[mbvh_tri_list[i] * 3]],
					t_vertices_gpu[t_indices_gpu[mbvh_tri_list[i] * 3 + 1]],
					t_vertices_gpu[t_indices_gpu[mbvh_tri_list[i] * 3 + 2]],
					intersection_point,
					current_t);

				if (!intersected_something || current_t < 0 || current_t >= ray->t) {
					continue;
				}

				ray->t = current_t;
				ray->intersection_point = intersection_point;

				vec3 normal = t_normals_gpu[mbvh_tri_list[i]];
				float dot_product = dot(ray->direction, normal);
				if (dot_product > 0.0f) {
					normal *= -1.0f;
				}

				ray->t = current_t;
				ray->intersected_material = t_mats_gpu[mbvh_tri_list[i]];
			}
		}

	}
	return (best_index != -1);
}

__device__ bool traverse_MBVH(vec3* t_vertices_gpu, vec3* t_normals_gpu, int* t_indices_gpu, Ray* ray, MBVHNode_CacheFriendly* mbvh, int* mbvh_tri_list, unsigned int seed) {

	int best_index = -1;
	float shortest_distance;

	//start from infinity
	shortest_distance = FLT_MAX;

	//create a stack of each ray
	//the stack is just a fixed size array of indices to BVH nodes
	int access_stack[64];

	int stack_index = 0;
	access_stack[0] = 0;
	stack_index++;

	vec3 intersection_point;

	// while the stack is not empty
	while (stack_index) {
		//pop from stack
		MBVHNode_CacheFriendly current_node = mbvh[access_stack[stack_index - 1]];
		stack_index--;

		if (!(current_node.u.leaf.count & 0x80000000)) { // if inner node
			//check all children for AABB collisions
			vec4 bounds_min_x = vec4(0.0f);
			vec4 bounds_max_x = vec4(0.0f);
			vec4 bounds_min_y = vec4(0.0f);
			vec4 bounds_max_y = vec4(0.0f);
			vec4 bounds_min_z = vec4(0.0f);
			vec4 bounds_max_z = vec4(0.0f);

			for (int i = 0; i < current_node.u.inner.child_count; i++) {
				bounds_min_x[i] = mbvh[current_node.u.inner.child_index + i].bounds[0].x;
				bounds_max_x[i] = mbvh[current_node.u.inner.child_index + i].bounds[1].x;

				bounds_min_y[i] = mbvh[current_node.u.inner.child_index + i].bounds[0].y;
				bounds_max_y[i] = mbvh[current_node.u.inner.child_index + i].bounds[1].y;

				bounds_min_z[i] = mbvh[current_node.u.inner.child_index + i].bounds[0].z;
				bounds_max_z[i] = mbvh[current_node.u.inner.child_index + i].bounds[1].z;
			}
			bvec4 result = bvec4();
			int* t_min_i = intersect_MBVH_node(ray->origin, ray->direction, 1.0f / ray->direction, ray->t, bounds_min_x, bounds_max_x, bounds_min_y, bounds_max_y, bounds_min_z, bounds_max_z, result);

			//add the hits to the stack
			if (any(result)) {
				for (int i = current_node.u.inner.child_count - 1; i >= 0; i--) {
					const int idx = (t_min_i[i] & 0b11);
					if (result[i] == 1) {
						access_stack[stack_index] = current_node.u.inner.child_index + i;
						stack_index++;
						if (stack_index > 64) {
							return false;
						}
					}
				}
			}
		}
		else { //else if leaf node
			for (int i = current_node.u.leaf.index_first_tri; i < current_node.u.leaf.index_first_tri + (current_node.u.leaf.count & 0x7fffffff); i++) {
				float current_t;
				vec3 intersection_point;

				bool intersected_something = intersect_triangle(*ray,
					t_vertices_gpu[t_indices_gpu[mbvh_tri_list[i] * 3]],
					t_vertices_gpu[t_indices_gpu[mbvh_tri_list[i] * 3 + 1]],
					t_vertices_gpu[t_indices_gpu[mbvh_tri_list[i] * 3 + 2]],
					intersection_point,
					current_t);

				if (!intersected_something || current_t < 0 || current_t >= ray->t) {
					continue;
				}

				ray->t = current_t;
				ray->intersection_point = intersection_point;

				vec3 normal = t_normals_gpu[mbvh_tri_list[i]];
				float dot_product = dot(ray->direction, normal);
				if (dot_product > 0.0f) {
					normal *= -1.0f;
				}

				ray->reflected_direction = normalize(cosineWeightedSample(normal, seed));
				ray->intersection_index = mbvh_tri_list[i];
				best_index = mbvh_tri_list[i];
			}
		}
	}
	return (best_index != -1);
}

__device__ bool traverse_BVH(vec3* t_vertices_gpu, vec3* t_normals_gpu, int* t_indices_gpu, Ray* ray, BVHNode_CacheFriendly* bvh, int* bvh_tri_list, unsigned int seed, int rayN, int debug = false) {
	int best_index = -1;
	float shortest_distance;

	// start from infinity
	shortest_distance = FLT_MAX;

	// create a stack for each ray
	// the stack is just a fixed size array of indices to BVH nodes
	int access_stack[64];

	int stack_index = 0;
	access_stack[0] = 0;
	stack_index++;

	vec3 intersection_point;

	// while the stack is not empty
	while (stack_index) {
		//pop from stack
		BVHNode_CacheFriendly current_node = bvh[access_stack[stack_index - 1]];
		stack_index--;

		if (!(current_node.u.leaf.count & 0x80000000)) { // if inner node
			if (intersect_AABB(ray->origin, (1.0f / ray->direction), current_node.bounds[0], current_node.bounds[1])) {
				access_stack[stack_index] = current_node.u.inner.index_right;
				stack_index++;
				access_stack[stack_index] = current_node.u.inner.index_left;
				stack_index++;
				if (stack_index > 64) {
					return false;
				}
			}
		}
		else { //else if leaf node
			for (int i = current_node.u.leaf.index_first_tri; i < current_node.u.leaf.index_first_tri + (current_node.u.leaf.count & 0x7fffffff); i++) {
				float current_t;
				vec3 intersection_point;

				bool intersected_something = intersect_triangle(*ray,
					t_vertices_gpu[t_indices_gpu[bvh_tri_list[i] * 3]],
					t_vertices_gpu[t_indices_gpu[bvh_tri_list[i] * 3 + 1]],
					t_vertices_gpu[t_indices_gpu[bvh_tri_list[i] * 3 + 2]],
					intersection_point,
					current_t);

				if (!intersected_something || current_t < 0 || current_t >= ray->t) {
					continue;
				}

				ray->t = current_t;
				ray->intersection_point = intersection_point;

				vec3 normal = t_normals_gpu[bvh_tri_list[i]];
				float dot_product = dot(ray->direction, normal);
				if (dot_product > 0.0f) {
					normal *= -1.0f;
				}
				ray->reflected_direction = normalize(cosineWeightedSample(normal, seed));
				ray->intersection_index = bvh_tri_list[i];
				best_index = bvh_tri_list[i];
			}
		}

	}

	return (best_index != -1);
}

__device__ vec3 GetTextureColour(int index, float u, float v, vec3* texture_buffer, vec3* texture_descriptors) {
	int width = texture_descriptors[index].g;
	int height = texture_descriptors[index].b;

	float x = fmod(u, 1.0f);
	float y = fmod(v, 1.0f);

	if (x < 0) x += 1.0f;
	if (y < 0) y += 1.0f;

	int ix = int(x * (width - 1));
	int iy = int(y * (height - 1));

	int tex_buffer_index = texture_descriptors[index].r + (ix + iy * width);
	return texture_buffer[tex_buffer_index];
}
//draw to cuda surface
__device__ void Draw(vec4& colour, int x, int y) {
	surf2Dwrite(colour, screen, x * sizeof(vec4), y);
}

/**DEBUG KERNELS/FUNCTIONS**/
//generic setup for printing things
__global__ void print_helper(const Scene scene, MBVHNode_CacheFriendly* mbvh, int* mbvh_tri_list, int frame) {
	printf("frame\n");
}

/**MAIN KERNELS**/
//reset kernel variables
__global__ void ResetAllVariables() {
	counter_primary = 0;
	counter_extend = 0;
	counter_shade = 0;
	counter_connect = 0;
	start_position = 0;
	count_shadow_ray = 0;
	connect_ray_index = 0;
	debug_count = 0;
}

__global__ void SetGlobalVariables(int ray_buffer_size) {

	const unsigned int last_frame_stop = ray_buffer_size - connect_ray_index;
	start_position += last_frame_stop;
	start_position = start_position % (SCRWIDTH * SCRHEIGHT);

	counter_primary = 0;
	counter_extend = 0;
	counter_shade = 0;
	counter_connect = 0;
	connect_ray_index = 0;
	debug_count = 0;
	y_count = 0;
	count_shadow_ray = 0;
}

//process and draw each pixel colour
__global__ void draw_frame(vec4* frame_buffer, const Scene scene) {

	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;
	if (x >= SCRWIDTH || y >= SCRHEIGHT) {
		return;
	}

	const int index = x + (y * SCRWIDTH);
	vec3 temp_colour = vec3();

	//sample counter is stored in .a 
	temp_colour.r = frame_buffer[index].r / frame_buffer[index].a;
	temp_colour.g = frame_buffer[index].g / frame_buffer[index].a;
	temp_colour.b = frame_buffer[index].b / frame_buffer[index].a;

	vec3 exponent = vec3(1.0f / 2.2f);
	vec4 colour = vec4(pow(temp_colour, exponent), 1.0f);
	Draw(colour, x, y);
}

//genereate kernel
__global__ void GeneratePrimaryRays(const Scene& scene, const vec3 topLeft, const vec3 stepH, const vec3 stepV, const  vec3 c_position, Ray* ray_buffer, int ray_buffer_size, const int frame, vec4* frame_buffer) {

	while (true) {
		int index = atomicAdd(&counter_primary, 1);
		int buffer_index = index + connect_ray_index;

		if (buffer_index >= ray_buffer_size) {
			return;
		}


		const int x = (start_position + index) % SCRWIDTH;
		const int y = ((start_position + index) / SCRWIDTH) % SCRHEIGHT;

		vec3 pixelPoint = vec3(topLeft.x, topLeft.y, topLeft.z) + (stepV * (y + 0.5f)) + (stepH * (x + 0.5f));
		vec3 rayDirection = vec3(pixelPoint.x - c_position.x, pixelPoint.y - c_position.y, pixelPoint.z - c_position.z);
		vec3 rayOrigin = vec3(c_position.x, c_position.y, c_position.z);
		rayDirection = normalize(rayDirection);

		
		ray_buffer[buffer_index] = Ray(rayOrigin, rayDirection, x  + (y * SCRWIDTH));

		atomicAdd(&(frame_buffer[ray_buffer[buffer_index].pixel_index].a), 1.0f);
	}
}

//extend kernel
__global__ void Extend(const Scene scene, MBVHNode_CacheFriendly* mbvh, int* mbvh_tri_list, Ray* ray_buffer, int ray_buffer_size, int triangle_count, int frame) {
	while (true) {
		int index = atomicAdd(&counter_extend, 1);
		unsigned int seed = (index + frame * 147565741) * 720898027 * index;
		if (index >= ray_buffer_size) {
			return;
		}

		bool hit = traverse_MBVH(scene.t_vertices_gpu, scene.t_normals_gpu, scene.t_indices_gpu, &ray_buffer[index], mbvh, mbvh_tri_list, seed);
	}
}

//shade kernel
__global__ void Shade(const Scene scene, Ray* shadow_ray_buffer, Ray* ray_buffer, Ray* ray_buffer_next, int ray_buffer_size, int triangle_count, int frame, vec4* frame_buffer) {
	while (true) {
		int index = atomicAdd(&counter_shade, 1);
		unsigned int seed = (index + frame * 147565741) * 720898027 * index;

		if (index >= ray_buffer_size) {
			return;
		}

		Ray* current_ray = &ray_buffer[index];
		float max_distance = MAXDISTANCE;
		if (ray_buffer[index].t < max_distance && ray_buffer[index].bounce <= MAXBOUNCE) {
			current_ray->intersected_material = scene.t_mats_gpu[current_ray->intersection_index];
		}

		//check material, if specular then decide if this ray is a diffuse or reflected
		if (current_ray->intersected_material.type == 3) {
			float r = random_float(seed);
			if (r > current_ray->intersected_material.specularity) {
				current_ray->intersected_material.type = 2;
			}
		}

		vec3 normal;
		switch (current_ray->intersected_material.type) {
			//Background
		case 0:
			current_ray->cumulative_colour = vec3(0);
			if (current_ray->bounce == 0 || current_ray->last_specular) {
				vec2 uv = { 1.0f + atan2f(current_ray->direction.x, -current_ray->direction.z) * glm::one_over_pi<float>() * 0.5f, 1.0f - acosf(current_ray->direction.y) * glm::one_over_pi<float>() };
				int index = uv.x + (uv.y * SCRWIDTH);
				vec3 skybox_colour = GetTextureColour(0, uv.x, uv.y, scene.texture_buffer_gpu, scene.texture_descriptors_gpu);
				current_ray->cumulative_colour = skybox_colour;
			}

			current_ray->terminate_flag = true;
			break;
			//Light
		case 1:
			current_ray->cumulative_colour = vec3(0);
			if (current_ray->bounce == 0 || current_ray->last_specular) {
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
			while (counter <= scene.light_tri_count) {
				split += scene.light_areas_gpu[counter - 1];
				float proportion = split / scene.total_light_area;

				if (proportion > rand) {

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

			normal = scene.t_normals_gpu[current_ray->intersection_index];
			if (dot(current_ray->reflected_direction, normal) < 0.0f) {
				normal *= -1.0f;
			}

			vec3 BRDF = vec3(current_ray->intersected_material.colour.r * INVPI, current_ray->intersected_material.colour.g * INVPI, current_ray->intersected_material.colour.b * INVPI);
			float inv_pdf_hemisphere_sample = PI / dot(current_ray->reflected_direction, normal);

			vec3 shadow_ray_direction = random_point - current_ray->intersection_point;
			float distance_sqared = dot(shadow_ray_direction, shadow_ray_direction);
			shadow_ray_direction = normalize(shadow_ray_direction);
			vec3 shadow_ray_origin = current_ray->intersection_point + (0.0001f * shadow_ray_direction);

			vec3 light_normal = scene.t_normals_gpu[scene.tri_count - counter];
			if (dot(-1.0f * shadow_ray_direction, light_normal) < 0.0f) {
				light_normal *= -1.0f;
			}

			float n_dot_l = dot(normal, shadow_ray_direction);
			float ln_dot_l = dot(light_normal, -1.0f * shadow_ray_direction);
			if (ln_dot_l > 0 && n_dot_l > 0) {

				float area = scene.light_areas_gpu[counter - 1];
				float inverse_area_pdf = scene.total_light_area / area;
				float solid_angle = (area * (ln_dot_l)) / distance_sqared;

				float pdf1 = 1 / solid_angle;
				float pdf2 = 1 / inv_pdf_hemisphere_sample;

				float combined_pdf = ((pdf1 / (pdf1 + pdf2)) * pdf1) + ((pdf2 / (pdf1 + pdf2)) * pdf2);


				vec3 shadow_colour = BRDF * scene.emission * inverse_area_pdf * solid_angle * n_dot_l * current_ray->cumulative_colour;
				int shadow_index = atomicAdd(&count_shadow_ray, 1);
				shadow_ray_buffer[shadow_index] = Ray(shadow_ray_origin, shadow_ray_direction, current_ray->pixel_index, shadow_colour);
			}
			vec3 addition = inv_pdf_hemisphere_sample * BRDF * (dot(current_ray->reflected_direction, normal));
			current_ray->cumulative_colour *= addition;
			current_ray->last_specular = false;
			break;
		case 3:
			normal = scene.t_normals_gpu[current_ray->intersection_index];
			if (dot(current_ray->reflected_direction, normal) < 0.0f) {
				normal *= -1.0f;
			}
			current_ray->reflected_direction = reflect(current_ray->direction, normal);
			current_ray->last_specular = true;
			break;
		case 4:
			normal = scene.t_normals_gpu[current_ray->intersection_index];
			if (dot(current_ray->reflected_direction, normal) < 0.0f) {
				normal *= -1.0f;
			}
			current_ray->reflected_direction = refract(current_ray->direction, normal, current_ray->intersected_material.specularity);
			current_ray->last_specular = true;
		default:
			break;
		}

		switch (current_ray->terminate_flag) {
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
			ray_buffer_next[e_index].last_specular = current_ray->last_specular;
			break;
		default:
			break;
		}
	}
}

__global__ void ShadeReference(const Scene& scene, Ray* ray_buffer, Ray* ray_buffer_next, int ray_buffer_size, vec4* frame_buffer) {
	while (true) {
		int index = atomicAdd(&counter_shade, 1);

		if (index >= ray_buffer_size) {
			return;
		}
		Ray* current_ray = &ray_buffer[index];
		float max_distance = MAXDISTANCE;

		if (ray_buffer[index].t < max_distance && ray_buffer[index].bounce <= MAXBOUNCE) {
			current_ray->intersected_material = scene.t_mats_gpu[current_ray->intersection_index];
		}

		switch (current_ray->intersected_material.type) {
			//Background
		case 0:
			if (current_ray->bounce == 0) {
				current_ray->cumulative_colour = scene.bg_colour;
			}
			else {
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
			if (dot_product < 0.0f) {
				normal *= -1.0f;
				dot_product = dot(current_ray->reflected_direction, normal);
			}
			vec3 BRDF = vec3(current_ray->intersected_material.colour.r * INVPI, current_ray->intersected_material.colour.g * INVPI, current_ray->intersected_material.colour.b * INVPI);
			float inv_PDF = PI / dot(current_ray->reflected_direction, normal);

			vec3 addition = inv_PDF * BRDF * dot_product;
			current_ray->cumulative_colour *= addition;
			break;
		default:
			break;
		}

		switch (current_ray->terminate_flag) {
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
__global__ void Connect(const Scene scene, MBVHNode_CacheFriendly* bvh, int* bvh_tri_list, Ray* shadow_ray_buffer, int triangle_count, vec4* frame_buffer) {
	while (true) {
		int index = atomicAdd(&counter_connect, 1);

		if (index >= count_shadow_ray) {
			return;
		}

		Ray* current_ray = &shadow_ray_buffer[index];

#if 1
		shadow_traverse_MBVH(scene.t_vertices_gpu, scene.t_normals_gpu, scene.t_mats_gpu, scene.t_indices_gpu, &shadow_ray_buffer[index], bvh, bvh_tri_list);
#else
		for (int i = 0; i < triangle_count; i++) {
			float current_t;
			vec3 intersection_point;

			bool intersected_something = intersect_triangle(shadow_ray_buffer[index],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3]],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 1]],
				scene.t_vertices_gpu[scene.t_indices_gpu[i * 3 + 2]],
				intersection_point,
				current_t);

			if (!intersected_something || current_t < 0 || current_t >= current_ray->t) {
				continue;
			}

			current_ray->t = current_t;
			current_ray->intersected_material = scene.t_mats_gpu[i];
		}
#endif

		if (current_ray->intersected_material.type == 1) {
			atomicAdd(&(frame_buffer[current_ray->pixel_index].r), current_ray->cumulative_colour.r);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].g), current_ray->cumulative_colour.g);
			atomicAdd(&(frame_buffer[current_ray->pixel_index].b), current_ray->cumulative_colour.b);
		}
	}
}

//Launcher function
cudaError launch_kernels(cudaArray_const_t array, vec4* frame_buffer, KernelParams& kernel_params, BVH* bvh, int ray_buffer_size, int frame, bool new_frame) {

	cudaError err = cudaAssert(BindSurfaceToArray(screen, array));
	if (err) {
		return err;
	}

	if (new_frame) {
		ResetAllVariables << <1, 1 >> > ();
	}else{
	}

	GeneratePrimaryRays << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.top_left, kernel_params.step_h, kernel_params.step_v, kernel_params.c_position, kernel_params.ray_buffer, ray_buffer_size, frame, frame_buffer);
	SetGlobalVariables << <1, 1 >> > (ray_buffer_size);
	Extend << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, bvh->cf_mbvh_gpu, bvh->mbvh_triangle_indices_gpu, kernel_params.ray_buffer, ray_buffer_size, kernel_params.scene.tri_count, frame);
#ifdef USE_REFERENCE	
	ShadeReference << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.ray_buffer, kernel_params.ray_buffer_next, ray_buffer_size, frame_buffer);
#else
	Shade << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, kernel_params.shadow_ray_buffer, kernel_params.ray_buffer, kernel_params.ray_buffer_next, ray_buffer_size, kernel_params.scene.tri_count, frame, frame_buffer);
	Connect << <kernel_params.sm_cores * 8, 128 >> > (kernel_params.scene, bvh->cf_mbvh_gpu, bvh->mbvh_triangle_indices_gpu, kernel_params.shadow_ray_buffer, kernel_params.scene.tri_count, frame_buffer);
#endif
	//cudaAssert(DeviceSynchronize());

	std::swap(kernel_params.ray_buffer, kernel_params.ray_buffer_next);
	dim3 threads = dim3(16, 16);
	dim3 blocks = dim3((SCRWIDTH + threads.x - 1) / threads.x, (SCRHEIGHT + threads.y - 1) / threads.y);
	draw_frame << <blocks, threads >> > (frame_buffer, kernel_params.scene);
	cudaAssert(DeviceSynchronize());

	cudaError c = cudaError();
	return c;
}