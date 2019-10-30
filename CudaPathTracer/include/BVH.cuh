#pragma once
#include <precomp.h>
struct AABB {
public:
	void CalculateCentre();

	glm::vec3 bounds[2]; // 1 = max, 0 = min
	glm::vec3 centre_point;
	int index = -1;
};

struct BVHNode {
public:
	void CalculateBounds(std::vector<AABB> t_AABBs);

	BVHNode* left;
	BVHNode* right;

	AABB aabb;

	std::vector<int> leaf_triangles;
	bool is_leaf;
};

//32 bit struct as described here:http://raytracey.blogspot.com/search?q=gpu+bvh
//storing bvh as a bunch of these structs is nice for chaching very cool very fast
struct BVHNode_CacheFriendly {
public:
	glm::vec3 bounds[2];

	union {
		struct {
			unsigned int index_left;
			unsigned int index_right;
		} inner;
		// leaf node: stores triangle count and starting index in triangle list
		struct {
			unsigned int count; // Top-most bit set, leafnode if set, innernode otherwise
			unsigned int index_first_tri;
		} leaf;
	} u;
};

class BVH {
public:
	//main bvh construction
	void ConstructBVH(std::vector<glm::vec3> t_vertices, std::vector<int> t_indices, int tri_count);
	void Subdivide(BVHNode* node, std::vector<AABB> AABBs, bool side);
	int Partition(BVHNode* node, std::vector<AABB> aabbs, float& splitPosition);

	//cache friendly
	void ConstructCacheFriendly(int tri_count);
	void PopulateCFBVH(unsigned int& current_index, unsigned int& triangle_index, BVHNode* current);

	//BVH variables
	BVHNode* root_node;
	BVHNode_CacheFriendly* cf_bvh;
	int* triangle_indices;

	BVHNode_CacheFriendly* cf_bvh_gpu;
	int* triangle_indices_gpu;
};