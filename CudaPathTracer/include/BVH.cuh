#pragma once
#include <precomp.h>
struct AABB {
public:
	void CalculateCentre();
	AABB() { 
		bounds[0] = glm::vec3(FLT_MAX); 
		bounds[1] = glm::vec3(-FLT_MAX);
		centre_point = glm::vec3(FLT_MAX);

	};
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

struct MBVHNode {
	void CalculateBounds();

	std::vector<MBVHNode*> children;
	std::vector<int> children_count;

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

struct MBVHNode_CacheFriendly {
	glm::vec3 bounds[2];
	union {
		struct {
			unsigned int child_index;
			unsigned int child_count;
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
	void ConstructBVH(std::vector<glm::vec3>& t_vertices, std::vector<int>& t_indices, int tri_count);
	void Subdivide(BVHNode* node, std::vector<AABB> AABBs, bool side, std::vector<glm::vec3>& t_vertices, std::vector<int>& t_indices);

	int ObjectPartition(BVHNode* node, std::vector<AABB> aabbs, float& splitPosition, float& split_cost);
	int SpatialPartition(BVHNode* node, std::vector<AABB> aabbs, float& split_position, float& split_cost, std::vector<glm::vec3>& t_vertices, std::vector<int>& t_indices);

	void MakeObjectSplit(std::vector<AABB>& left_AABBs, std::vector<AABB>& right_AABBs, int split_object, float split_axis_position, BVHNode* node, std::vector<AABB>& AABBs);
	void MakeSpatialSplit(std::vector<AABB>& left_AABBs, std::vector<AABB>& right_AABBs, int split_spatial, float split_axis_position, BVHNode* node, std::vector<AABB>& AABBs);

	void Collapse();

	//cache friendly
	void ConstructCacheFriendly(int tri_count);
	void PopulateCFBVH(unsigned int& current_index, unsigned int& triangle_index, BVHNode* current);
	void PopulateCFBVH(unsigned int& current_index,  unsigned int& cumulative_index, unsigned int& triangle_index, MBVHNode* current);

	//BVH variables
	BVHNode* root_node;
	MBVHNode* root_node_mbvh;
	BVHNode_CacheFriendly* cf_bvh;
	int* triangle_indices;
	MBVHNode_CacheFriendly* cf_mbvh;
	int* mbvh_triangle_indices;

	BVHNode_CacheFriendly* cf_bvh_gpu;
	int* triangle_indices_gpu;
	MBVHNode_CacheFriendly* cf_mbvh_gpu;
	int* mbvh_triangle_indices_gpu;

	int test_b = 0;
	int test_m = 0;
	//pointer to vertices and indices
	std::vector<int>* t_indices;
	std::vector<glm::vec3>* t_vertices;

	//tester
	int count_nodes(BVHNode* root);
	int count_nodes(MBVHNode* root);
	bool intersect_AABB(glm::vec3& origin, const glm::vec3& direction_inverse, const glm::vec3& aabb_min, const glm::vec3& aabb_max);
	int* intersect_MAABB(const glm::vec3& origin, const glm::vec3& direction_inverse, const float& t,
		const glm::vec4& bounds_min_x, const glm::vec4& bounds_max_x,
		const glm::vec4& bounds_min_y, const glm::vec4& bounds_max_y,
		const glm::vec4& bounds_min_z, const glm::vec4& bounds_max_z,
		glm::bvec4& result);
	bool traverse_BVH(std::vector<glm::vec3> t_vertices, std::vector<glm::vec3> t_normals, std::vector<int> t_indices, glm::vec3 origin, glm::vec3 direction, BVHNode_CacheFriendly* bvh, int* bvh_tri_list, int x, int y);
	bool traverse_MBVH(std::vector<glm::vec3> t_vertices, std::vector<glm::vec3> t_normals, std::vector<int> t_indices, glm::vec3 origin, glm::vec3 direction, float t, MBVHNode_CacheFriendly* mbvh, int* mbvh_tri_list, int x, int y);

};