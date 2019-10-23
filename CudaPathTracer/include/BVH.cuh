#pragma once
#include <precomp.h>
class AABB {
public:
	void CalculateCentre();

	glm::vec3 bounds[2];
	glm::vec3 centre_point;
	glm::vec3 color = glm::vec3((rand() % 255) / 255.0f, (rand() % 255) / 255.0f, (rand() % 255) / 255.0f);
};

class BVHNode {
public:
	void CalculateBounds(std::vector<AABB> t_AABBs, int first, int count, int* indices, int tri_count);

	BVHNode* left;
	BVHNode* right;
	AABB aabb;
	bool is_leaf;
	int first;
	int count;
};

class BVH {
public:
	BVH();
	BVH(glm::vec3* t_vertices_gpu, glm::vec3* t_normals_gpu, int* t_indices_gpu);
	void ConstructBVH(std::vector<glm::vec3> t_vertices, std::vector<int> t_indices, int tri_count);
	void Subdivide(BVHNode* node, std::vector<glm::vec3> t_vertices, std::vector<int> t_indices, int tri_count);
	int Partition(BVHNode* node);

	//BVH variables
	BVHNode* root_node;
	int* indices;
	std::vector<AABB> t_AABBs;

	//pointers to scene vertices, indices and normals
	glm::vec3* t_vertices_gpu;
	glm::vec3* t_normals_gpu;
	int* t_indices_gpu;
};