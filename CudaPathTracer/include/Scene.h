#pragma once
#include <precomp.h>
class Scene{
	public:
		void init();
		//pointers to vectors for CUDA
		glm::vec3* t_vertices_gpu;
		int* t_indices_gpu;
		Material* t_mats_gpu;
		glm::vec3* t_normals_gpu;
		int tri_count;
		float emission = 100.0f;

	private:
		void cornell_setup();
		//lists for all the triangle data
		glm::vec3 bg_colour;
		std::vector<glm::vec3> t_vertices;
		std::vector<int> t_indices;
		std::vector<Material> t_mats;
		std::vector<glm::vec3> t_normals;
};

