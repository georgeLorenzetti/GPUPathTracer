#pragma once
#include <precomp.h>
class Scene{
	public:
		void Init();
		//pointers to vectors for CUDA
		Material* t_mats_gpu;
		glm::vec3* t_vertices_gpu;
		glm::vec3* t_normals_gpu;
		int* t_indices_gpu;
		float* light_areas_gpu;


		//extra public variables **should probs move emission to Material**
		int tri_count;
		int light_tri_count;
		float total_light_area;
		float emission = 10.0f;
		glm::vec3 bg_colour;

	private:
		//set up the cornell box scene
		void CornellSetup();

		//lists for all the triangle data
		std::vector<glm::vec3> t_vertices;
		std::vector<int> t_indices;
		std::vector<Material> t_mats;
		std::vector<glm::vec3> t_normals;
		std::vector<float> light_areas;

};

