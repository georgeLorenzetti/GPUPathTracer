#include <precomp.h>

using namespace glm;

void Scene::cornell_setup(){
	//scene
	this->bg_colour = vec3(66, 134, 244);

	this->t_vertices = {

		//floor
		vec3(-3.5f, -3.0f, -8.0f), //0
		vec3(3.5f, -3.0f, -8.0f),  //1
		vec3(-3.5f, -3.0f, -3.0f), //2
		vec3(3.5f, -3.0f, -3.0f),  //3

		//ceiling
		vec3(-3.5f, 3.0f, -8.0f), //4
		vec3(3.5f, 3.0f, -8.0f),  //5
		vec3(-3.5f, 3.0f, -3.0f), //6
		vec3(3.5f, 3.0f, -3.0f),  //7

		//left wall
		vec3(-3.5f, 3.0f, -8.0f), //8
		vec3(-3.5f, -3.0f, -8.0f),//9
		vec3(-3.5f, 3.0f, -3.0f), //10
		vec3(-3.5f, -3.0f, -3.0f),//11

		////right wall
		vec3(3.5f, 3.0f, -8.0f), //12
		vec3(3.5f, -3.0f, -8.0f),//13
		vec3(3.5f, 3.0f, -3.0f), //14
		vec3(3.5f, -3.0f, -3.0f),//15

		//back wall
		vec3(-3.5f, 3.0f, -8.0f), //16
		vec3(-3.5f, -3.0f, -8.0f),//17
		vec3(3.5f, 3.0f, -8.0f), //18
		vec3(3.5f, -3.0f, -8.0f),//19

		//big box
		vec3(0.0f, -3.0f, -6.0f), //20
		vec3(0.0f, -3.0f, -5.0f),//21
		vec3(-1.0f, -3.0f, -6.0f), //22
		vec3(-1.0f, -3.0f, -5.0f),//23
		vec3(0.0f, -1.0f, -6.0f), //24
		vec3(0.0f, -1.0f, -5.0f),//25
		vec3(-1.0f, -1.0f, -6.0f), //26
		vec3(-1.0f, -1.0f, -5.0f),//27

		//lights
		vec3(0.5f, 2.9f, -4.0f), //28
		vec3(0.5f, 2.9f, -5.0f), //29
		vec3(-0.5f, 2.9f, -4.0f), //30
		vec3(-0.5f, 2.9f, -5.0f), //31
	};

	this->t_indices = {

		//floor
		0, 1, 2,
		2, 1, 3,

		//ceiling
		4, 5, 6,
		6, 5, 7,

		//left wall
		8, 9, 10,
		10, 9, 11,

		//right wall
		12, 13, 14,
		14, 13, 15,

		//back wall
		16, 17, 18,
		18, 17, 19,

		////big box
		21, 22, 23,
		21, 20, 22,

		21, 20, 25,
		20, 24, 25,

		27, 25, 21,
		27, 23, 21,

		26, 24, 20,
		26, 24, 22,

		27, 26, 22,
		27, 23, 22,

		26, 24, 25,
		26, 27, 25,


		//lights
		28, 29, 30,
		30, 29, 31
	};

	this->t_mats = std::vector<Material>();

	//*scene triangles*//
	//floor
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));

	//ceiling
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));

	//left wall
	this->t_mats.push_back(Material(2, vec3(255.0f, 0.0f, 0.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 0.0f, 0.0f)));

	//right wall
	this->t_mats.push_back(Material(2, vec3(0.0f, 255.0f, 0.0f)));
	this->t_mats.push_back(Material(2, vec3(0.0f, 255.0f, 0.0f)));

	//back wall
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));

	//big box
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));
	this->t_mats.push_back(Material(2, vec3(255.0f, 255.0f, 255.0f)));

	//lights
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));

	this->t_normals = std::vector<vec3>();
	for (int i = 0; i < t_mats.size(); i++){
		vec3 vector1 = (t_vertices[t_indices[i * 3 + 1]] - t_vertices[t_indices[i * 3]]);
		vec3 vector2 = (t_vertices[t_indices[i * 3 + 2]] - t_vertices[t_indices[i * 3]]);
		vec3 normal = cross(vector1, vector2);
		normal = normalize(normal);
		normal *= -1;
		t_normals.push_back(normal);
	}
}

void Scene::init(){
	cornell_setup();
	this->tri_count = this->t_mats.size();
	std::cout << tri_count << std::endl;
	cudaAssert(Malloc(&(this->t_vertices_gpu), this->t_vertices.size() * sizeof(vec3)));
	cudaAssert(Malloc(&(this->t_indices_gpu), this->t_indices.size() * sizeof(int)));
	cudaAssert(Malloc(&(this->t_mats_gpu), this->t_mats.size() * sizeof(Material)));
	cudaAssert(Malloc(&(this->t_normals_gpu), this->t_normals.size() * sizeof(vec3)));

	cudaAssert(Memcpy(this->t_vertices_gpu, this->t_vertices.data(), this->t_vertices.size() * sizeof(vec3), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->t_indices_gpu, this->t_indices.data(), this->t_indices.size() * sizeof(int), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->t_mats_gpu, this->t_mats.data(), this->t_mats.size() * sizeof(Material), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->t_normals_gpu, this->t_normals.data(), this->t_normals.size() * sizeof(vec3), cudaMemcpyHostToDevice));

}
