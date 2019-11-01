#include <precomp.h>
#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>

using namespace glm;

void Scene::CornellSetup(){
	//scene
	this->bg_colour = vec3(66.0f/255.0f, 134.0f/255.0f, 244.0f/255.0f);

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
		vec3(-0.0f, -3.0f, -7.0f), //20
		vec3(-0.0f, -3.0f, -5.0f), //21
		vec3(-2.0f, -3.0f, -7.0f),//22
		vec3(-2.0f, -3.0f, -5.0f),//23
		vec3(-0.0f, -1.0f, -7.0f), //24
		vec3(-0.0f, -1.0f, -5.0f), //25
		vec3(-2.0f, -1.0f, -7.0f),//26
		vec3(-2.0f, -1.0f, -5.0f),//27

		//small box
		vec3(1.0f, -3.0f, -5.0f), //28
		vec3(1.5f, -3.0f, -4.5f), //29
		vec3(2.0f, -3.0f, -5.0f),//30
		vec3(1.5f, -3.0f, -5.5f),//31
		vec3(1.0f, -2.0f, -5.0f), //32
		vec3(1.5f, -2.0f, -4.5f), //33
		vec3(2.0f, -2.0f, -5.0f),//34
		vec3(1.5f, -2.0f, -5.5f),//35

		//lights
		vec3(1.5f, 2.99999f, -5.0f), //36
		vec3(1.5f, 2.99999f, -6.5f), //37
		vec3(-1.5f, 2.99999f, -5.0f), //38
		vec3(-1.5f, 2.99999f, -6.5f), //39
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
		27, 25, 23,
		23, 21, 25,

		25, 24, 20,
		25, 21, 20,

		26, 24, 22,
		24, 20, 22,

		27, 26, 22,
		27, 23, 22,

		26, 24, 25,
		26, 27, 25,
		
		22, 23, 21,
		23, 21, 20,

		//small box
		28, 29, 33,
		28, 32, 33,

		33, 34, 30,
		33, 29, 30,

		35, 34, 30,
		35, 31, 30,

		32, 35, 31,
		31, 28, 32,

		32, 33, 34,
		32, 35, 34,

		28, 29, 30,
		30, 31, 28,

		//lights
		36, 37, 38,
		38, 37, 39
	};

	this->t_mats = std::vector<Material>();

	//*scene triangles*//
	//floor
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));

	//ceiling
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));

	//left wall
	this->t_mats.push_back(Material(2, vec3(1.0f, 0.0f, 0.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 0.0f, 0.0f)));

	//right wall
	this->t_mats.push_back(Material(2, vec3(0.0f, 1.0f, 0.0f)));
	this->t_mats.push_back(Material(2, vec3(0.0f, 1.0f, 0.0f)));

	//back wall
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));

	//big box
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));

	//small box
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));
	this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 1.0f)));

	//lights
	this->light_tri_count = 0;
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->light_tri_count++;
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->light_tri_count++;

	this->t_normals = std::vector<vec3>();
	for (int i = 0; i < t_mats.size(); i++){
		vec3 vector1 = (t_vertices[t_indices[i * 3 + 1]] - t_vertices[t_indices[i * 3]]);
		vec3 vector2 = (t_vertices[t_indices[i * 3 + 2]] - t_vertices[t_indices[i * 3]]);
		vec3 normal = cross(vector1, vector2);
		normal = normalize(normal);
		normal *= -1;
		t_normals.push_back(normal);
	}

	this->light_areas = std::vector<float>();
	this->tri_count = this->t_mats.size();
	this->total_light_area = 0;
	for (int i = 1; i <= light_tri_count; i++){
		vec3 va = t_vertices[t_indices[t_mats.size()*3 - (i*3)]];
		vec3 vb = t_vertices[t_indices[t_mats.size()*3 - (i*3) + 1]];
		vec3 vc = t_vertices[t_indices[t_mats.size()*3 - (i*3) + 2]];

		vec3 ab = vb - va;
		vec3 ac = vc - va;

		vec3 cr = cross(ab, ac);
		float area = length(cr) * 0.5f;
		this->light_areas.push_back(area);
		this->total_light_area += area;
	}
}

void Scene::DeerSceneSetup() {
	this->bg_colour = vec3(66.0f / 255.0f, 134.0f / 255.0f, 244.0f / 255.0f);
	this->tri_count = 0;

	this->t_vertices = std::vector<vec3>();
	this->t_indices = std::vector<int>();
	this->t_normals = std::vector<vec3>();
	this->t_mats = std::vector<Material>();
	LoadObject("assets/dragon.obj");

	for (int i = 0; i < this->tri_count; i++) {
		this->t_mats.push_back(Material(2, vec3(1.0f, 0.0f, 1.0f)));
	}

	//floor
	//this->t_indices.push_back(this->t_vertices.size());
	//this->t_vertices.push_back(vec3(1.0f, -1.0f, -2.0f));
	//this->t_indices.push_back(this->t_vertices.size());
	//this->t_vertices.push_back(vec3(1.0f, -1.0f, 0.0f));
	//this->t_indices.push_back(this->t_vertices.size());
	//this->t_vertices.push_back(vec3(-1.0f, -1.0f, -2.0f));
	//this->t_mats.push_back(Material(2, vec3(1.0f, 0.0f, 0.0f)));
	//this->tri_count++;

	//this->t_indices.push_back(this->t_vertices.size());
	//this->t_vertices.push_back(vec3(-1.0f, -1.0f, -2.0f));
	//this->t_indices.push_back(this->t_vertices.size());
	//this->t_vertices.push_back(vec3(-1.0f, -1.0f, 0.0f));
	//this->t_indices.push_back(this->t_vertices.size());
	//this->t_vertices.push_back(vec3(1.0f, -1.0f, 0.0f));
	//this->t_mats.push_back(Material(2, vec3(1.0f, 1.0f, 0.0f)));
	//this->tri_count++;

	//lights
	this->light_tri_count = 0;

	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(1.5f, 3.99999f, -5.0f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(1.5f, 3.99999f, -1.5f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(-1.5f, 3.99999f, -5.0f));
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->light_tri_count++;
	this->tri_count++;

	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(-1.5f, 3.99999f, -5.0f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(1.5f, 3.99999f, -1.5f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(-1.5f, 3.99999f, -1.5f));
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->light_tri_count++;
	this->tri_count++;

	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(1.5f, 3.99999f,  5.0f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(1.5f, 3.99999f, 1.5f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(-1.5f, 3.99999f, 5.0f));
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->light_tri_count++;
	this->tri_count++;

	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(-1.5f, 3.99999f, 5.0f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(1.5f, 3.99999f, 1.5f));
	this->t_indices.push_back(this->t_vertices.size());
	this->t_vertices.push_back(vec3(-1.5f, 3.99999f, 1.5f));
	this->t_mats.push_back(Material(1, vec3(0.0f, 0.0f, 0.0f)));
	this->light_tri_count++;
	this->tri_count++;

	this->t_normals = std::vector<vec3>();
	for (int i = 0; i < t_mats.size(); i++) {
		vec3 vector1 = (t_vertices[t_indices[i * 3 + 1]] - t_vertices[t_indices[i * 3]]);
		vec3 vector2 = (t_vertices[t_indices[i * 3 + 2]] - t_vertices[t_indices[i * 3]]);
		vec3 normal = cross(vector1, vector2);
		normal = normalize(normal);
		normal *= -1;
		t_normals.push_back(normal);
	}

	this->light_areas = std::vector<float>();
	this->total_light_area = 0;
	for (int i = 1; i <= light_tri_count; i++) {
		vec3 va = t_vertices[t_indices[t_mats.size() * 3 - (i * 3)]];
		vec3 vb = t_vertices[t_indices[t_mats.size() * 3 - (i * 3) + 1]];
		vec3 vc = t_vertices[t_indices[t_mats.size() * 3 - (i * 3) + 2]];

		vec3 ab = vb - va;
		vec3 ac = vc - va;

		vec3 cr = cross(ab, ac);
		float area = length(cr) * 0.5f;
		this->light_areas.push_back(area);
		this->total_light_area += area;
	}

}

void Scene::LoadObject(std::string filename) {
	tinyobj::attrib_t attribute;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn, err;
	std::string path = filename;

	if (!tinyobj::LoadObj(&attribute, &shapes, &materials, &warn, &err, path.c_str())) {
		throw std::runtime_error(warn + err);
	}

	std::vector<tinyobj::index_t> ind = shapes[0].mesh.indices;
	for (int i = 0; i < ind.size(); i++) {
		this->t_indices.push_back(ind[i].vertex_index);
	}

	this->tri_count = ind.size() / 3;
	for (int i = 0; i < attribute.vertices.size() / 3; i++) {
		tinyobj::real_t vx = attribute.vertices[3 * i + 0] * 0.5f;
		tinyobj::real_t vy = attribute.vertices[3 * i + 1] * 0.5f;
		tinyobj::real_t vz = attribute.vertices[3 * i + 2] * 0.5f;

		vy -= 1.0f;
		vz -= 1.0f;
		this->t_vertices.push_back(vec3((float)vx, (float)vy, (float)vz));
	}
}

void Scene::Init(){
	DeerSceneSetup();
	//CornellSetup();
	this->tri_count = this->t_mats.size();

	cudaAssert(Malloc(&(this->t_vertices_gpu), this->t_vertices.size() * sizeof(vec3)));
	cudaAssert(Malloc(&(this->t_indices_gpu), this->t_indices.size() * sizeof(int)));
	cudaAssert(Malloc(&(this->t_mats_gpu), this->t_mats.size() * sizeof(Material)));
	cudaAssert(Malloc(&(this->t_normals_gpu), this->t_normals.size() * sizeof(vec3)));
	cudaAssert(Malloc(&(this->light_areas_gpu), this->light_areas.size() * sizeof(float)));

	cudaAssert(Memcpy(this->t_vertices_gpu, this->t_vertices.data(), this->t_vertices.size() * sizeof(vec3), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->t_indices_gpu, this->t_indices.data(), this->t_indices.size() * sizeof(int), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->t_mats_gpu, this->t_mats.data(), this->t_mats.size() * sizeof(Material), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->t_normals_gpu, this->t_normals.data(), this->t_normals.size() * sizeof(vec3), cudaMemcpyHostToDevice));
	cudaAssert(Memcpy(this->light_areas_gpu, this->light_areas.data(), this->light_areas.size() * sizeof(float), cudaMemcpyHostToDevice));
}

void Scene::ColourBVH(BVHNode* node) {
	
	if (node->is_leaf) {
		vec3 colour = vec3(static_cast <float> (rand()) / static_cast <float> (RAND_MAX), static_cast <float> (rand()) / static_cast <float> (RAND_MAX), static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		for (int i = 0; i < node->leaf_triangles.size(); i++) {
			t_mats[node->leaf_triangles[i]].colour = colour;
		}
		return;
	}

	ColourBVH(node->left);
	ColourBVH(node->right);
}

void Scene::ColourBVH2(BVH* bvh, unsigned int current_node, unsigned int current_tri) {
	std::cout << " current node " << (bvh->cf_bvh[current_node].u.leaf.count) << " " << (bvh->cf_bvh[current_node].u.leaf.count & 0x7fffffff) << std::endl;
	if ((bvh->cf_bvh[current_node].u.leaf.count & 0x80000000)) {
		vec3 colour = vec3(static_cast <float> (rand()) / static_cast <float> (RAND_MAX), static_cast <float> (rand()) / static_cast <float> (RAND_MAX), static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		unsigned int count = (bvh->cf_bvh[current_node].u.leaf.count & 0x7fffffff);
		unsigned int first = bvh->cf_bvh[current_node].u.leaf.index_first_tri;
		for (int i = first; i < first + count; i++) {
			t_mats[bvh->triangle_indices[i]].colour = colour;

		}
		return;
	}

	ColourBVH2(bvh, bvh->cf_bvh[current_node].u.inner.index_left, current_tri);
	ColourBVH2(bvh, bvh->cf_bvh[current_node].u.inner.index_right, current_tri);

}

void Scene::UpdateMatsGPU() {
	cudaAssert(Memcpy(this->t_mats_gpu, this->t_mats.data(), this->t_mats.size() * sizeof(Material), cudaMemcpyHostToDevice));
}
