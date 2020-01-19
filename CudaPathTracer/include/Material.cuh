#pragma once
#include <precomp.h>

/* Material IDs
0 = Background
1 = Light
2 = Lambertian
*/
class Material{
public:
	__host__ __device__ Material(){};
	__host__ __device__ Material(int _type, glm::vec3 _colour){
		this->type = _type;
		this->colour = _colour;
		this->specularity = 0.7f;
	};
	__host__ __device__ Material(int _type, glm::vec3 _colour, float _specularity) {
		this->type = _type;
		this->colour = _colour;
		this->specularity = _specularity;
	}

	int type;
	glm::vec3 colour;
	float specularity;
};

