#pragma once
#include <precomp.h>

class Camera{
public:
	Camera();
	Camera(glm::vec4 p, glm::mat4 m, int f);

	glm::vec3 c_position;
	glm::mat4 c_matrix;
	glm::vec3 translation;
	int fov;
	float translation_speed = 2.5f;
	float rotation_theta = 0.5f;
};

