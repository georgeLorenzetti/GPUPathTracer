#pragma once
#include <precomp.h>

class Camera{
public:
	Camera();
	Camera(glm::vec4 p, glm::mat4 m, int f);

	glm::vec4 c_position;
	glm::mat4 c_matrix;
	int fov;

	float translation_speed = 2.0f;
	float rotation_theta = 1.0f;
};

