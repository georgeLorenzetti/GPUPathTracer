#pragma once
#include <precomp.h>

class Camera{
public:
	Camera();
	Camera(glm::vec3 p, glm::mat4 m, int f);

	glm::vec3 GetPosition(){ return c_position; };
	glm::mat4 GetMatrix(){ return c_matrix; };
	int GetFOV(){ return fov; };

	void SetPosition(glm::vec3 p){ c_position = p; };
	void SetMatrix(glm::mat4 m){ c_matrix = m; };
	void SetFOV(int i){ fov = i; };
private:
	glm::vec3 c_position;
	glm::mat4 c_matrix;
	int fov;
};

