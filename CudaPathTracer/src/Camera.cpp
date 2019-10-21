#pragma once
#include <precomp.h>

using namespace glm;

Camera::Camera(){
	this->c_matrix = mat4(1.0f);
	this->fov = 90;
	this->c_position = vec3(0, 0, 0);
	this->translation = glm::vec3(0);
}

Camera::Camera(vec4 p, mat4 m, int f){
	this->c_position = p;
	this->c_matrix = m;
	this->fov = f;
}