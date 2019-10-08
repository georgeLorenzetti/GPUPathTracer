#pragma once
#include <precomp.h>

using namespace glm;

Camera::Camera(){
	this->c_position = vec3(0, 0, 0);
	this->c_matrix = mat4(1.0f);
	this->fov = 90;
}

Camera::Camera(vec3 p, mat4 m, int f){
	this->c_position = p;
	this->c_matrix = m;
	this->fov = f;
}