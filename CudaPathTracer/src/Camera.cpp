#pragma once
#include <precomp.h>

using namespace glm;

Camera::Camera(){
	this->c_position = vec4(0, 0, 0, 1);
	this->c_matrix = mat4(1.0f);
	this->fov = 90;
}

Camera::Camera(vec4 p, mat4 m, int f){
	this->c_position = p;
	this->c_matrix = m;
	this->fov = f;
}