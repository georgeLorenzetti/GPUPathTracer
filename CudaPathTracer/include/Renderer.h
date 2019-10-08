#pragma once
#include "precomp.h"

class Renderer {
	public:
		Renderer(int w, int h);
		~Renderer();
		cudaError set_size(const int width, const int height);
		void draw();

		int width;
		int height;
		cudaGraphicsResource* graphics_resource = nullptr;
		cudaArray* cuda_array = nullptr;
		GLuint render_buffer;
		GLuint frame_buffer;

};

