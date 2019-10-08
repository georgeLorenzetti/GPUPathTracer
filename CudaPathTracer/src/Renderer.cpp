#pragma once
#include <precomp.h>

Renderer::Renderer(int _width, int _height){
	glCreateRenderbuffers(1, &render_buffer);
	glCreateFramebuffers(1, &frame_buffer);

	glNamedFramebufferRenderbuffer(frame_buffer, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, render_buffer);

	set_size(_width, _height);
}
Renderer::~Renderer(){
	cudaError cuda_err;

	// unregister CUDA resources
	if (graphics_resource != nullptr)
		cuda_err = cudaAssert(GraphicsUnregisterResource(graphics_resource));

	glDeleteRenderbuffers(1, &render_buffer);
	glDeleteFramebuffers(1, &frame_buffer);
	
}

cudaError Renderer::set_size(int _width, int _height){
	cudaError cuda_err = cudaSuccess;

	width = _width;
	height = _height;

	if (graphics_resource != nullptr)
		cuda_err = cudaAssert(GraphicsUnregisterResource(graphics_resource));

	glNamedRenderbufferStorage(render_buffer, GL_RGBA32F, width, height);
	cuda_err = cudaAssert(GraphicsGLRegisterImage(&graphics_resource, render_buffer, GL_RENDERBUFFER, cudaGraphicsRegisterFlagsSurfaceLoadStore | cudaGraphicsRegisterFlagsWriteDiscard));
	cuda_err = cudaAssert(GraphicsMapResources(1, &graphics_resource, 0));
	cuda_err = cudaAssert(GraphicsSubResourceGetMappedArray(&cuda_array, graphics_resource, 0, 0));
	cuda_err = cudaAssert(GraphicsUnmapResources(1, &graphics_resource, 0));

	return cuda_err;
}

void Renderer::draw(){
	glBlitNamedFramebuffer(frame_buffer, 0, 0, 0, width, height, 0, height, width, 0, GL_COLOR_BUFFER_BIT, GL_NEAREST);
}