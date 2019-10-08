#pragma once
#include "precomp.h"

static void glfw_error_callback(int error, const char* description){
	std::cout << description << "\n";
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window){
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height){
	// make sure the viewport matches the new window dimensions; note that width and 
	// height will be significantly larger than specified on retina displays.
	glViewport(0, 0, width, height);
}

int main(){
	// glfw: initialize and configure
	// ------------------------------
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

	// glfw window creation
	// --------------------
	GLFWwindow* window = glfwCreateWindow(SCRWIDTH, SCRHEIGHT, "CudaTracer", NULL, NULL);
	if (window == NULL){
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwShowWindow(window);
	glewExperimental = GL_TRUE;
	GLenum error = glewInit();

	if (error != GLEW_OK){
		std::cout << "BAD BAD GLEW NOT WORK: " << glewGetErrorString(error) << std::endl;
		exit(EXIT_FAILURE);
	}

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// glad: load all OpenGL function pointers
	// ---------------------------------------
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	//CUDA setup
	cudaError err;
	int opengl_device_id;
	unsigned int device_count;
	err = cudaAssert(GLGetDevices(&device_count, &opengl_device_id, 1u, cudaGLDeviceListAll));
	int cuda_device_id = 0;
	err = cudaAssert(SetDevice(cuda_device_id));

	const bool multi_gpu = opengl_device_id != cuda_device_id;
	struct cudaDeviceProp properties;

	err = cudaAssert(GetDeviceProperties(&properties, opengl_device_id));
	printf("GL   : %-24s (%2d)\n", properties.name, properties.multiProcessorCount);

	err = cudaAssert(GetDeviceProperties(&properties, cuda_device_id));
	printf("CUDA : %-24s (%2d)\n", properties.name, properties.multiProcessorCount);

	//Set up cuda interop
	Renderer cuda_interop = Renderer(SCRWIDTH, SCRHEIGHT);
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);


	//Pathtracer variables
	cudaAssert(DeviceSynchronize());
	glm::vec4* frame_buffer;

	//allocate CUDA memory for variables
	cudaAssert(Malloc(&frame_buffer, SCRWIDTH * SCRHEIGHT * sizeof(glm::vec4)));
	cudaAssert(Memset(frame_buffer, 0, SCRWIDTH * SCRHEIGHT * sizeof(glm::vec4)));

	//set up path tracer object
	PathTracer path_tracer = PathTracer(properties.multiProcessorCount);
	// render loop
	// -----------
	while (!glfwWindowShouldClose(window)){
		// input
		// -----
		processInput(window);

		path_tracer.Trace(&cuda_interop, frame_buffer);


		// render
		// ------
		cuda_interop.draw();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// glfw: terminate, clearing all previously allocated GLFW resources.
	// ------------------------------------------------------------------
	glfwTerminate();
	return 0;
}