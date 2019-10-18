#pragma once

//define values
#define NOMINMAX
#define PI 3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f
#define SCRWIDTH 512	
#define SCRHEIGHT 512
#define INVPI 1 / PI
#define DOUBLEPI 2 * PI
#define MAXBOUNCE 5
#define BUFFERSIZE (SCRWIDTH*SCRHEIGHT)
#define MAXDISTANCE 1e20f;


//variables

//Dependancy headers
#include <GL/glew.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <GLFW/glfw3.h>
#include <glad/glad.h>
// C++ headers
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

// Namespaced C headers:
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>

//CUDA headers
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <device_launch_parameters.h>
#include <CudaAssert.h>

#include <Game.h>
#include <Material.cuh>
#include <Scene.h>
#include <Kernel.cuh>
#include <Renderer.h>
#include <Camera.h>
#include <PathTracer.h>

#include <Definitions.h>
