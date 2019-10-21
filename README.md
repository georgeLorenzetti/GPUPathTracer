# GPUPathTracer

A Wavefront path tracer written in C++/CUDA.

Current Features:
 - Next Event Estimation
 - Lambertian materials

Features to add:
 - Specular/Refractive/Frensel materials
 - Cosine weighted importance sampling
 - Russian Roulette termination
 - BVH Acceleration structue/.obj loading
 - Texturing
 - Change it so all paths are finished before moving to next frame
 - Probs extra stuff that I haven't though of yet
 
 ![alt text](https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_IPgWP2Y2Lt.png)
 
 Bugs (to remind myself what I need to fix)
 - Different window sizes act strange
