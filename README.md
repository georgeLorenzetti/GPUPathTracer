# GPUPathTracer

A Wavefront path tracer written in C++/CUDA.

Current Features:
 - Next event estimation
 - Lambertian materials
 - Cosine weighted hemisphere sampling
 - BVH Acceleration structue/.obj loading

To Do:
 - Specular/Refractive/Frensel materials
 - Russian Roulette termination
 - Texturing
 - Skybox
 - Probs extra stuff that I haven't though of yet
 
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_IPgWP2Y2Lt.png" height="300" width="300"></img>
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_9OfNaPH2vr.png" height="300" width="300"></img>
 
 Bugs (to remind myself what I need to fix)
 - Different window sizes act strange
