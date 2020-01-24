# GPUPathTracer

A Wavefront path tracer written in C++/CUDA.

Current Features:
 - Next event estimation
 - Lambertian materials
 - Cosine weighted hemisphere sampling
 - (Quad) MBVH acceleration structure with 32 byte nodes.

To Do:
 - More BRDFs
 - Russian Roulette termination
 - Texturing
 - Skybox
 - Probs extra stuff that I haven't though of yet
 
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_IPgWP2Y2Lt.png"></img>
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_Ehj1jNW2Ks.png"></img>
 
 Bugs (to remind myself what I need to fix)
 - Different window sizes act strange
