# GPUPathTracer

A real-time graphics renderer utilising the Wavefront Path Tracer algorithm written in C++/CUDA.

Current Features:
 - Next event estimation
 - Cosine weighted hemisphere sampling
 - Snell/Fresnel Dialectric, Lambertian & purely specular materials
 - (Quad) MBVH acceleration structure with 32 byte nodes
 - Skybox

To Do:
 - More BSDFs
 - Russian Roulette termination
 - Object Texturing
 
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_AYH7Zq44ZS.png"></img>
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_MPbI74MEJx.png"></img>
 <img src="https://github.com/georgeLorenzetti/GPUPathTracer/blob/master/CudaPathTracer/screenies/CudaPathTracer_6k8vMpKmfA.png"></img>
Dragon scene (~800k triangles) renders in real-time at a constant >60fps with an RTX2060.
