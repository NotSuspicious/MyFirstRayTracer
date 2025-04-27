# My First Ray Tracer
Completed as an assignment in CSCI 420 at USC.

Thank you Jernej Barbic, USC.

# Samples
Samples 1 and 2 show Super Sampling Anti-Aliasing (SSAA) and Soft Shadows.

Samples 3 and 4 have shadow acne and no SSAA because it's slow to render.

<p align="center">
  <img src="https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/001.jpg" width="300"/>
  <img src="https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/002.jpg" width="300"/>
  <img src="https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/003.jpg" width="300"/>
  <img src="https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/004.jpg" width="300"/>
</p>

# Features
1. Sphere and Triangle support
2. Phong Shading
3. Image saving
4. Soft Shadows with Area Lights
- To reduce shadow acne, increase ```LIGHT_SAMPLES``` in ```hw3.cpp```
5. Super Sampling Anti-Aliasing
- To adjust sample size, change ```SUPER_SAMPLING``` to 2 in ```hw3.cpp```

# Instructions
1. Compile using the Makefile in ```MyFirstRayTracer/RayTracer```
2. Usage ```./hw3 <input scenefile> [output jpegname]```
   - Scenes are in ```MyFirstRayTracer/RayTracer/scenes```
3. Render ```test1.scene``` and ```spheres.scene``` first, to get a sense of render time.
   - If nothing shows up for a while, try reducing LIGHT_SAMPLES to 1 and SUPER_SAMPLING to 1
