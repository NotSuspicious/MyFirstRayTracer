# Created by William Zhao
Completed as an assignment in CSCI 420 at USC.
Thank you Jernej Barbic, USC.

# Samples
Samples 001.jpg and 002.jpg show features 7 & 8.
![Shapes](https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/001.jpg)
![Spheres](https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/002.jpg)

Samples 002.jpg and 004.jpg have shadow acne and no AA because it's slow to render
![Shapes](https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/003.jpg)
![Spheres](https://raw.githubusercontent.com/NotSuspicious/MyFirstRayTracer/refs/heads/main/RayTracer/samples/004.jpg)

# Features
1. Sphere and Triangle support
2. Phong Shading
3. Image saving
4. Soft Shadows with Area Lights
- To reduce shadow acne, increase LIGHT_SAMPLES in hw3.cpp
5. Super Sampling Anti-Aliasing
- To adjust sample size, change SUPER_SAMPLING to 2 in hw3.cpp

# Instructions
1. Compile using the Makefile in ```MyFirstRayTracer/RayTracer```
2. Usage ```./hw3 <input scenefile> [output jpegname]```
   - Scenes are in MyFirstRayTracer/RayTracer/scenes
3. Render test1.scene and spheres.scene first, to get a sense of render time.
   - If nothing shows up for a while, try reducing LIGHT_SAMPLES to 1 and SUPER_SAMPLING to 1
