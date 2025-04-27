Assignment #3: Ray tracing

FULL NAME: William Zhao


MANDATORY FEATURES
------------------

Feature:                                 
-------------------------------------    
1) Ray tracing triangles

2) Ray tracing sphere

3) Triangle Phong Shading

4) Sphere Phong Shading

5) Shadows rays

6) Still images
   
7) Soft Shadows using Area Lights
- To reduce shadow acne, increase LIGHT_SAMPLES in hw3.cpp

8) Super Sampling Anti Aliasing
- To adjust sample size, change SUPER_SAMPLING to 2 in hw3.cpp

Samples
-------------------------------------
Samples 001.jpg and 002.jpg show features 7 & 8.

Samples 002.jpg and 004.jpg have shadow acne and no AA because it's slow to render

Instructions
------------------
Render Test1.scene and Spheres.scene first, to get a sense of render time.
If nothing shows up for a while, try reducing LIGHT_SAMPLES to 1 and SUPER_SAMPLING to 1
