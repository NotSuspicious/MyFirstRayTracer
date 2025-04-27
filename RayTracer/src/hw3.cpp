/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: William Zhao
 * *************************
*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Math.h"
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100
#define SMALL_VALUE 1e-4
#define LIGHT_SAMPLES 64
#define SUPER_SAMPLING 2

char * filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2
#define MAX_DEPTH 1

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera, in degrees.
#define fov 60.0
#define BACKGROUND_COLOR Color(0.0f, 0.0f, 0.0f)

// Buffer to store the image when saving it to a JPEG.
unsigned char buffer[HEIGHT][WIDTH][3];

unsigned char colors[HEIGHT*SUPER_SAMPLING][WIDTH*SUPER_SAMPLING][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Color {
  float r,g,b,a;
  Color() : r(0), g(0), b(0), a(0) {}
  Color(float _r, float _g, float _b, float _a=1.0f)
    : r(_r), g(_g), b(_b), a(_a) {}
  Color operator+(const Color& o) const {
    return Color(r+o.r, g+o.g, b+o.b, a+o.a);
  }
  Color operator*(float o) const {
    return Color(r*o, g*o, b*o, a);
  }
};

struct Triangle
{
  Vertex v[3];
  Vector3 planeNormal;
  void CalculatePlaneNormal(){
    Vector3 A (v[0].position[0], v[0].position[1], v[0].position[2]);
    Vector3 B(v[1].position[0], v[1].position[1], v[1].position[2]);
    Vector3 C(v[2].position[0], v[2].position[1], v[2].position[2]);
    Vector3 AB = B - A;
    Vector3 AC = C - A;
    planeNormal = Vector3::Cross(AB, AC);
    planeNormal.Normalize();
  }

  Vector3 GetInterpolatedNormal(Vector3 p) {
    Vector3 bary = CalculateBarycentricCoordinates(*this, p);
    Vector3 normal;
    normal = normal + bary.x * Vector3(v[0].normal[0], v[0].normal[1], v[0].normal[2]);
    normal = normal + bary.y * Vector3(v[1].normal[0], v[1].normal[1], v[1].normal[2]);
    normal = normal + bary.z * Vector3(v[2].normal[0], v[2].normal[1], v[2].normal[2]);
    normal.Normalize();
    return normal;
  }

  Vector3 CalculateBarycentricCoordinates(Triangle triangle, Vector3 p){
    Vector3 barycentric;
    float ax, ay, bx, by, cx, cy, px, py;
    ProjectTriangleAndPointTo2D(triangle, p, ax, ay, bx, by, cx, cy, px, py);

    float areaABC = GetArea2D(ax, ay, bx, by, cx, cy);
    if (fabs(areaABC) < SMALL_VALUE) {
        barycentric.x = barycentric.y = barycentric.z = -1.0f;
        return barycentric;
    }

    float areaPBC = GetArea2D(px, py, bx, by, cx, cy);
    float areaAPC = GetArea2D(ax, ay, px, py, cx, cy);
    float areaABP = GetArea2D(ax, ay, bx, by, px, py);

    barycentric.x = areaPBC / areaABC;
    barycentric.y = areaAPC / areaABC;
    barycentric.z = areaABP / areaABC;

    return barycentric;
  }

  static float GetArea2D(float ax, float ay, float bx, float by, float cx, float cy){
    return (0.5f) * ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay));
  }

  void ProjectPointTo2D(const Vector3& normal, const Vector3& p, float& u, float& v) {
    float nx = fabs(normal.x);
    float ny = fabs(normal.y);
    float nz = fabs(normal.z);
    if (nx > ny && nx > nz) {
        // Project to YZ
        u = p.y;
        v = p.z;
    } else if (ny > nz) {
        // Project to XZ
        u = p.x;
        v = p.z;
    } else {
        // Project to XY
        u = p.x;
        v = p.y;
    }
  }

  void ProjectTriangleAndPointTo2D(const Triangle& tri, const Vector3& p,
    float& ax, float& ay,
    float& bx, float& by,
    float& cx, float& cy,
    float& px, float& py)
  {
    Vector3 normal = tri.planeNormal;
    ProjectPointTo2D(normal, Vector3(tri.v[0].position[0], tri.v[0].position[1], tri.v[0].position[2]), ax, ay);
    ProjectPointTo2D(normal, Vector3(tri.v[1].position[0], tri.v[1].position[1], tri.v[1].position[2]), bx, by);
    ProjectPointTo2D(normal, Vector3(tri.v[2].position[0], tri.v[2].position[1], tri.v[2].position[2]), cx, cy);
    ProjectPointTo2D(normal, p, px, py);
  }
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
  double radius = 0.5f;
  Vector3 GetRandomPointOnLight() {
    Vector3 lightPos(position[0], position[1], position[2]);
    Vector3 randomPoint = lightPos + Vector3::RandomInUnitSphere() * radius;
    return randomPoint;
  }
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);


Color CastShadowRayFromTriangle(Vector3 point, Triangle source, Ray view, float view_t, int depth = 0);
Color CastShadowRayFromSphere(Vector3 point, Sphere source, Ray view, float view_t, int depth = 0);

Ray cameraRays[HEIGHT * SUPER_SAMPLING][WIDTH * SUPER_SAMPLING];

void GenerateCameraRays() {
  float aspect = (float)WIDTH / (float)HEIGHT;
  float scale = tan(fov * 0.5f * M_PI / 180.0f);
  for (int y = 0; y < HEIGHT * SUPER_SAMPLING; y++) {
    for (int x = 0; x < WIDTH * SUPER_SAMPLING; x++) {
      float px = (2.0f * ((x + 0.5f) / (float)(WIDTH * SUPER_SAMPLING)) - 1.0f) * aspect * scale;
      float py = (1.0f - 2.0f * ((y + 0.5f) / (float)(HEIGHT * SUPER_SAMPLING))) * scale;
      Vector3 dir(px, py, -1.0f);
      cameraRays[HEIGHT * SUPER_SAMPLING - y - 1][x] = Ray(Vector3(0, 0, 0), dir);
    }
  }
}

bool HasHitSphere(Sphere sphere, Ray ray, float& t){
  Vector3 d = ray.d;
  Vector3 p = ray.p;
  Vector3 c(sphere.position[0], sphere.position[1], sphere.position[2]);
  Vector3 oc = p - c;

  float a = Vector3::Dot(d, d); 
  float b = 2.0f * Vector3::Dot(d, oc);
  float c_term = Vector3::Dot(oc, oc) - sphere.radius * sphere.radius;

  float discriminant = b * b - 4 * a * c_term;
  if (discriminant <= 0.0f) {
      return false;
  }

  float sqrt_disc = sqrt(discriminant);
  float t0 = (-b - sqrt_disc) / (2.0f * a);
  float t1 = (-b + sqrt_disc) / (2.0f * a);

  if (t0 > SMALL_VALUE) {
      t = t0;
      return true;
  } else if (t1 > SMALL_VALUE) {
      t = t1;
      return true;
  }
  return false;
}

bool HasHitTriangle(Triangle triangle, Ray ray, float& t){
  float ray_triplanenorm_dot = Vector3::Dot(triangle.planeNormal, ray.d);
  if (fabs(ray_triplanenorm_dot) < SMALL_VALUE){
      return false; //ray parallel to triangle
  }

  //intersection with triangle plane
  Vector3 v0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  float d = Vector3::Dot(triangle.planeNormal, v0);
  float t_temp = (d - Vector3::Dot(triangle.planeNormal, ray.p)) / ray_triplanenorm_dot;

  if (t_temp < SMALL_VALUE) {
    return false; // Intersection behind ray origin
  }

  Vector3 p = ray.p + ray.d * t_temp;
  Vector3 bary = triangle.CalculateBarycentricCoordinates(triangle, p);

  // Check if point is inside triangle
  if (bary.x >= 0.0f && bary.y >= 0.0f && bary.z >= 0.0f &&
    bary.x <= 1.0f && bary.y <= 1.0f && bary.z <= 1.0f &&
    fabs(bary.x + bary.y + bary.z - 1.0f) < SMALL_VALUE) {
    t = t_temp;
    return true;
  }

  return false;
}

bool CheckSphereIntersections(Ray ray, Sphere& hitSphere, float& t){
  bool hasSphereIntersect = false;
  float minDistSqr = MAXFLOAT;
  t = -1.0f;
  for (int i = 0 ; i < num_spheres ; i++){
    if (spheres[i].radius < SMALL_VALUE){
      continue;
    }
    float rayT;
    if (HasHitSphere(spheres[i], ray, rayT)){
      float distSqr = ray.at(rayT).LengthSquared();
      if (distSqr < minDistSqr){
        hasSphereIntersect = true;
        minDistSqr = distSqr;
        t = rayT;
        hitSphere = spheres[i];
      }
    }
  }
  return hasSphereIntersect;
}

bool CheckTriangleIntersections(Ray ray, Triangle& hitTriangle, float& t)
{
  bool hasTriangleIntersect = false;
    float minDistSqr = MAXFLOAT;
    t = -1.0f;
    for (int i = 0; i < num_triangles; i++) {
        float rayT;
        if (HasHitTriangle(triangles[i], ray, rayT)) {
            Vector3 hitPoint = ray.p + ray.d * rayT;
            float distSqr = (hitPoint - ray.p).LengthSquared();
            if (distSqr < minDistSqr) {
                hasTriangleIntersect = true;
                minDistSqr = distSqr;
                t = rayT;
                hitTriangle = triangles[i];
            }
        }
    }
  return hasTriangleIntersect;
}

Color CalculateSphereColor(Sphere sphere, Light light, Vector3 l, Vector3 n, Vector3 v)
{
  Vector3 r = Vector3::Reflect(-1 * l, n);
  float kd[3] = {sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]};
  float ks[3] = {sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]};
  float shininess = sphere.shininess;

  float l_dot_n = Vector3::Dot(l, n);
  if (l_dot_n < 0){ l_dot_n = 0; }
  float r_dot_v = Vector3::Dot(r, v);
  if (r_dot_v < 0){  r_dot_v = 0;  }
  
  Color sphereColor;

  sphereColor.r = light.color[0] * (kd[0] * l_dot_n + ks[0] * pow(r_dot_v, shininess));
  sphereColor.g = light.color[1] * (kd[1] * l_dot_n + ks[1] * pow(r_dot_v, shininess));
  sphereColor.b = light.color[2] * (kd[2] * l_dot_n + ks[2] * pow(r_dot_v, shininess));
  sphereColor.a = 1;

  return sphereColor;
}

Color CalculateTriangleColor(Triangle triangle, Light light, Vector3 point, Vector3 l, Vector3 n, Vector3 v)
{
  Vector3 r = Vector3::Reflect(-1 *l, n);
  Vector3 bary = triangle.CalculateBarycentricCoordinates(triangle, point);
  float kd[3] = {bary.x * triangle.v[0].color_diffuse[0] + bary.y * triangle.v[1].color_diffuse[0] + bary.z * triangle.v[2].color_diffuse[0],
              bary.x * triangle.v[0].color_diffuse[1] + bary.y * triangle.v[1].color_diffuse[1] + bary.z * triangle.v[2].color_diffuse[1],
              bary.x * triangle.v[0].color_diffuse[2] + bary.y * triangle.v[1].color_diffuse[2] + bary.z * triangle.v[2].color_diffuse[2]}; 
  float ks[3] = {bary.x * triangle.v[0].color_specular[0] + bary.y * triangle.v[1].color_specular[0] + bary.z * triangle.v[2].color_specular[0],
              bary.x * triangle.v[0].color_specular[1] + bary.y * triangle.v[1].color_specular[1] + bary.z * triangle.v[2].color_specular[1],
              bary.x * triangle.v[0].color_specular[2] + bary.y * triangle.v[1].color_specular[2] + bary.z * triangle.v[2].color_specular[2]};
  float shininess = bary.x * triangle.v[0].shininess + bary.y * triangle.v[1].shininess + bary.z * triangle.v[2].shininess;

  float l_dot_n = Vector3::Dot(l, n);
  if (l_dot_n < 0){ l_dot_n = 0; }
  float r_dot_v = Vector3::Dot(r, v);
  if (r_dot_v < 0){  r_dot_v = 0;  }

  Color triangleColor;

  triangleColor.r = light.color[0] * (kd[0] * l_dot_n + ks[0] * pow(r_dot_v, shininess));
  triangleColor.g = light.color[1] * (kd[1] * l_dot_n + ks[1] * pow(r_dot_v, shininess));
  triangleColor.b = light.color[2] * (kd[2] * l_dot_n + ks[2] * pow(r_dot_v, shininess));
  triangleColor.a = 1;

  return triangleColor;
}

Color CastShadowRayFromTriangle(Vector3 point, Triangle source, Ray view, float view_t, int depth)
{
  if (depth >= MAX_DEPTH){
    return Color();
  }
  Color color;
  for (int i = 0; i < num_lights; i++) {
    Color total;
    for (int s=0; s<LIGHT_SAMPLES; s++){
      Vector3 lightPos = lights[i].GetRandomPointOnLight();
      Vector3 L = lightPos - point; float maxT = L.Length(); L.Normalize();
      Vector3 N = source.GetInterpolatedNormal(point);
      Ray ray(point + SMALL_VALUE * N, L);

      float sphere_t, tri_t;
      Sphere hitSphere;
      Triangle hitTriangle;
      bool sphere_intersect = CheckSphereIntersections(ray, hitSphere, sphere_t);
      bool tri_intersect = CheckTriangleIntersections(ray, hitTriangle, tri_t);

      if ((sphere_intersect && sphere_t > SMALL_VALUE && sphere_t < maxT) ||
          (tri_intersect && tri_t > SMALL_VALUE && tri_t < maxT)) {
        continue;
      }
      else{
        Vector3 l = ray.d; l.Normalize();
        Vector3 n = source.GetInterpolatedNormal(point);
        Vector3 v = -1.0 * view.d; v.Normalize();
        Color newCol = CalculateTriangleColor(source, lights[i], point, l, n, v);
        total = total + newCol;
      }
    }
    Color softshadowColor = total * (1.0f / LIGHT_SAMPLES);
    color = color + softshadowColor;
  }
  // ambient
  Color amb(ambient_light[0], ambient_light[1], ambient_light[2]);
  color = color + amb;
  return color;
}

Color CastShadowRayFromSphere(Vector3 point, Sphere source, Ray view, float view_t, int depth)
{
  if (depth >= MAX_DEPTH){
    return Color();
  }
  Color color;
  for (int i = 0; i < num_lights; i++) {
    Color total;
    for(int s=0; s<LIGHT_SAMPLES; s++){
      Vector3 lightPos = lights[i].GetRandomPointOnLight();
      Vector3 L = lightPos - point; float maxT = L.Length(); L.Normalize();
      Vector3 N = (point - Vector3(source.position[0], source.position[1], source.position[2])) / source.radius;
      Ray ray(point + SMALL_VALUE * N, L);

      float sphere_t, tri_t;
      Sphere hitSphere;
      Triangle hitTriangle;
      bool sphere_intersect = CheckSphereIntersections(ray, hitSphere, sphere_t);
      bool tri_intersect = CheckTriangleIntersections(ray, hitTriangle, tri_t);

      if ((sphere_intersect && sphere_t > 0 && sphere_t < maxT) ||
          (tri_intersect && tri_t > 0 && tri_t < maxT)) {
        continue;
      }
      else{
        Vector3 l = ray.d; l.Normalize();
        Vector3 n = (point - Vector3(source.position[0], source.position[1], source.position[2])) / source.radius;
        Vector3 v = -1.0 * view.d; v.Normalize();
        Color newCol = CalculateSphereColor(source, lights[i], l, n, v);
        total = total + newCol;
      }
    }
    Color softshadowColor = total * (1.0f / LIGHT_SAMPLES);
    color = color + softshadowColor;
  }
  // ambient
  Color amb(ambient_light[0], ambient_light[1], ambient_light[2]);
  color = color + amb;
  return color;
}

Color CastCameraRay(Ray ray, int depth = 0)
{
  if (depth >= MAX_DEPTH){
    return Color();
  }
  Color color;
  float sphere_t, tri_t;
  Sphere hitSphere;
  Triangle hitTriangle;
  bool sphere_intersect = CheckSphereIntersections(ray, hitSphere, sphere_t);
  bool tri_intersect = CheckTriangleIntersections(ray, hitTriangle, tri_t);

  if (sphere_intersect && tri_intersect) {
    if (sphere_t < tri_t) {
      Vector3 hitPosition = ray.at(sphere_t);
      return CastShadowRayFromSphere(hitPosition, hitSphere, ray, sphere_t);
    } else {
      Vector3 hitPosition = ray.at(tri_t);
      return CastShadowRayFromTriangle(hitPosition, hitTriangle, ray, tri_t);
    }
  }
  else if (sphere_intersect) {
    Vector3 hitPosition = ray.at(sphere_t);
    // return Color(hitSphere.color_diffuse[0], hitSphere.color_diffuse[1], hitSphere.color_diffuse[2]);
    return CastShadowRayFromSphere(hitPosition, hitSphere, ray, sphere_t);
  } else if (tri_intersect) {
    Vector3 hitPosition = ray.at(tri_t);
    // return Color(hitTriangle.v[0].color_diffuse[0], hitTriangle.v[0].color_diffuse[1], hitTriangle.v[0].color_diffuse[2]);
    return CastShadowRayFromTriangle(hitPosition, hitTriangle, ray, tri_t);
  } else {
    return BACKGROUND_COLOR;
  }
}

bool IsInColors(int x, int y)
{
  if (x < 0 || x >= WIDTH * SUPER_SAMPLING || y < 0 || y >= HEIGHT * SUPER_SAMPLING) {
    return false;
  }
  return true;
}

void draw_scene()
{
  GenerateCameraRays();
  for (int y = 0; y < HEIGHT * SUPER_SAMPLING; y++) {
    for (int x = 0; x < WIDTH * SUPER_SAMPLING; x++) {
      Ray ray = cameraRays[y][x];
      Color color = CastCameraRay(ray);
      unsigned char r = (unsigned char)(color.r * 255);
      unsigned char g = (unsigned char)(color.g * 255);
      unsigned char b = (unsigned char)(color.b * 255);
      colors[y][x][0] = r;
      colors[y][x][1] = g;
      colors[y][x][2] = b;
    }
  }
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    // Do not worry about this usage of OpenGL. This is here just so that we can draw the pixels to the screen,
    // after their R,G,B colors were determined by the ray tracer.
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      Color superSampledColor;
      for (int i = 0; i < SUPER_SAMPLING; i++) {
        for (int j = 0; j < SUPER_SAMPLING; j++) {
          if (IsInColors(x*SUPER_SAMPLING+i, y*SUPER_SAMPLING+j)) {
            unsigned char g = colors[y*SUPER_SAMPLING+i][x*SUPER_SAMPLING+j][1];
            unsigned char r = colors[y*SUPER_SAMPLING+i][x*SUPER_SAMPLING+j][0];
            unsigned char b = colors[y*SUPER_SAMPLING+i][x*SUPER_SAMPLING+j][2];
            superSampledColor = superSampledColor + Color(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
          }
          
        }
      }
      superSampledColor = superSampledColor * (1.0f / pow(SUPER_SAMPLING, 2));
      unsigned char r = (unsigned char)(superSampledColor.r * 255);
      unsigned char g = (unsigned char)(superSampledColor.g * 255);
      unsigned char b = (unsigned char)(superSampledColor.b * 255);
      plot_pixel(x, y, r, g, b);
    }

    glEnd();
    glFlush();
  }

  printf("Ray tracing completed.\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      t.CalculatePlaneNormal();
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

