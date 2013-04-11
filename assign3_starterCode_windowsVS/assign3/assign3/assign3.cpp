/*
CSCI 480
Assignment 3 Raytracer

Name: Tim Righettini
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

// Include other libraries
#include <iostream>
#include <sstream> 
#include <vector>
#include <stdio.h>
#include <string>
#include <cmath>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

#define PI 3.14159265 // Needed for the trig calculations involved in this assignment

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

/*My Variables START*/

// Struct used to hold a point -- will be used for both vectors and points since they can be represented the same way
struct point {
	double x;
	double y;
	double z;
};

// Will hold the componets of a ray
struct ray {
	point origin; // Where the ray comes from
	point vectorDirection; // Will hold the normalized direction of the ray
	double t; // The t multiple of the ray -- used to check for collisions
};

double aspectRatio; // Will hold the aspect ratio as a decimal value -- used for the calculation of the image plane

// For the four corner points of the image plane
point cornerPoints[4];

// For holding all of the rays
ray rays[WIDTH][HEIGHT];

double pixelTraversalValueX; // Will hold how far across the image plane is needed to traverse from one pixel center to the next horizonally
double pixelTraversalValueY; // Will hold how far across the image plane is needed to traverse from one pixel center to the next vertically
// Note that half of each of these values accounts for the pixel center offset

point cameraOrigin; // Will hold the camera origin, which in the base program will be (0,0,0)

/*My Variables END*/

/*My Math Functions BEGIN*/

point getCrossProduct(point a, point b) {
	point c; // The result of the cross product

	c.x = (a.y * b.z) - (a.z * b.y);
	c.y = (a.z * b.x) - (a.x * b.z);
	c.z = (a.x * b.y) - (a.y * b.x);

	return c;
}

point getUnitVector(point a) { // Unitize a vector
	// Get the magnitude
	float aDist = pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2);

	if (aDist == 0) {
		return a; // Cannot divide by zero
	}

	aDist = sqrt(aDist);

	// Normalize the vector
	a.x /= aDist;
	a.y /= aDist;
	a.z /= aDist;
		
	return a;
}

double getDotProduct(point a, point b) { // Return the dot product between two vectors
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

double getDistance(point a, point b) { // Will use the distance formula to return a distance
	return sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y), 2));
}

/*My Math Functions END*/

/* RAY TRACING FUNCTIONS START*/

// Declare function prototypes:
void doStepOne(); // This will complete all of the raycasting, which includes the following three functions
void calculateFourCornerPoints(); // Will calculate the four corner points of the image plane
void calculatePixelTraversalIntervals(); // Will calculate how much distance needs to be traversed per pixel with ray creation
void calculateRays(); // When all of the appropriate values are found, actually calculate the rays

void doStepOne() {
	std::cout << "-----Completing Step One-----"<< std::endl;
	calculateFourCornerPoints();
	calculatePixelTraversalIntervals();
	calculateRays();
}

void calculateFourCornerPoints() {
	// Set the origin for the camera
	cameraOrigin.x = 0; cameraOrigin.y = 0;	cameraOrigin.z = 0;

	// Set the aspect ratio value
	aspectRatio = (double)WIDTH/HEIGHT;

	// Convert FOV/2 to radians for tan() 
	double fovRadiansDivTwo = fov * (PI/180) / 2;

	// Now calculate the base x, y, and z values for the image plane coordinates
	point baseImagePlanePoint;
	baseImagePlanePoint.x = aspectRatio * tan(fovRadiansDivTwo);
	baseImagePlanePoint.y = tan(fovRadiansDivTwo);
	baseImagePlanePoint.z = -1;

	// Print base image plane corner point
	std::cout << "Base Image Plane Corner x: " << baseImagePlanePoint.x << std::endl;
	std::cout << "Base Image Plane Corner y: " << baseImagePlanePoint.y << std::endl;
	std::cout << "Base Image Plane Corner z: " << baseImagePlanePoint.z << std::endl;

	// Based upon the base values above, create the four corner points
	cornerPoints[0] = baseImagePlanePoint;
	cornerPoints[1] = baseImagePlanePoint;
	cornerPoints[2] = baseImagePlanePoint;
	cornerPoints[3] = baseImagePlanePoint;

	// Now make adjustments, depending on what the point is 

	// Top Left Point
	cornerPoints[0].x *= -1;
	cornerPoints[0].y *= 1;

	// Top Right Point
	cornerPoints[1].x *= 1;
	cornerPoints[1].y *= 1;

	// Bottom Left Point
	cornerPoints[2].x *= -1;
	cornerPoints[2].y *= -1;

	// Bottom Right Point
	cornerPoints[3].x *= 1;
	cornerPoints[3].y *= -1;

	// Now print out the four corner points
	std::cout << "-----Resulting Corner Points-----"<< std::endl;
	for (int i = 0; i < 4; i++) {
		std::cout << "Point " << i << ": Image Plane Corner x: " << cornerPoints[i].x << std::endl;
		std::cout << "Point " << i << ": Image Plane Corner y: " << cornerPoints[i].y << std::endl;
		std::cout << "Point " << i << ": Image Plane Corner z: " << cornerPoints[i].z << std::endl;
	}
}

void calculatePixelTraversalIntervals() {

}

void calculateRays() {

}


/* RAY TRACING FUNCTIONS END */

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      plot_pixel(x,y,x%256,y%256,(x+y)%256);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
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
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
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
      else if(stricmp(type,"light")==0)
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

  // Do the raytracing calculations
  doStepOne(); // Uniformly send out rays from one location
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
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
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
