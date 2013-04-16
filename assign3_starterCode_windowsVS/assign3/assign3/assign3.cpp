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
#include <iomanip>
#include <cmath>
#include <time.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10000

#define PI 3.14159265 // Needed for the trig calculations involved in this assignment

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_JPEG;

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
int numRandomLights = 0; // Number of satellite lights that will be added around a main light, for the purposes of soft shadows
const int sampleNumber = 1; // How many extra rays will be casted per x and y -- 1 is normal sampling, two is double sampling, etc

// Struct used to hold a point -- will be used for both vectors and points since they can be represented the same way
struct point {
	double x;
	double y;
	double z;
};

// Struct that will hold the components of a color from 0 to 255
struct color {
	unsigned char red; 
	unsigned char green;
	unsigned char blue;
}; 

// Struct that will hold the components of a color from 0 to 1 double
struct colorFloat {
	double red; 
	double green;
	double blue;
}; 

// Will hold the componets of a ray

enum collisionWithShape {TRIANGLE, SPHERE};

struct ray {
	point origin; // Where the ray comes from
	point direction; // Will hold the normalized direction of the ray
	double t; // The t multiple of the ray -- used to check for collisions
	bool isSetT; // This value will initially be set to false, because t has not been set upon instantiation
	collisionWithShape collisionShape; // What type of saphe did the ray collide with?
	point collisionNormal; // Will hold the result of a normal calculation, normalized
	color collisionColor; // Will hold the color in terms of 0 to 255
	colorFloat collisionColorFloat; // Will hold the color from the Phong calculations
	int collisionIndex; // Which shape did the ray collide with, along with the shape type aforementioned?
	point collisionPoint; // Where the collision actually occurred for this ray, will be easier to reference this way, as it will be used more than once
	point barycentricRatios; // Will hold the alpha (x), beta (y), and gamma, or charlie, (z) values determined from the barycentric calculations
};

double aspectRatio; // Will hold the aspect ratio as a decimal value -- used for the calculation of the image plane

// For the four corner points of the image plane
point cornerPoints[4];

// For holding all of the rays
ray rays[WIDTH * sampleNumber][HEIGHT * sampleNumber];
colorFloat averagedColors[WIDTH][HEIGHT]; // This array will used to hold the averaged values of the array above (if it holds supersampled data)

double pixelTraversalValueX; // Will hold how far across the image plane is needed to traverse from one pixel center to the next horizonally
double pixelTraversalValueY; // Will hold how far across the image plane is needed to traverse from one pixel center to the next vertically
// Note that half of each of these values accounts for the pixel center offset

point cameraOrigin; // Will hold the camera origin, which in the base program will be (0,0,0)

point attenuationValues; // Will be used for the calculation of softening light

/*My Variables END*/

/*My Math Functions BEGIN*/

// Function protoTypes:
point getCrossProduct(point a, point b);
point getUnitVector(point a);
double getDotProduct(point a, point b);
double getDistance(point a, point b);
double getQuadraticFormula(boolean isLowerVal, double b, double c);
double getTriangleArea(point a, point b, point c);

point getCrossProduct(point a, point b) {
	point c; // The result of the cross product

	c.x = (a.y * b.z) - (a.z * b.y);
	c.y = (a.z * b.x) - (a.x * b.z);
	c.z = (a.x * b.y) - (a.y * b.x);

	return c;
}

point getUnitVector(point a) { // Unitize a vector
	// Get the magnitude
	double aDist = getDotProduct(a, a);

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
	return sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y), 2) + + pow((b.z - a.z), 2));
}

double getQuadraticFormula(boolean isLowerVal, double b, double c) { // Will use the quadratic value to return the two points of the sphere intersection
	double t = 0;
	if (isLowerVal)
		t = (-b - ( sqrt(pow(b, 2) - (4*c)) )) / 2; // Get the lower value of the quadratic formula
	else
		t = (-b + ( sqrt(pow(b, 2) - (4*c)) )) / 2; // Get the higher value of the quadratic formula
	return t;
}

double getTriangleArea(point a, point b, point c) {
	// Compute the areas with respect to the three planes	
	double d0 = 0.5 * ( ((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y)) );
	double d1 = 0.5 * ( ((b.x - a.x) * (c.z - a.z)) - ((c.x - a.x) * (b.z - a.z)) );
	double d2 = 0.5 * ( ((b.y - a.y) * (c.z - a.z)) - ((c.y - a.y) * (b.z - a.z)) );

	// Return the area that is the largest
	if (abs(d0) > abs(d1)) {
		if (abs(d0) > abs(d2)) {
			return d0;
		}
		else {
			return d2;
		}
	}
	else {
		if (abs(d1) > abs(d2)) {
			return d1;
		}
		else {
			return d2;
		}
	}
}

/*My Math Functions END*/

/* RAY TRACING FUNCTIONS START*/

// Declare function prototypes:
void doStepOne(); // This will complete all of the raycasting, which includes the following three functions
void calculateFourCornerPoints(); // Will calculate the four corner points of the image plane
void calculatePixelTraversalIntervals(); // Will calculate how much distance needs to be traversed per pixel with ray creation
void calculateRays(); // When all of the appropriate values are found, actually calculate the rays

void doStepTwo(); // This will complete all of the collosion detection
void checkCollisionsSpheres(); // This will loop through all of the spheres and find the closest point of intersection for this ray
void checkCollisionsPolygons(); // This will loop through all of the triangles and find the closest intersection point for this ray -- may overwrite what was found in the sphere for obviousl reasons.

void doStepThree(); // This will complete all of the normal and lighting calculations, depending on what was collided with
void calculateNormals(); // This will calculate the normals for both rays that intersected with spheres and triangles -- different lighting equations will be derived depending on the "shape" parameter of a ray
void calculateColor(); // This will calculate phong lighting on every ray that collided with something, else, the background color will be used
void calculateShadowRay(ray *collisionRay, double lightCollisionPoint[], ray *shadowRay); // Will return whether the shadow ray collided with an object or not, may need a check for when the origin of a ray is inside a sphere
bool checkShadowCollisionsSpheres(ray *shadowRay, double lightCollisionPoint[]); // Will check all of the spheres to see if there is a collision with a shadow ray
bool checkShadowCollisionsTriangles(ray *shadowRay, double lightCollisionPoint[]); // Will check all of the triangles to see if there is a collision with a shadow ray
void convertColorValues(); // Will convert the color values of doubles to chars, and will do any averaging if necessary

/*STEP ONE FUNCTIONS START*/
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
	// Here is where the corner points map out to for the image plane:
	/*

	0---------------1
	|               |
	|               |
	|               |
	2---------------3

	*/
}

void calculatePixelTraversalIntervals() {
	// Find the x distance between the two x points for the image plane (using the actual distance here)
	// Using corner points 0 and 1
	pixelTraversalValueX = getDistance(cornerPoints[0], cornerPoints[1]);

	// Find the y distance between the two x points for the image plane (using the actual distance here)
	// Using corner points 0 and 2
	pixelTraversalValueY = getDistance(cornerPoints[0], cornerPoints[2]);

	// Now print these values

	std::cout << "-----Resulting X/Y Distances-----" << std::endl;
	std::cout << "xDistance: " << pixelTraversalValueX << std::endl;
	std::cout << "yDistance: " << pixelTraversalValueY << std::endl;
	// Now divide these values by the width and Height, respectively, and then the actual travseral distance per pixel should be attained

	pixelTraversalValueX = (double)pixelTraversalValueX/WIDTH/sampleNumber;
	pixelTraversalValueY = (double)pixelTraversalValueY/HEIGHT/sampleNumber;

	// Now Print these values
	std::cout << "-----Resulting X/Y Distance Intervals-----" << std::endl;
	std::cout << "pixelTraversalValueX: " << pixelTraversalValueX << std::endl;
	std::cout << "pixelTraversalValueY: " << pixelTraversalValueY << std::endl;
}

void calculateRays() {
	// Get the top left pixel center by taking each pixel traveral value and dividing it by two, make the result a point
	point pixelCenter;

	pixelCenter.x = pixelTraversalValueX/2;
	pixelCenter.y = pixelTraversalValueY/2;
	pixelCenter.z = -1; // Remember, Z is always -1 for the image plane

	// Print this value, make sure that it prints correctly
	std::cout << "-----PixelCenter X/Y/Z Values-----" << std::endl;
	std::cout << "pixelCenter.x " << pixelCenter.x << std::endl;
	std::cout << "pixelCenter.y " << pixelCenter.y << std::endl;
	std::cout << "pixelCenter.z " << pixelCenter.z << std::endl;
	// Now use this value as a starting point to iterate through all of the pixels in the array

	// Need index values for storing the rays
	int x = 0; // The 0th pixel on the x-axis, will be used to store rays into the array structure
	int y = 0; // The 0th pixel on the y-axis, will be used to store rays into the array structure

	int numRaysCreated = 0; // Will keep track of how many rays have been created within the loop, to make sure that this is correct

	// Set the loop starting point for ray shooting begin at the top left corner plus the pixel center offset
	point loopPoint;
	loopPoint.x = cornerPoints[0].x + pixelCenter.x;
	loopPoint.y = cornerPoints[0].y - pixelCenter.y;
	loopPoint.z = pixelCenter.z;

	// Create the ray shooting loop
	while (loopPoint.y > cornerPoints[2].y) {
		while (loopPoint.x < cornerPoints[1].x) {
			// Create the ray direction from the loopPoint, and the origin from cameraOrigin, then normalize the direction
			ray ray;
			ray.origin = cameraOrigin;
			ray.direction = loopPoint;

			// Normalize the direction of the ray
			ray.direction = getUnitVector(ray.direction);

			// Set t to 0, initially
			ray.t = 0;

			ray.isSetT = false; // This values needs to be set to false initially, or t will never get set

			// Initialize the floats for phong lighting to the ambient light upon initialization to save some effort later.
			// Red calculations
			ray.collisionColorFloat.red = ambient_light[0];
			// Green calculations
			ray.collisionColorFloat.green = ambient_light[1];
			// Blue Calculations
			ray.collisionColorFloat.blue = ambient_light[2];

			// Now store this ray into the array data structure, which maps to the image plane (i.e., the first data index [0][0] is the ray shooting through the top left pixel
			rays[x][y] = ray;

			// Increment the two incrementors
			x++; // Increment the x array index counter
			loopPoint.x += pixelTraversalValueX; // Decrement the y pixel counter	// Increment the x pixel counter

			// Increment numRaysCreated
			numRaysCreated++;
		}
		y++; // Increment the y array index counter
		loopPoint.y -= pixelTraversalValueY; // Decrement the y pixel counter
		loopPoint.x = cornerPoints[0].x + pixelCenter.x; // Reset the x pixel center
		x = 0; // Reset the x array index counter
	}

	// Print this value, make sure that it prints correctly
	std::cout << "-----RESULTS OF STEP ONE-----" << std::endl;
	std::cout << "Number of rays created: " << numRaysCreated << std::endl;

	// Test, iterate through the array to make sure that there are no empty spots
	
	//std::cout << "Ray Stats: " << std::endl;
	//for (x = 0; x < WIDTH; x++) {
	//	for (y = 0; y < HEIGHT; y++) {				
	//		std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].direction.x " << rays[x][y].direction.x << std::endl;
	//		std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].direction.y " << rays[x][y].direction.y << std::endl;
	//		std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].direction.z " << rays[x][y].direction.z << std::endl;
	//		std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].vectD  Magnitude: " << getDotProduct(rays[x][y].direction, rays[x][y].direction) << std::endl;
	//		std::cout << "rays[" << x << "][" << y << "].origin.x " << rays[x][y].origin.x << std::endl;
	//		std::cout << "rays[" << x << "][" << y << "].origin.y " << rays[x][y].origin.y << std::endl;
	//		std::cout << "rays[" << x << "][" << y << "].origin.z " << rays[x][y].origin.z << std::endl;
	//		std::cout << "rays[" << x << "][" << y << "].t " << rays[x][y].t << std::endl;
	//	}
	//}
}

/*STEP ONE FUNCTIONS END*/

/*STEP TWO FUNCTIONS START*/

void doStepTwo() {
	std::cout << "-----Completing Step Two-----"<< std::endl;
	checkCollisionsSpheres();
	checkCollisionsPolygons();
}
void checkCollisionsSpheres() {

	int numSphereCollisions = 0; // Value that tallies all of the sphere collisions

	// Loop through all of the spheres
	for (int i = 0; i < num_spheres; i++) {
		for (int x = 0; x < WIDTH * sampleNumber; x++) {
			for (int y = 0; y < HEIGHT * sampleNumber; y++) {				
				// Get xd^2 + yd^2 + zd^2 for the ray
				//double a = getDotProduct(rays[x][y].direction, rays[x][y].direction); 

				// Get 2 * ( (xd(x0 - xc)) + ((yd(y0 - yc)) + ((zd(z0 - zc)) )
				double b = 2.0 * ( (rays[x][y].direction.x * (rays[x][y].origin.x - spheres[i].position[0])) + 
								   (rays[x][y].direction.y * (rays[x][y].origin.y - spheres[i].position[1])) + 
								   (rays[x][y].direction.z * (rays[x][y].origin.z - spheres[i].position[2])) ); 
				
				// Get (x0 - xc)^2 + (y0 - yc)^2 + (z0 - zc)^2 + (r)^2
				double c = pow((rays[x][y].origin.x - spheres[i].position[0]), 2) + 
						   pow((rays[x][y].origin.y - spheres[i].position[1]), 2) + 
						   pow((rays[x][y].origin.z - spheres[i].position[2]), 2) - 
						   pow((spheres[i].radius), 2); 

				//std::cout << (pow(b, 2) - (4*c) < 0) << std::endl;

				// Now that the values have been initialized, let's make sure that b^2 -4c is NOT negative
				if (pow(b, 2) - (4*c) > 0) { // Do the calculation, else, continue to the next loop iteration without any calculation					
					double t_0 = getQuadraticFormula(true,b,c); 
					double t_1 = getQuadraticFormula(false,b,c);
					// std::cout << "Intersection with ray [" << x << "][" << y << "] and sphere[" << i << "] with t values " << t_0 << " & " << t_1 << "!" << std::endl;
					if ((t_0 < t_1)) { // If t0 is the first intersection, store that value within the ray to easily find where the collision occurred in space
						if (t_0 < rays[x][y].t || !rays[x][y].isSetT) { // If the new t is closer to the camera than the old t, replace the old t with the new one, or the ray's t needs to be set for the first time
							rays[x][y].t = t_0;
							rays[x][y].isSetT = true;
							rays[x][y].collisionShape = SPHERE;
							rays[x][y].collisionIndex = i;
						}
					}
					else { // If t1 is the first intersection, store that value within the ray to easuly find where the collision occurred in space
						if (t_1 < rays[x][y].t || !rays[x][y].isSetT) { // If the new t is closer to the camera than the old t, replace the old t with the new one, or the ray's t needs to be set for the first time
							rays[x][y].t = t_1;
							rays[x][y].isSetT = true;
							rays[x][y].collisionShape = SPHERE;
							rays[x][y].collisionIndex = i;
						}
					}
					// Calculate the collision point for reference
					rays[x][y].collisionPoint.x = rays[x][y].origin.x + (rays[x][y].direction.x) * rays[x][y].t;
					rays[x][y].collisionPoint.y = rays[x][y].origin.y + (rays[x][y].direction.y) * rays[x][y].t;
					rays[x][y].collisionPoint.z = rays[x][y].origin.z + (rays[x][y].direction.z) * rays[x][y].t;

					// Note that these values may be overwritten in the triangle test if there is indeed a value smaller than the current t that comes from that test
					numSphereCollisions++;
				}
			}
		}
	}

	// Print this value, make sure that it prints correctly
	std::cout << "-----RESULTS OF STEP TWO: PART ONE-----" << std::endl;
	std::cout << "Number of Sphere Collisions: " << numSphereCollisions << std::endl;

	// Test, iterate through the array to make sure that there are no empty spots
	
	//std::cout << "Ray Stats (t): " << std::endl;
	//for (int x = 0; x < WIDTH; x++) {
	//	for (int y = 0; y < HEIGHT; y++) {				
	//		if (rays[x][y].t > 0) {
	//			std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].direction.x " << rays[x][y].direction.x << std::endl;
	//			std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].direction.y " << rays[x][y].direction.y << std::endl;
	//			std::cout << std::setprecision(16) << "rays[" << x << "][" << y << "].direction.z " << rays[x][y].direction.z << std::endl;				
	//			std::cout << "rays[" << x << "][" << y << "].t " << rays[x][y].t << std::endl;
	//		}
	//	}
	//}

}
void checkCollisionsPolygons() {
	int numTriangleCollisions = 0; // Value that tallies all of the sphere collisions

	// Loop through all of the spheres
	for (int i = 0; i < num_triangles; i++) {
		for (int x = 0; x < WIDTH * sampleNumber; x++) {
			for (int y = 0; y < HEIGHT * sampleNumber; y++) {
				///*
				// First, find the d value for the plane that the triangle exists on
				// Create vector from Point A->B (B-A) and Point A->C (C-A)
				point ab;
				ab.x = triangles[i].v[1].position[0] - triangles[i].v[0].position[0];
				ab.y = triangles[i].v[1].position[1] - triangles[i].v[0].position[1];
				ab.z = triangles[i].v[1].position[2] - triangles[i].v[0].position[2];

				point ac;
				ac.x = triangles[i].v[2].position[0] - triangles[i].v[0].position[0];
				ac.y = triangles[i].v[2].position[1] - triangles[i].v[0].position[1];
				ac.z = triangles[i].v[2].position[2] - triangles[i].v[0].position[2];
				//*/

				// Get the Cross of these two points, which contains the A,B,C values for the plane
				point normal = getCrossProduct(ab,ac);

				// Now let's find -d
				double d = (normal.x * triangles[i].v[0].position[0]) + (normal.y * triangles[i].v[0].position[1]) + (normal.z * triangles[i].v[0].position[2]);
				d *= -1; // Make the d go to the appropriate side of the equation

				/*
				// Method adapted from http://nehe.gamedev.net/tutorial/shadows/16010/ -- may be used if my current method does not work with all Tris
				// Find the equation of a plane
				point planeXYZ; 
				planeXYZ.x = v0.y * (v1.z - v2.z) + v1.y * (v2.z - v0.z) + v2.y * (v0.z - v1.z);
				planeXYZ.y = v0.z * (v1.x - v2.x) + v1.z * (v2.x - v0.x) + v2.z * (v0.x - v1.x);
				planeXYZ.z = v0.x * (v1.y - v2.y) + v1.x * (v2.y - v0.y) + v2.x * (v0.y - v1.y);
				
				
				// Now let's find d
				d = -( v0. x *( v1.y * v2.z - v2.y * v1.z ) + v1.x * (v2.y * v0.z - v0.y * v2.z) + v2.x * (v0.y * v1.z - v1.y * v0.z) );
				*/

				// Since all triangles appear within one plane, it is safe to assume that the normals at each vertex (and throughout the test of the triangle) are the same
				// Create points for the plane equation and the normal
				// Now attempt to find a t that satisfies the triangle collision equation

				if (abs(getDotProduct(normal, rays[x][y].direction)) < 0.0001) { // Calculation needs to abort to prevent divide by 0 error if true
					continue;
				}

				double t = -(getDotProduct(normal, rays[x][y].origin) + d)/(getDotProduct(normal, rays[x][y].direction));

				// Now that the values have been initialized, let's make sure that b^2 -4c is NOT negative
				if (t > 0) { // Do the calculation, else, continue to the next loop iteration without any calculation				

					// Calculate the collision point for reference
					point planeCollisionPoint;
					planeCollisionPoint.x = rays[x][y].origin.x + (rays[x][y].direction.x) * t;
					planeCollisionPoint.y = rays[x][y].origin.y + (rays[x][y].direction.y) * t;
					planeCollisionPoint.z = rays[x][y].origin.z + (rays[x][y].direction.z) * t;

					// Since the collision point lies within the plane, check to see if this point actually resides within the triangle
				
					point v0; 
					point v1; 
					point v2; 					

					v0.x = triangles[i].v[0].position[0];
					v0.y = triangles[i].v[0].position[1];
					v0.z = triangles[i].v[0].position[2];

					v1.x = triangles[i].v[1].position[0];
					v1.y = triangles[i].v[1].position[1];
					v1.z = triangles[i].v[1].position[2];

					v2.x = triangles[i].v[2].position[0];
					v2.y = triangles[i].v[2].position[1];
					v2.z = triangles[i].v[2].position[2];

					// Get all of the areas for the barycentric coordinates
					double alpha   = getTriangleArea(planeCollisionPoint, v1, v2)/getTriangleArea(v0,v1,v2);
					double beta    = getTriangleArea(v0, planeCollisionPoint, v2)/getTriangleArea(v0,v1,v2);
					double charlie = 1 - alpha - beta; // Hack out the third coordinate from the other two to save on FP divisions

					if (alpha > 0 && beta > 0 && charlie > 0 && (abs(1 - alpha - beta - charlie) < 0.0001)) { // If all points are positive
						if (t < rays[x][y].t || rays[x][y].isSetT == false) { // If t < current t, or if t is not set
							rays[x][y].t = t;
							rays[x][y].isSetT = true;
							rays[x][y].collisionShape = TRIANGLE;
							rays[x][y].collisionIndex = i;

							// Calculate the collision point for reference
							rays[x][y].collisionPoint.x = planeCollisionPoint.x;
							rays[x][y].collisionPoint.y = planeCollisionPoint.y;
							rays[x][y].collisionPoint.z = planeCollisionPoint.z;

							// Hold the BaryCentric ratios for further use with shading
							rays[x][y].barycentricRatios.x = alpha;
							rays[x][y].barycentricRatios.y = beta;
							rays[x][y].barycentricRatios.z = charlie;

							// Note that these values may be overwritten in the triangle test if there is indeed a value smaller than the current t that comes from that test
							numTriangleCollisions++;
						}
					}
				}
			}
		}
	}
	// Print this value, make sure that it prints correctly
	std::cout << "-----RESULTS OF STEP TWO: PART TWO-----" << std::endl;
	std::cout << "Number of Triangle Collisions: " << numTriangleCollisions << std::endl;
}

/*STEP TWO FUNCTIONS END*/

/*STEP THREE FUNCTIONS START*/

void doStepThree() { 
	std::cout << "-----Completing Step Three-----"<< std::endl;
	calculateNormals();
	calculateColor();
	convertColorValues();
}

void calculateNormals() {
	int numSphereNormals = 0;
	int numTriangleNormals = 0;
	
	for (int x = 0; x < WIDTH * sampleNumber; x++) {
		for (int y = 0; y < HEIGHT * sampleNumber; y++) {
			// Calculate the normals for each rays that has collided with something (check the isSetT bool)
			if (rays[x][y].isSetT) { // Then a collision at point t was recorded
				if (rays[x][y].collisionShape == SPHERE) { // Then this ray collided with a sphere, do the sphere normal calculation
					rays[x][y].collisionNormal.x = (rays[x][y].collisionPoint.x - spheres[rays[x][y].collisionIndex].position[0]) / spheres[rays[x][y].collisionIndex].radius;
					rays[x][y].collisionNormal.y = (rays[x][y].collisionPoint.y - spheres[rays[x][y].collisionIndex].position[1]) / spheres[rays[x][y].collisionIndex].radius;
					rays[x][y].collisionNormal.z = (rays[x][y].collisionPoint.z - spheres[rays[x][y].collisionIndex].position[2]) / spheres[rays[x][y].collisionIndex].radius;
					// Normalize the collisonNormal
					rays[x][y].collisionNormal = getUnitVector(rays[x][y].collisionNormal);
					numSphereNormals++;
				}
				else if (rays[x][y].collisionShape == TRIANGLE) { // Then this ray collided with a triangle, do the triangle normal calculation
					// Use barycentric coordinates to find the interpolated norms of the vertices
					rays[x][y].collisionNormal.x = (rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].normal[0]) +
						                           (rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].normal[0]) +
												   (rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].normal[0]);

					rays[x][y].collisionNormal.y = (rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].normal[1]) +
												   (rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].normal[1]) +
				    							   (rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].normal[1]);

					rays[x][y].collisionNormal.z = (rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].normal[2]) +
						                           (rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].normal[2]) +
												   (rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].normal[2]);
					// Normalize the collisonNormal
					rays[x][y].collisionNormal = getUnitVector(rays[x][y].collisionNormal);
					numTriangleNormals++;
				}
			}
		}
	}
	// Print this value, make sure that it prints correctly
	std::cout << "-----RESULTS OF STEP THREE: PART ONE-----" << std::endl;
	std::cout << "Number of Sphere Normal Calculations: " << numSphereNormals << std::endl;
	std::cout << "Number of Triangle Normal Calculations: " << numTriangleNormals << std::endl;
}

void calculateColor() {
	for (int x = 0; x < WIDTH * sampleNumber; x++) {
		for (int y = 0; y < HEIGHT * sampleNumber; y++) {
			// Calculate the normals for each rays that has collided with something (check the isSetT bool)
			if (rays[x][y].isSetT) { // Then a collision at point t was recorded
				// Run through all of the light sources and cast shadow rays
				for (int i = 0; i < num_lights; i++) {
					// Make a shadowRay pointer, it will be needed for the lighting calculations
					ray *shadowRay = new ray;
					shadowRay->isSetT = false; // Make sure that this value is set to false, or else things will not work
					calculateShadowRay(&rays[x][y], lights[i].position, shadowRay);
					if (shadowRay->isSetT == false) { // If there were no collisions with the shadow ray, then calculate lighting for this light
						// Do light calculations here
						point v; // Set up a point to be used as the v vector of the Phong lighting equation (which is really just the reverse direction of rays[x][y])
						v.x = -rays[x][y].direction.x;
						v.y = -rays[x][y].direction.y;
						v.z = -rays[x][y].direction.z;

						// Get the reflected vector from the shadow ray, which is l
						point r;
						// Step one, calculate l dot n
						double dotResultLN = getDotProduct(shadowRay->direction, rays[x][y].collisionNormal);

						// Step two, calculate n - l
						point pointN_L;
						pointN_L.x = rays[x][y].collisionNormal.x - shadowRay->direction.x;
						pointN_L.y = rays[x][y].collisionNormal.y - shadowRay->direction.y;
						pointN_L.z = rays[x][y].collisionNormal.z - shadowRay->direction.z;

						// Step three, multiply everything together
						r.x = (2 * (dotResultLN * rays[x][y].collisionNormal.x)) - shadowRay->direction.x;
						r.y = (2 * (dotResultLN * rays[x][y].collisionNormal.y)) - shadowRay->direction.y;
						r.z = (2 * (dotResultLN * rays[x][y].collisionNormal.z)) - shadowRay->direction.z;

						// Just to make sure r is normalized
						r = getUnitVector(r);

						if (getDotProduct(r,r) == 0) {
							r = rays[x][y].collisionNormal;
						}

						// Get the attenuations calculation before calculating the colors themselves
						// Convert the light position array into a point
						point lightPosition;
						lightPosition.x = lights[i].position[0];
						lightPosition.y = lights[i].position[1];
						lightPosition.z = lights[i].position[2];

						double distanceToLight = getDistance(lightPosition, shadowRay->origin);
						double attenuationPhong = 1.0 / (attenuationValues.x + (attenuationValues.y * distanceToLight) + (attenuationValues.z * (pow(distanceToLight, 2))) );

						// Calculate the two dot products (l * n) && (r * v) and clamp them if they go below 0 or above 1
						double dotLN = dotResultLN;

						if (dotLN < 0) { // Clamp to zero for phong lighting
							dotLN = 0;
						}

						else if (dotLN > 1) { // Clamp to one for phong lighting
							dotLN = 1;
						}

						double dotRV = getDotProduct(r,v);

						if (dotRV < 0) { // Clamp to zero for phong lighting
							dotRV = 0;
						}

						else if (dotRV > 1) { // Clamp to one for phong lighting
							dotRV = 1;
						}

						// Now do color calculations
						if (rays[x][y].collisionShape == SPHERE) {
							// Red calculations
							rays[x][y].collisionColorFloat.red += /*attenuationPhong **/ lights[i].color[0] * ( 
								(spheres[rays[x][y].collisionIndex].color_diffuse[0] * dotLN) + 
								(spheres[rays[x][y].collisionIndex].color_specular[0] * pow(dotRV, spheres[rays[x][y].collisionIndex].shininess)
								));

							// Green calculations
							rays[x][y].collisionColorFloat.green += /*attenuationPhong **/ lights[i].color[1] * ( 
								(spheres[rays[x][y].collisionIndex].color_diffuse[1] * dotLN) + 
								(spheres[rays[x][y].collisionIndex].color_specular[1] * pow(dotRV, spheres[rays[x][y].collisionIndex].shininess)
								));

							// Blue Calculations
							rays[x][y].collisionColorFloat.blue += /*attenuationPhong **/ lights[i].color[2] * ( 
								(spheres[rays[x][y].collisionIndex].color_diffuse[2] * dotLN) + 
								(spheres[rays[x][y].collisionIndex].color_specular[2] * pow(dotRV, spheres[rays[x][y].collisionIndex].shininess)
								));
						}
						else if (rays[x][y].collisionShape == TRIANGLE) {
							// Red calculations
							rays[x][y].collisionColorFloat.red += /*attenuationPhong **/ lights[i].color[0] * ( 
								// Diffuse Lighting
								((
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].color_diffuse[0]) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].color_diffuse[0]) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].color_diffuse[0])
								) * 
								dotLN ) +
								// Specular Lighting
								((
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].color_specular[0]) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].color_specular[0]) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].color_specular[0])								
								) *
								pow(dotRV, (
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].shininess) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].shininess) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].shininess)	
								))) );
								
							// Green Calculations
							rays[x][y].collisionColorFloat.green += /*attenuationPhong **/ lights[i].color[1] * ( 
								// Diffuse Lighting
								((
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].color_diffuse[1]) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].color_diffuse[1]) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].color_diffuse[1])
								) * 
								dotLN ) +
								// Specular Lighting
								((
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].color_specular[1]) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].color_specular[1]) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].color_specular[1])								
								) *
								pow(dotRV, (
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].shininess) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].shininess) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].shininess)	
								))) );

							// Blue Calculations
							rays[x][y].collisionColorFloat.blue += /*attenuationPhong **/ lights[i].color[2] * ( 
								// Diffuse Lighting
								((
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].color_diffuse[2]) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].color_diffuse[2]) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].color_diffuse[2])
								) * 
								dotLN ) +
								// Specular Lighting
								((
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].color_specular[2]) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].color_specular[2]) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].color_specular[2])								
								) *
								pow(dotRV, (
									(rays[x][y].barycentricRatios.x * triangles[rays[x][y].collisionIndex].v[0].shininess) + 
									(rays[x][y].barycentricRatios.y * triangles[rays[x][y].collisionIndex].v[1].shininess) +
									(rays[x][y].barycentricRatios.z * triangles[rays[x][y].collisionIndex].v[2].shininess)	
								))) );
						}
					}
					else {
						rays[x][y].collisionColorFloat.red += 0;
						rays[x][y].collisionColorFloat.blue += 0;
						rays[x][y].collisionColorFloat.green += 0;
					}
					// At the end of all of this, make sure to delete the shadow ray to prevent a memory leak
					delete shadowRay;
				}
				// End Clamping
				// Make sure to clamp the phong values if they are below 0 or above 1
				if (rays[x][y].collisionColorFloat.red < 0) { // Clamp to 0
					rays[x][y].collisionColorFloat.red = 0;
				}
				else if (rays[x][y].collisionColorFloat.red > 1) { // Clamp to 1
					rays[x][y].collisionColorFloat.red = 1;
				}

				if (rays[x][y].collisionColorFloat.green < 0) { // Clamp to 0
					rays[x][y].collisionColorFloat.green = 0;
				}
				else if (rays[x][y].collisionColorFloat.green > 1) { // Clamp to 1
					rays[x][y].collisionColorFloat.green = 1;
				}

				if (rays[x][y].collisionColorFloat.blue < 0) { // Clamp to 0
					rays[x][y].collisionColorFloat.blue = 0;
				}
				else if (rays[x][y].collisionColorFloat.blue > 1) { // Clamp to 1
					rays[x][y].collisionColorFloat.blue = 1;
				}
			}
			else {
				rays[x][y].collisionColorFloat.red = 1;
				rays[x][y].collisionColorFloat.green = 1;
				rays[x][y].collisionColorFloat.blue = 1;
			}
		}
	}
	std::cout << "Step Three: Finished calculating colors!" << std::endl;
}

void calculateShadowRay(ray *collisionRay, double lightCollisionPoint[], ray *shadowRay) {

	// Set the origin of the shadowRay to the origin of the collision
	shadowRay->origin.x = collisionRay->collisionPoint.x;
	shadowRay->origin.y = collisionRay->collisionPoint.y;
	shadowRay->origin.z = collisionRay->collisionPoint.z;

	// Get the two unit rays that goes from:
	// a) The camera to the light source -- needs to be calculated here
	// b) The camera to the collision point -- attained through the first function argument

	ray *cameraToLight = new ray;

	cameraToLight->origin = cameraOrigin;
	cameraToLight->direction.x = lightCollisionPoint[0] - cameraOrigin.x;
	cameraToLight->direction.y = lightCollisionPoint[1] - cameraOrigin.y;
	cameraToLight->direction.z = lightCollisionPoint[2] - cameraOrigin.z;

	// Normalize the direction vector
	cameraToLight->direction = getUnitVector(cameraToLight->direction);

	// Set the direction of the ray towards the light source (through (b) - (a)), and then normalize it again just to be sure
	shadowRay->direction.x = cameraToLight->direction.x - collisionRay->direction.x;
	shadowRay->direction.y = cameraToLight->direction.y - collisionRay->direction.y;
	shadowRay->direction.z = cameraToLight->direction.z - collisionRay->direction.z;

	// Make sure to delete the cameraToLight ray, since it is not needed anymore
	delete cameraToLight;

	// Set the direction of the ray towards the light source (through (b) - (a)), and then normalize it again just to be sure
	shadowRay->direction.x = lightCollisionPoint[0] - shadowRay->origin.x;
	shadowRay->direction.y = lightCollisionPoint[1] - shadowRay->origin.y;
	shadowRay->direction.z = lightCollisionPoint[2] - shadowRay->origin.z;

	// Normalize the shadowRay direction
	shadowRay->direction = getUnitVector(shadowRay->direction);

	// Now check to see if this shadow ray collided with anything
	checkShadowCollisionsSpheres(shadowRay, lightCollisionPoint); 
	checkShadowCollisionsTriangles(shadowRay, lightCollisionPoint);

	return;
}

bool checkShadowCollisionsSpheres(ray *shadowRay, double lightCollisionPoint[]) {

	for (int i = 0; i < num_spheres; i++) {
		// Get 2 * ( (xd(x0 - xc)) + ((yd(y0 - yc)) + ((zd(z0 - zc)) )
		double b = 2.0 * ( (shadowRay->direction.x * (shadowRay->origin.x - spheres[i].position[0])) + 
							(shadowRay->direction.y * (shadowRay->origin.y - spheres[i].position[1])) + 
							(shadowRay->direction.z * (shadowRay->origin.z - spheres[i].position[2])) ); 
				
		// Get (x0 - xc)^2 + (y0 - yc)^2 + (z0 - zc)^2 + (r)^2
		double c = pow((shadowRay->origin.x - spheres[i].position[0]), 2) + 
					pow((shadowRay->origin.y - spheres[i].position[1]), 2) + 
					pow((shadowRay->origin.z - spheres[i].position[2]), 2) - 
					pow((spheres[i].radius), 2); 

		// Now that the values have been initialized, let's make sure that b^2 -4c is NOT negative
		if (pow(b, 2) - (4*c) > 0) { // If there is an intersection, return true				
			// Calculate the point of collision
			double t_0 = getQuadraticFormula(true,b,c); 
			double t_1 = getQuadraticFormula(false,b,c);
			// If t_0 or t_1 are NOT > 0, then go into the next loop iteration
			if (t_0 <= 0.00001 && t_1 <= 0.00001) { // There will always be a tangental collision with one of the spheres at t = 0 for the shadow ray, and this case must be nullified, since !(t > 0)
				continue;
			}

			// Check to make sure that at least one object is in front of the light source
			point planeCollisionPoint_0;
			planeCollisionPoint_0.x = shadowRay->origin.x + (shadowRay->direction.x) * t_0;
			planeCollisionPoint_0.y = shadowRay->origin.y + (shadowRay->direction.y) * t_0;
			planeCollisionPoint_0.z = shadowRay->origin.z + (shadowRay->direction.z) * t_0;

			point planeCollisionPoint_1;
			planeCollisionPoint_1.x = shadowRay->origin.x + (shadowRay->direction.x) * t_1;
			planeCollisionPoint_1.y = shadowRay->origin.y + (shadowRay->direction.y) * t_1;
			planeCollisionPoint_1.z = shadowRay->origin.z + (shadowRay->direction.z) * t_1;

			// Convert the light position array into a point
			point lightPosition;
			lightPosition.x = lightCollisionPoint[0];
			lightPosition.y = lightCollisionPoint[1];
			lightPosition.z = lightCollisionPoint[2];
			
			// Make sure to check that object collision in question occurred IN FRONT of the light source instead of behind it.  
			// Check to see if the distance of the collision with this current sphere is closer that the distance of the light
			if (!(getDistance(planeCollisionPoint_0, shadowRay->origin) < getDistance(lightPosition, shadowRay->origin))
				&& 
				!(getDistance(planeCollisionPoint_1, shadowRay->origin) < getDistance(lightPosition, shadowRay->origin))
			) continue;		

			// If one of these collisions is at t > 0, though, then there really was a collision, and this needs to be accounted for
			
			if ((t_0 < t_1) && t_0 > 0.00001) { // If t0 is the first intersection, store that value within the ray to easily find where the collision occurred in space
				if (t_0 < shadowRay->t || !shadowRay->isSetT) { // If the new t is closer to the camera than the old t, replace the old t with the new one, or the ray's t needs to be set for the first time
					shadowRay->t = t_0;
					shadowRay->isSetT = true;
					shadowRay->collisionShape = SPHERE;
					shadowRay->collisionIndex = i;
				}
			}
			else if (t_1 > 0.00001) { // If t1 is the first intersection, store that value within the ray to easuly find where the collision occurred in space
				if (t_1 < shadowRay->t || !shadowRay->isSetT) { // If the new t is closer to the camera than the old t, replace the old t with the new one, or the ray's t needs to be set for the first time
					shadowRay->t = t_1;
					shadowRay->isSetT = true;
					shadowRay->collisionShape = SPHERE;
					shadowRay->collisionIndex = i;
				}
			}

		
		}	
	}
	return false; // No collisions
}

bool checkShadowCollisionsTriangles(ray *shadowRay, double lightCollisionPoint[]) {
	// Loop through all of the spheres
	for (int i = 0; i < num_triangles; i++) {
		// First, find the d value for the plane that the triangle exists on
		// Create vector from Point A->B (B-A) and Point A->C (C-A)
		point ab;
		ab.x = triangles[i].v[1].position[0] - triangles[i].v[0].position[0];
		ab.y = triangles[i].v[1].position[1] - triangles[i].v[0].position[1];
		ab.z = triangles[i].v[1].position[2] - triangles[i].v[0].position[2];

		point ac;
		ac.x = triangles[i].v[2].position[0] - triangles[i].v[0].position[0];
		ac.y = triangles[i].v[2].position[1] - triangles[i].v[0].position[1];
		ac.z = triangles[i].v[2].position[2] - triangles[i].v[0].position[2];

		// Get the Cross of these two points, which contains the A,B,C values for the plane
		point normal = getCrossProduct(ab,ac);

		// Now let's find -d
		double d = (normal.x * triangles[i].v[0].position[0]) + (normal.y * triangles[i].v[0].position[1]) + (normal.z * triangles[i].v[0].position[2]);
		d *= -1; // Make the d go to the appropriate side of the equation

		// Since all triangles appear within one plane, it is safe to assume that the normals at each vertex (and throughout the test of the triangle) are the same
		// Create points for the plane equation and the normal
		// Now attempt to find a t that satisfies the triangle collision equation

		if (abs(getDotProduct(normal, shadowRay->direction)) < 0.0001) { // Calculation needs to abort to prevent divide by 0 error if true
			continue;
		}

		double t = -(getDotProduct(normal, shadowRay->origin) + d)/(getDotProduct(normal, shadowRay->direction));

		// Now that the values have been initialized, let's make sure that b^2 -4c is NOT negative
		if (t > 0.0001) { // Do the calculation, else, continue to the next loop iteration without any calculation				

			// Check to make sure that at least object is in front of the light source
			point planeCollisionPoint;
			planeCollisionPoint.x = shadowRay->origin.x + (shadowRay->direction.x) * t;
			planeCollisionPoint.y = shadowRay->origin.y + (shadowRay->direction.y) * t;
			planeCollisionPoint.z = shadowRay->origin.z + (shadowRay->direction.z) * t;

			// Since the collision point lies within the plane, check to see if this point actually resides within the triangle				
			point v0; 
			point v1; 
			point v2;

			v0.x = triangles[i].v[0].position[0];
			v0.y = triangles[i].v[0].position[1];
			v0.z = triangles[i].v[0].position[2];

			v1.x = triangles[i].v[1].position[0];
			v1.y = triangles[i].v[1].position[1];
			v1.z = triangles[i].v[1].position[2];

			v2.x = triangles[i].v[2].position[0];
			v2.y = triangles[i].v[2].position[1];
			v2.z = triangles[i].v[2].position[2];

			// Get all of the areas for the barycentric coordinates
			double alpha   = getTriangleArea(planeCollisionPoint, v1, v2)/getTriangleArea(v0,v1,v2);
			double beta    = getTriangleArea(v0, planeCollisionPoint, v2)/getTriangleArea(v0,v1,v2);
			double charlie = 1 - alpha - beta; // Hack out the third coordinate from the other two to save on FP divisions

			if (alpha > 0 && beta > 0 && charlie > 0 && (abs(1 - alpha - beta - charlie) < 0.0001) && t > 0.00001) { // If all points are positive
				// Convert the light position array into a point
				point lightPosition;
				lightPosition.x = lightCollisionPoint[0];
				lightPosition.y = lightCollisionPoint[1];
				lightPosition.z = lightCollisionPoint[2];

				point lightToOrigin;
				lightToOrigin.x = lightPosition.x - shadowRay->origin.x;
				lightToOrigin.y = lightPosition.y - shadowRay->origin.y;
				lightToOrigin.z = lightPosition.z - shadowRay->origin.z;

				point collisionToOrigin;
				collisionToOrigin.x = planeCollisionPoint.x - shadowRay->origin.x;
				collisionToOrigin.y = planeCollisionPoint.y - shadowRay->origin.y;
				collisionToOrigin.z = planeCollisionPoint.z - shadowRay->origin.z;

				point c1 = getUnitVector(collisionToOrigin);

				point c2 = getUnitVector(lightToOrigin);

				if ( getDotProduct(c1,c1) == getDotProduct(c2,c2) ) {
					//std::cout << std::endl;
				}					
			
				// Make sure to check that object collision in question occurred IN FRONT of the light source instead of behind it.  
				// Check to see if the distance of the collision with this current sphere is closer that the distance of the light
				if (!(getDotProduct(collisionToOrigin, collisionToOrigin) < getDotProduct(lightToOrigin, lightToOrigin))) // If so, return true
					continue;
				
				if (t < shadowRay->t || shadowRay->isSetT == false) { // If t < current t, or if t is not set
					shadowRay->t = t;
					shadowRay->isSetT = true;
					shadowRay->collisionShape = TRIANGLE;
					shadowRay->collisionIndex = i;

					// Calculate the collision point for reference
					shadowRay->collisionPoint.x = planeCollisionPoint.x;
					shadowRay->collisionPoint.y = planeCollisionPoint.y;
					shadowRay->collisionPoint.z = planeCollisionPoint.z;

					// Hold the BaryCentric ratios for further use with shading
					shadowRay->barycentricRatios.x = alpha;
					shadowRay->barycentricRatios.y = beta;
					shadowRay->barycentricRatios.z = charlie;
				}
			}
		}
	}
	return false;
}

void convertColorValues() {
	if (sampleNumber == 1) {
		for (int x = 0; x < WIDTH; x++) {
			for (int y = 0; y < HEIGHT; y++) {
				if (rays[x][y].isSetT) {
					// Finally, convert these floats into the char format
					rays[x][y].collisionColor.red = rays[x][y].collisionColorFloat.red * 255;
					rays[x][y].collisionColor.green = rays[x][y].collisionColorFloat.green * 255;
					rays[x][y].collisionColor.blue = rays[x][y].collisionColorFloat.blue * 255;
				}
				else {
					rays[x][y].collisionColor.red = 255;
					rays[x][y].collisionColor.green = 255;
					rays[x][y].collisionColor.blue = 255;
				}
			}
		}
	}
	else { // Get all of the rays from a block of sampleNumber size, add up their colors, and then divide by sampleNum squared
		for (int x = 0; x < WIDTH * sampleNumber; x+=sampleNumber) {
			for (int y = 0; y < HEIGHT * sampleNumber; y+=sampleNumber) {
				point colorSum; // Will hold the sum of all the colors for the next two loops
				colorSum.x = 0;
				colorSum.y = 0;
				colorSum.z = 0;
				for (int xNum = x; xNum < x + sampleNumber; xNum++) {
					for (int yNum = y; yNum < y + sampleNumber; yNum++) {
						colorSum.x += rays[xNum][yNum].collisionColorFloat.red;
						colorSum.y += rays[xNum][yNum].collisionColorFloat.green;
						colorSum.z += rays[xNum][yNum].collisionColorFloat.blue;
					}
				}
				// Then average the values with sampleNum squared
				colorSum.x /= (sampleNumber*sampleNumber);
				colorSum.y /= (sampleNumber*sampleNumber);
				colorSum.z /= (sampleNumber*sampleNumber);

				// Make sure to clamp the phong values if they are below 0 or above 1
				if (colorSum.x < 0) { // Clamp to 0
					colorSum.x = 0;
				}
				else if (colorSum.x > 1) { // Clamp to 1
					colorSum.x = 1;
				}
				if (colorSum.y < 0) { // Clamp to 0
					colorSum.y = 0;
				}
				else if (colorSum.y > 1) { // Clamp to 1
					colorSum.y = 1;
				}
				if (colorSum.z < 0) { // Clamp to 0
					colorSum.z = 0;
				}
				else if (colorSum.z > 1) { // Clamp to 1
					colorSum.z = 1;
				}

				// Initialize this color to the ambient color
				averagedColors[x/sampleNumber][y/sampleNumber].red = ambient_light[0];
				averagedColors[x/sampleNumber][y/sampleNumber].green = ambient_light[1];
				averagedColors[x/sampleNumber][y/sampleNumber].blue = ambient_light[2];

				// Now put the results of this into the averageRays array
				averagedColors[x/sampleNumber][y/sampleNumber].red = colorSum.x;
				averagedColors[x/sampleNumber][y/sampleNumber].green = colorSum.y;
				averagedColors[x/sampleNumber][y/sampleNumber].blue = colorSum.z;
			}
		}
	}
}

/*STEP THREE FUNCTIONS END*/


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
		if (sampleNumber == 1) {
			plot_pixel(x, HEIGHT-y-1, rays[x][y].collisionColor.red, rays[x][y].collisionColor.green, rays[x][y].collisionColor.blue);
		}
		else {
			char red = averagedColors[x][y].red * 255;
			char green = averagedColors[x][y].green * 255;
			char blue = averagedColors[x][y].blue * 255;
			plot_pixel(x, HEIGHT-y-1, red, green, blue);
		}
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

	  // Add in the method to create the sub lights for soft shadows
	  for (int i = 0; i < numRandomLights; i++) {
		  Light light = l;
		  // Randomly perturb the position of this new light
		  /*
		  light.position[0] += rand() % 1000;
		  light.position[0] -= 500;
		  light.position[0] /= 100000;
		  
		  light.position[1] += rand() % 1000;
		  light.position[1] -= 500;
		  light.position[1] /= 100000;
		  
		  light.position[2] += rand() % 1000;
		  light.position[2] -= 500;
		  light.position[2] /= 100000;
		  */

		  ///*
		  for (int j = 0; j < 3; j++) {
			  int randOffset = (rand() % 1001) - 500;
			  double offset = (double)randOffset/(double)10000;
			  light.position[j] += offset;
		  }
		  //*/
		  

		  // Weaken the intensity of the light proportionally to the number of satellite lights being created
		  //light.color[0] /= (numRandomLights + 1);
		  //light.color[1] /= (numRandomLights + 1);
		  //light.color[2] /= (numRandomLights + 1);

		  lights[num_lights++] = light;
	  }

	  /*
	  // Divide the intensity of tha main light by the number of random lights + 1 after the after loop finishes executing
	  l.color[0] /= (numRandomLights + 1);
	  l.color[1] /= (numRandomLights + 1);
	  l.color[2] /= (numRandomLights + 1);
	  */

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

  // Set any variables that needs to be set for rayTracing here
  // Set the attenuation variables for the phong calculations
  attenuationValues.x = 1;
  attenuationValues.y = 1;
  attenuationValues.z = 1;

  // Seed rand
  srand(time(NULL));

  // Do the raytracing calculations
  doStepOne(); // Uniformly send out rays from one location
  doStepTwo(); // Check for intersections with triangles and spheres
  doStepThree(); // Complete all of the lighting functionality of the raytracer here
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