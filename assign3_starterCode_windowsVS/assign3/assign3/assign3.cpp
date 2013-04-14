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

// Struct that will hold the components of a color from 0 to 255
struct color {
	unsigned char red; 
	unsigned char green;
	unsigned char blue;
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
	color collisionColor; // Will hold the color 
	int collisionIndex; // Which shape did the ray collide with, along with the shape type aforementioned?
	point collisionPoint; // Where the collision actually occurred for this ray, will be easier to reference this way, as it will be used more than once
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

// Function protoTypes:
point getCrossProduct(point a, point b);
point getUnitVector(point a);
double getDotProduct(point a, point b);
double getDistance(point a, point b);
double getQuadraticFormula(boolean isLowerVal, double b, double c);

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
	return sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y), 2));
}

double getQuadraticFormula(boolean isLowerVal, double b, double c) { // Will use the quadratic value to return the two points of the sphere intersection
	double t = 0;
	if (isLowerVal)
		t = (-b - ( sqrt(pow(b, 2) - (4*c)) )) / 2; // Get the lower value of the quadratic formula
	else
		t = (-b + ( sqrt(pow(b, 2) - (4*c)) )) / 2; // Get the higher value of the quadratic formula
	return t;
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
bool calculateShadowRay(ray *collisionRay, double lightCollisionPoint[], ray *shadowRay); // Will return whether the shadow ray collided with an object or not, may need a check for when the origin of a ray is inside a sphere
bool checkShadowCollisionsSpheres(ray *shadowRay, double lightCollisionPoint[]); // Will check all of the spheres to see if there is a collision with a shadow ray
bool checkShadowCollisionsTriangles(ray *shadowRay, double lightCollisionPoint[]); // Will check all of the triangles to see if there is a collision with a shadow ray

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

	pixelTraversalValueX = (double)pixelTraversalValueX/WIDTH;
	pixelTraversalValueY = (double)pixelTraversalValueY/HEIGHT;

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
		for (int x = 0; x < WIDTH; x++) {
			for (int y = 0; y < HEIGHT; y++) {
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

}

/*STEP TWO FUNCTIONS END*/

/*STEP THREE FUNCTIONS START*/

void doStepThree() { 
	std::cout << "-----Completing Step THREE-----"<< std::endl;
	calculateNormals();
	calculateColor();
}

void calculateNormals() {
	int numSphereNormals = 0;
	int numTriangleNormals = 0;
	
	for (int x = 0; x < WIDTH; x++) {
		for (int y = 0; y < HEIGHT; y++) {
			// Calculate the normals for each rays that has collided with something (check the isSetT bool)
			if (rays[x][y].isSetT) { // Then a collision at point t was recorded
				if (rays[x][y].collisionShape == SPHERE) { // Then this ray collided with a sphere, do the sphere normal calculation
					rays[x][y].collisionNormal.x = (rays[x][y].collisionPoint.x - spheres[rays[x][y].collisionIndex].position[0]) / spheres[rays[x][y].collisionIndex].radius;
					rays[x][y].collisionNormal.y = (rays[x][y].collisionPoint.y - spheres[rays[x][y].collisionIndex].position[1]) / spheres[rays[x][y].collisionIndex].radius;
					rays[x][y].collisionNormal.z = (rays[x][y].collisionPoint.z - spheres[rays[x][y].collisionIndex].position[2]) / spheres[rays[x][y].collisionIndex].radius;
					numSphereNormals++;
				}
				else if (rays[x][y].collisionShape == TRIANGLE) { // Then this ray collided with a triangle, do the triangle normal calculation
					continue; // Will add in content when the time comes
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
	for (int x = 0; x < WIDTH; x++) {
		for (int y = 0; y < HEIGHT; y++) {
			// Calculate the normals for each rays that has collided with something (check the isSetT bool)
			if (rays[x][y].isSetT) { // Then a collision at point t was recorded
				// Run through all of the light sources and cast shadow rays
				for (int i = 0; i < num_lights; i++) {
					// Make a shadowRay pointer, it will be needed for the lighting calculations
					ray *shadowRay = new ray;
					if (calculateShadowRay(&rays[x][y], lights[i].position, shadowRay)) { // If there were no collisions with the shadow ray, then calculate lighting for this light
						// Do light calculations here
					}

					// At the end of all of this, make sure to delete the shaow ray to prevent a memory leak
					delete shadowRay;
				}
			}
		}
	}
}

bool calculateShadowRay(ray *collisionRay, double lightCollisionPoint[], ray *shadowRay) {

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

	shadowRay->direction = getUnitVector(shadowRay->direction);

	// Now check to see if this shadow ray collided with anything
	if (checkShadowCollisionsSpheres(shadowRay, lightCollisionPoint)) return false; // If there is a collision with a sphere, return false
	if (checkShadowCollisionsTriangles(shadowRay, lightCollisionPoint)) return false; // If there is a collision with a Triangle, return false

	return true;
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
			shadowRay->collisionPoint.x = shadowRay->origin.x + (shadowRay->direction.x) * shadowRay->t;
			shadowRay->collisionPoint.y = shadowRay->origin.y + (shadowRay->direction.y) * shadowRay->t;
			shadowRay->collisionPoint.z = shadowRay->origin.z + (shadowRay->direction.z) * shadowRay->t;

			// Convert the light position array into a point
			point lightPosition;
			lightPosition.x = lightCollisionPoint[0];
			lightPosition.y = lightCollisionPoint[1];
			lightPosition.z = lightCollisionPoint[2];
			
			// Make sure to check that object collision in question occurred IN FRONT of the light source instead of behind it.  
			// Check to see if the distance of the collision with this current sphere is closer that the distance of the light
			if (getDistance(shadowRay->collisionPoint, shadowRay->origin) < getDistance(lightPosition, shadowRay->origin)) // If so, return true
				return true;
		}
	}

	return false; // No collisions
}

bool checkShadowCollisionsTriangles(ray *shadowRay, double lightCollisionPoint[]) {
	return false;
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