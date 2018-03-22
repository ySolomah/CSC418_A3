/***********************************************************
	
	Starter code for Assignment 3

	Implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade(Ray3D& ray) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	Vector3D norm = ray.intersection.normal;
	norm.normalize();

	Vector3D lightVec = pos - ray.intersection.point;
	lightVec.normalize();

	Vector3D viewVec = -ray.dir;
	viewVec.normalize();

	Vector3D reflectVec = 2 * lightVec.dot(norm) * norm - lightVec;
	reflectVec.normalize();

	Color amb = ray.intersection.mat->ambient * col_ambient;
	Color diff = fmax(0, norm.dot(lightVec)) * col_diffuse * ray.intersection.mat->diffuse;
	Color spec = fmax(0, pow(viewVec.dot(reflectVec), ray.intersection.mat->specular_exp)) * col_specular * ray.intersection.mat->specular;

	ray.col = ray.col + amb + diff + spec;
	ray.col.clamp();

}

