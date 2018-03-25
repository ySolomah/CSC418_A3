/***********************************************************
	
	Starter code for Assignment 3

	Implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade(Ray3D& ray, float attFactor) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	// assume we get world position of light source
	Vector3D light = pos - ray.intersection.point;
	light.normalize(); 

	Vector3D normal = ray.intersection.normal;
	normal.normalize();

	double angle = normal.dot(light);
	if (angle < 0){
		angle = 0;
	}
	// calculate specular component
	double specular_exp = ray.intersection.mat->specular_exp;

	Vector3D reflected_light = 2.0*(light.dot(normal))*normal-light;
	reflected_light.normalize();

	Vector3D view = -ray.dir;
	view.normalize();

	double specAngle = reflected_light.dot(view);
	if (specAngle < 0){
		specAngle = 0.0;
	}
	specAngle = pow(specAngle, specular_exp);
	
	Color ambient = ray.intersection.mat->ambient;
	Color diffuse = ray.intersection.mat->diffuse;
	Color specular = ray.intersection.mat->specular;
	
	ray.col = ray.col + attFactor * ambient;
	ray.col = ray.col + attFactor *  diffuse*(angle*col_diffuse);
	ray.col = ray.col + attFactor *  specular*(specAngle*col_specular);


	ray.col.clamp();
}

