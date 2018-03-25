/***********************************************************
	
	Starter code for Assignment 3

	Implementations of functions in raytracer.h, 
	and the main function which specifies the scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

void Raytracer::traverseScene(Scene& scene, Ray3D& ray)  {
	for (size_t i = 0; i < scene.size(); ++i) {
		SceneNode* node = scene[i];

		if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
			// if we intersect with something, we want secondary reflection
			ray.intersection.mat = node->mat;
		}
	}
}

void Raytracer::computeTransforms(Scene& scene) {
	// right now this method might seem redundant. But if you decide to implement 
	// scene graph this is where you would propagate transformations to child nodes
		
	for (size_t i = 0; i < scene.size(); ++i) {
		SceneNode* node = scene[i];

		node->modelToWorld = node->trans;
		node->worldToModel = node->invtrans; 
	}
}

void Raytracer::computeShading(Ray3D& ray, LightList& light_list, Scene& scene) {
	for (size_t  i = 0; i < light_list.size(); ++i) {
		LightSource* light = light_list[i];
		
		int in_shadow = 0;
		// Each lightSource provides its own shading function.
		// Implement shadows here if needed.
		// TODO: hard shadowing
		// look from ray intersection point to light source, if we intersect
		// an object before hitting light source then shadow
		Point3D origin = ray.intersection.point;
		Vector3D to_light = light->get_position() - origin;
		to_light.normalize();
		// offset so we don't intersect the same object
		Ray3D shadow(origin+0.000000001*to_light, to_light);
		traverseScene(scene, shadow);
		// if we hit something before the light
		if (!shadow.intersection.none){
			ray.col = Color(0.0, 0.0, 0.0);
			return;
		}
		light->shade(ray);
	}
}

Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list, int depth) {
	Color col(0.0, 0.0, 0.0); 
	if (depth == 0){
		return col;
	}
	traverseScene(scene, ray); 

	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray, light_list, scene); 
		col = ray.col;  
		Ray3D reflected_ray;
		Vector3D normal = ray.intersection.normal;
		Vector3D dir = -ray.dir;
		dir.normalize();
		// add offset so it doesn't intersect the same object again
		float EPS = 0.000000001;
		reflected_ray.origin = ray.intersection.point+EPS*dir;
		reflected_ray.dir = 2*(dir.dot(normal))*normal-dir;
		reflected_ray.dir.normalize();

		col = col+0.25*shadeRay(reflected_ray, scene, light_list, depth-1);
	}
	col.clamp();
	return col; 
}	

void Raytracer::render(Camera& camera, Scene& scene, LightList& light_list, Image& image) {
	computeTransforms(scene);

	Matrix4x4 viewToWorld;
	double factor = (double(image.height)/2)/tan(camera.fov*M_PI/360.0);

	viewToWorld = camera.initInvViewMatrix();

	// Construct a ray for each pixel.
	for (int i = 0; i < image.height; i++) {
		for (int j = 0; j < image.width; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
			imagePlane[0] = (-double(image.width)/2 + 0.5 + j)/factor;
			imagePlane[1] = (-double(image.height)/2 + 0.5 + i)/factor;
			imagePlane[2] = -1;

			
			
			Ray3D ray;
		    ray.origin = viewToWorld * origin;
            ray.dir = viewToWorld * (imagePlane - origin);
            ray.dir.normalize();

			Color col = shadeRay(ray, scene, light_list, 3); 
			image.setColorAtPixel(i, j, col);			
		}
	}
}

