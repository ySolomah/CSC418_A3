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

void Raytracer::computeShading(Scene& scene, Ray3D& ray, LightList& light_list) {
	for (size_t  i = 0; i < light_list.size(); ++i) {
		LightSource* light = light_list[i];
		
		// Each lightSource provides its own shading function.
		// Implement shadows here if needed.
		Point3D lightPos = light->get_position();
		Point3D shadePos = ray.intersection.point;
		Vector3D dir = lightPos-shadePos;
		Ray3D intersectShadowRay(shadePos, lightPos-shadePos);
		traverseScene(scene, intersectShadowRay);
		if(intersectShadowRay.intersection.none) {
			light->shade(ray);
		}
	}
}

Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list, int reflDepth) {
	Color col(0.0, 0.0, 0.0); 
	traverseScene(scene, ray); 
	int maxDepth = 4;

	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		if(reflDepth < maxDepth) {
			computeShading(scene, ray, light_list); 
			Vector3D direction = ray.dir;
			Vector3D norm = ray.intersection.normal;
			Vector3D reflectionDirection = -2*direction.dot(norm) * norm + direction;
			Ray3D reflectionRay(ray.intersection.point, reflectionDirection);
			reflectionRay.origin = ray.intersection.point;
			reflectionRay.dir = reflectionDirection;
			Color refCol = shadeRay(reflectionRay, scene, light_list, reflDepth + 1);
			col = ray.col + ray.intersection.mat->reflection_amnt * refCol;
		}
	}

	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.  
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
			Color finalcol (0.0, 0.0, 0.0);
			for(float m = -0.25; m < 0.5; m = m+0.25) {
				for(float n = -0.25; n < 0.5; n = n+0.25) {
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					imagePlane[0] = (-double(image.width)/2 + 0.5 + j + m)/factor;
					imagePlane[1] = (-double(image.height)/2 + 0.5 + i + n)/factor;
					imagePlane[2] = -1;

					
					
					Ray3D ray;
					// TODO: Convert ray to world space  

					Vector3D direction = imagePlane - origin;
					direction = viewToWorld * direction;
					origin = viewToWorld * origin;
					ray = Ray3D(origin, direction);
					ray.origin = origin;
					ray.dir = direction;
					
					Color col = shadeRay(ray, scene, light_list, 0); 
					finalcol = finalcol + 0.11111 * col;
				}
			}
			image.setColorAtPixel(i, j, finalcol);			
		}
	}
}

