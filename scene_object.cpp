/***********************************************************
	
	Starter code for Assignment 3
	
	Implements scene_object.h

***********************************************************/

#include <cmath>
#include "scene_object.h"

bool UnitSquare::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	Vector3D direction;
	Point3D origin;

	direction = worldToModel * ray.dir;
	origin = worldToModel * ray.origin;

	double t = 0.0 - origin[2] / direction[2];

	if( t <= 0 ) {
		return(false);
	}

	Ray3D objSpaceRay = Ray3D(origin, direction);
	objSpaceRay.dir = direction;
	objSpaceRay.origin = origin;

	Point3D p (origin[0] + t * direction[0], origin[1] + t * direction[1], 0.0);
	if( p[0] >= -0.5 && p[0] <= 0.5 && p[1] >= -0.5 && p[1] <= 0.5 ) {
		if( ray.intersection.none || t < ray.intersection.t_value ) {
			Intersection intersect;
			//ray.intersection = intersect;
			ray.intersection.none = false;
			ray.intersection.t_value = t;
			ray.intersection.point = modelToWorld * p;
			ray.intersection.normal = transNorm(worldToModel, Vector3D(0, 0, 1));
			ray.intersection.normal.normalize();
			return(true);
		}
	}

	return(false);
}

bool UnitSphere::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	Vector3D dir;
	Point3D origin;

	dir = worldToModel * ray.dir;
	origin = worldToModel * ray.origin;

	Point3D sphereOrigin = Point3D(0, 0, 0);
	Ray3D objSpaceRay = Ray3D(origin, dir);
	objSpaceRay.dir = dir;
	objSpaceRay.origin = origin;
	Vector3D differentialOrigin = objSpaceRay.origin - sphereOrigin;

	double R = 1;
	double Rsquared = R * R;
	double A = dir.dot(dir);
	double B = 2 * dir.dot(differentialOrigin);
	double C = differentialOrigin.dot(differentialOrigin) - Rsquared;

	double denom = 2 * A;
	double descrim = B * B - 4 * A * C;
	if( descrim < 0 ) {
		return(false);
	}

	double t_1 = (-B + sqrt(descrim)) / denom;
	double t_2 = (-B - sqrt(descrim)) / denom;
	double intersect_t = t_1;
	if (t_2 < t_1){
		intersect_t = t_2;
	}
	// intersection must be in front of image plane
	if (intersect_t <= 0){
		return false;
	}

	if( ray.intersection.none || intersect_t < ray.intersection.t_value ) {
		Intersection intersect;
		//Point3D intersectPoint (origin + intersect_t * dir);
		double t (intersect_t);
		Point3D intersectPoint (origin[0] + t * dir[0], origin[1] + t * dir[1], origin[2] + t * dir[2]);

		//ray.intersection = intersect;
		ray.intersection.none = false;
		ray.intersection.t_value = intersect_t;
		ray.intersection.point = modelToWorld * intersectPoint;
		ray.intersection.normal = transNorm(worldToModel, Vector3D((origin + intersect_t * dir)[0],(origin + intersect_t * dir)[1],(origin + intersect_t * dir)[2]));
		ray.intersection.normal.normalize();
		return(true);
	}



	return(false);
}

void SceneNode::rotate(char axis, double angle) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
			this->trans = this->trans*rotation; 	
			angle = -angle;
		} 
		else {
			this->invtrans = rotation*this->invtrans; 
		}	
	}
}

void SceneNode::translate(Vector3D trans) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	this->trans = this->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	this->invtrans = translation*this->invtrans; 
}

void SceneNode::scale(Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	this->trans = this->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	this->invtrans = scale*this->invtrans; 
}


