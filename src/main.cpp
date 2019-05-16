#include <aurora/Math.h>
#include <aurora/Utility.h>
#include <aurora/Vector.h>
#include <aurora/Color.h>
#include <aurora/Matrix.h>
#include <aurora/Image.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>

using namespace aurora;
using namespace std;

float uniformRandom1D(){
	return rand()/(RAND_MAX + 1.0);
}

Vector2 uniformRandom2D(){
	return Vector2(uniformRandom1D(), uniformRandom1D());
}

void stratifiedSamples(int count, vector<Vector2> & samples){
	samples.resize(count);
	
	int size = sqrt(count);
	float inverseSize = 1.0 / size;
	
	for(int i = 0; i<count; i++){
		Vector2 offset(i/size, i%size);
		samples[i] = (offset + uniformRandom2D())* inverseSize;
	}
}

float gaussian1D(float sample, float width){
	float radius = width / 2;
	return fmax(0, exp(-sample * sample) - exp(-radius * radius));
}

float gaussian2D(Vector2 & sample, float width){
	return gaussian1D(sample.x, width) * gaussian1D(sample.y, width);
}


class Shape
{

public:
	Intersection *intersection;
	intersection = new Intersection();
	BSDF *bsdf;
	bsdf = new BSDF();
	Shape()
	{
	}
	Shape(BSDF *bsdf)
	{
		this->bsdf = bsdf;
	}
	virtual bool intersects(const Ray &ray, Intersection &intersection) = 0;
	virtual void calculateShaderGlobals(
		const Ray &ray, const Intersection &intersection,
		ShaderGlobals &shaderGlobals) = 0;
	virtual float surfaceArea() = 0;
};

class Film
{

public:
	float width;
	float height;
	Film()
	{
	}
	Film(float width, float height)
	{
		this->width = width;
		this->height = height;
	}
	float aspectRatio()
	{
		return width / height;
	}
};

class Camera
{

public:
	Film *film;
	float fieldOfView;
	Matrix4 worldMatrix;
	Camera()
	{
	}
	Camera(float fieldOfView, Film &film, Matrix4 worldOfMatrix)
	{
		this->fieldOfView = fieldOfView;
		this->film = film;
		this->worldOfMatrix = worldOfMatrix;
	}
	void lookAt(Vector3 &position, Vector3 &target, Vector3 &up)
	{
		Vector3 W = position - target;
		W = W / W.length();

		Vector3 U = up.cross(W);
		U = U / U.length();

		Vector3 V = W.cross(U);

		worldMatrix[0][0] = U.x;
		worldMatrix[0][1] = U.y;
		worldMatrix[0][2] = U.z;

		worldMatrix[1][0] = V.x;
		worldMatrix[1][1] = V.y;
		worldMatrix[1][2] = V.z;

		worldMatrix[2][0] = W.x;
		worldMatrix[2][1] = W.y;
		worldMatrix[2][2] = W.z;
	}

	Ray generateRay(float x, float y, Vector2 &sample) const
	{

		float scale = tan(fieldOfView / 2);

		Vector3 Pc = new Vector3();
		Pc.x = (2.0 * (x + sample.x + 0.5)/ film.width - 1.0-) * scale * film.aspectratio();
		Pc.y = (1.0 - 2.0 * (y+sample.y + 0.5) / film.height) * scale;
		Pc.z = -1.0;

		Vector3 Pl = Pc * worldMatrix;
		Vector3 P = new Vector3(worldMatrix[3][0], worldMatrix[3][1], worldMatrix[3][2]);
		Vector3 D = Pl - P;
		D = D / D.length();

		return new Ray(P, D);
	}
};

class BSDF
{

public:
	BSDFType type;
	Color3 color;
	BSDF()
	{
	}
	BSDF(BSDFType type, Color3 color)
	{
		this->type = type;
		this->color = color;
	}
};

enum BSDFType {

	public :
		int Light = 0;
	int Diffuse = 1;
	int Specular = 2;
	int None = 3;
};

class RenderOptions{
	public:
	
	int width;
	int height;
	int maximumDepth;
	int cameraSamples;
	int lightSamples;
	int diffuseSamples;
	float filterWidth;
	float gamma;
	float exposure;
	
	RenderOptions(){}
	RenderOptions(int width, int height, int maximumDepth, int cameraSamples, int lightSamples, int diffuseSamples, float filterWidth, float gamma, float exposure){
		this->width width;
		this->height = height;
		this->maximumDepth = maximumDepth;
		this->cameraSamples =  cameraSamples;
		this->lightSamples = lightSamples;
		this->diffuseSamples = diffuseSamples;
		this->filterWidth = filterWidth;
		this->gamma = gamma;
		this->exposure = exposure;
	}
	
};

class Renderer{
	
	public:
	RenderOptions * options;
	Camera * camera;
	Scene * scene;
	Renderer(){}
	Renderer(RenderOptions * options, Camera * camera, Scene * scene){
		this->options = options;
		this->camera = camera;
		this->scene = scene;
	}
	Color3 computeDirectIllumination(const BSDF & bsdf; shaderGlobals & shaderGlobals) const{
		
		return Color3();
	}
	Color3 computeIndirectIllumination(BSDF & bsdf; shaderGlobals & shaderGlobals, int depth){
	
		return Color3();
	}
	Color3 trace(const Ray & ray, int depth){
		
		Intersection intersection;
		if(scene->intersects(ray, intersection)) return Color3(1.0, 1.0, 1.0);
		
		return Color3;
		
		
	}
	void render(Image3 * image)const{
		
		const Vector2 half(0.5, 0.5);
		
		for(int i = 0; i < options->width; i++){
			for(int j = 0; j < options->height; j++){
				vector<Vector2> samples;
				stratifiedSamples(options->cameraSamples, samples);
				
				Color3 color;
				float weight = 0;
				
				for(int k=0; k<options->cameraSamples; k++){
					Vector2 sample = samples[k] - half) * options->filterWidth;
					Ray ray = camera->generateRay(i, j, sample);
					
					float w = gaussian2D(sample, options->filterWidth);
					
					color += trace(ray, 0) * w;
					weight += w;
				}
				
				color /= weight;
				color.applyExposure();
				color.applyGamma();
			}
		}
	}

};

class Intersection
{

public:
	bool hit;
	float distance;
	int index;
	Intersection()
	{
		hit = false;
		distance = AURORA_INFINITY;
		index = -1;
	}
	Intersection(bool hit, float distance, int index)
	{
		this->hit = hit;
		this->distance = distance;
		this->index = index;
	}
};

class Scene
{

	std::vector<Shape *> shapes;

	Scene();
	Scene(const std::vector<Shape *> &shapes)
	{
		this->shapes = shapes;
	}

public:
	vector<Shape> shapes;
	Scene()
	{
	}
	Scene::Scene(vector<Shapes> shapes)
	{

		this->vector<Shapes> = shapes;
	}
	bool intersects(const Ray &ray, Intersection &intersection)
	{
		for (int i = 0; i < shapes.size(); i++)
		{
			Shape *shape = shapes[i];

			Intersection temp;
			shape->intersects(ray, temp);

			if (temp.hit && temp.distance < intersection.distance)
			{
				intersection.hit = temp.hit;
				intersection.distance = temp.distance;
				intersection.index = i;
			}
		}

		return intersection.hit;
	}
};

class Ray
{

public:
	Vector3 origin;
	Vector3 direction;
	Ray()
	{
	}
	Ray(const Vector3 &origin, const Vector3 &direction)
	{

		this->origin = origin;
		this->direction = direction;
	}
	Vector3 point(float distance)
	{

		return origin + direction * distance;
	}
};

class ShaderGlobals
{

public:
	Vector3 point;
	Vector3 normal;
	Vector2 uv;
	Vector3 tangentU;
	Vector3 tangentV;
	Vector3 viewDirection;
	Vector3 lightDirection;
	Vector3 lightPoint;
	Vector3 lightNormal;

	ShaderGlobals()
	{
	}
	ShaderGlobals(const Vector3 &point, const Vector3 &normal, const Vector2 &uv, const Vector3 &tangentU,
				  const Vector3 &tangentV, const Vector3 &viewDirection,
				  const Vector3 &lightDirection, const Vector3 &lightPoint, const Vector3 &lightNormal)
	{

		this->point = point;
		this->normal = normal;
		this->uv = uv;
		this->tangentU = tangentU;
		this->tangentV = tangentV;
		this->viewDirection = viewDirection;
		this->lightDirection = lightDirection;
		this->lightPoint = lightPoint;
		this->lightNormal = lightNormal;
	}
};

class Sphere : public Shape
{

public:
	Vector3 position;
	float radius;

	Sphere()
	{
	}
	Sphere(Vector3 position, float radius, BSDF bsdf)
	{
		this->position = position;
		this->radius = radius;
	}

	virtual bool intersects(const Ray &ray, Intersection &intersection)
	{

		Vector3 l = position - ray.origin;
		float t = l.dot(ray.direction);

		if (t < 0)
			return false;

		float d2 = l.dot(l) - t * t;
		float r2 = radius * radius;

		if (d2 > r2)
			return false;

		float dt = std::sqrt(r2 - d2);

		float t0 = t - dt;
		float t1 = t + dt;

		if (t0 > t1)
			std::swap(t0, t1);

		if (t0 < 0)
		{
			t0 = t1;

			if (t0 < 0)
				return false;
		}

		intersection.hit = true;
		intersection.distance = t0;

		return true;
	}
	virtual void calculateShaderGlobals(
		const Ray &ray, const Intersection &intersection,
		ShaderGlobals &shaderGlobals)
	{
		shaderGlobals.point = ray.point(intersection.distance);
		shaderGlobals.normal = (shaderGlobals.point - position).normalize();

		float theta = std::atan2(shaderGlobals.normal.x, shaderGlobals.normal.z);
		float phi = std::acos(shaderGlobals.normal.y);

		shaderGlobals.uv.x = theta * AURORA_INV_PI * 0.5;
		shaderGlobals.uv.y = phi * AURORA_INV_PI;

		shaderGlobals.tangentU.x = std::cos(theta);
		shaderGlobals.tangentU.y = 0;
		shaderGlobals.tangentU.z = -std::sin(theta);

		shaderGlobals.tangentV.x = std::sin(theta) * std::cos(phi);
		shaderGlobals.tangentV.y = -std::sin(phi);
		shaderGlobals.tangentV.z = std::cos(theta) * std::cos(phi);
	}
	virtual float surfaceArea()
	{
		return 4.0 * AURORA_PI * radius * radius;
	}
};

int main(int argc, char **argv)
{
	RenderOptions options(500, 250, 1, 4, 1, 1, 2.0, 2.2, 0);
	
		
	Renderer Renderer(&options, &camera, &scene);
	BSDF *bsdf = new BSDF(BSDFType::Diffuse, Color3(1.0, 1.0, 1.0));
	Shape *shape = new Sphere(Vector3(0, 0, 0), 1.0, bsdf);

	std::vector<Shape *> shapes;
	shapes.push_back(shape);

	Scene scene(shapes);

	Ray ray(Vector3(0, 0, 10.0), Vector3(0, 0, -1.0));
	Intersection intersection;

	scene.intersects(ray, intersection);

	std::cout << "Hit: " << intersection.hit << std::endl;
	std::cout << "Distance: " << intersection.distance << std::endl;
	std::cout << "Index: " << intersection.index << std::endl;

	delete bsdf;
	delete shape;

	return 0;
}
