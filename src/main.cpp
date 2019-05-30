#include <aurora/Math.h>
#include <aurora/Utility.h>
#include <aurora/Vector.h>
#include <aurora/Color.h>
#include <aurora/Matrix.h>
#include <aurora/Image.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>

using namespace aurora;

float uniformRandom1D() {
    return std::rand() / (RAND_MAX + 1.0);
}
Vector2 uniformRandom2D() {
    return Vector2(uniformRandom1D(), uniformRandom1D());
}

void stratifiedSample(int count, std::vector<Vector2> & samples) {
    samples.reserve(count);
    
    int size = std::sqrt(count);
    float inverseSize = 1.0 / size;
    
    for (int i = 0; i < count; i++) {
        Vector2 offset(i / size, i % size);
        Vector2 point = (offset + uniformRandom2D()) * inverseSize;
        
        samples.push_back(point);
    }
}

float gaussian1D(float sample, float width) {
    float radius = width * 0.5;
    return std::fmax(0, std::exp(-sample * sample) - std::exp(-radius * radius));
}
float gaussian2D(const Vector2 & sample, float width) {
    return gaussian1D(sample.x, width) * gaussian1D(sample.y, width);
}

struct Ray {
    Vector3 origin;
    Vector3 direction;
    
    Ray() {}
    Ray(const Vector3 & origin, const Vector3 & direction) {
        this->origin = origin;
        this->direction = direction;
    }
    
    Vector3 point(float distance) const {
        return origin + direction * distance;
    }
};

struct Intersection {
    bool hit;
    float distance;
    int index;
    
    Intersection() {
        hit = false;
        distance = AURORA_INFINITY;
        index = -1;
    }
    Intersection(bool hit, float distance, int index) {
        this->hit = hit;
        this->distance = distance;
        this->index = index;
    }
};

struct ShaderGlobals {
    Vector3 point;
    Vector3 normal;
    Vector2 uv;
    Vector3 tangentU;
    Vector3 tangentV;
    Vector3 viewDirection;
    Vector3 lightDirection;
    Vector3 lightPoint;
    Vector3 lightNormal;
    
    ShaderGlobals() {}
    ShaderGlobals(
            const Vector3 & point, const Vector3 & normal, const Vector2 & uv,
            const Vector3 & tangentU, const Vector3 & tangentV,
            const Vector3 & viewDirection, const Vector3 & lightDirection,
            const Vector3 & lightPoint, const Vector3 & lightNormal) {
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

enum BSDFType {
    Light = 0,
    Diffuse,
    Specular,
    None
};

struct BSDF {
    BSDFType type;
    Color3 color;
    
    BSDF() {}
    BSDF(BSDFType type, const Color3 & color) {
        this->type = type;
        this->color = color;
    }
};

struct Shape {
    BSDF * bsdf;
    
    Shape() {}
    Shape(BSDF * bsdf) {
        this->bsdf = bsdf;
    }
    
    virtual bool intersects(const Ray & ray, Intersection & intersection) const = 0;
    virtual void calculateShaderGlobals(
            const Ray & ray, const Intersection & intersection,
            ShaderGlobals & shaderGlobals) const = 0;
    virtual float surfaceArea() const = 0;
};

struct Sphere : Shape {
    Vector3 position;
    float radius;
    
    Sphere() : Shape() {}
    Sphere(const Vector3 & position, float radius, BSDF * bsdf) : Shape(bsdf) {
        this->position = position;
        this->radius = radius;
    }
    
    virtual bool intersects(const Ray & ray, Intersection & intersection) const {
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
        
        if (t0 < 0) {
            t0 = t1;
            
            if (t0 < 0)
                return false;
        }
        
        intersection.hit = true;
        intersection.distance = t0;
        
        return true;
    }
    virtual void calculateShaderGlobals(
            const Ray & ray, const Intersection & intersection,
            ShaderGlobals & shaderGlobals) const {
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
        
        shaderGlobals.viewDirection = -ray.direction;
    }
    virtual float surfaceArea() const {
        return 4.0 * AURORA_PI * radius * radius;
    }
};

struct Scene {
    std::vector<Shape *> shapes;
    
    Scene();
    Scene(const std::vector<Shape *> & shapes) {
        this->shapes = shapes;
    }
    
    bool intersects(const Ray & ray, Intersection & intersection) const {
        for (int i = 0; i < shapes.size(); i++) {
            Shape * shape = shapes[i];
            
            Intersection temp;
            shape->intersects(ray, temp);
            
            if (temp.hit && temp.distance < intersection.distance) {
                intersection.hit = temp.hit;
                intersection.distance = temp.distance;
                intersection.index = i;
            }
        }
        
        return intersection.hit;
    }
};

struct Film {
    float width;
    float height;
    
    Film() {}
    Film(float width, float height) {
        this->width = width;
        this->height = height;
    }
    
    float aspectRatio() const {
        return width / height;
    }
};

struct Camera {
    float fieldOfView;
    Film film;
    Matrix4 worldMatrix;
    
    Camera() {}
    Camera(float fieldOfView, const Film & film, const Matrix4 & worldMatrix) {
        this->fieldOfView =fieldOfView;
        this->film = film;
        this->worldMatrix = worldMatrix;
    }
    
    void lookAt(const Vector3 & position, const Vector3 & target, const Vector3 & up) {
        Vector3 w = (position - target).normalize();
        Vector3 u = up.cross(w).normalize();
        Vector3 v = w.cross(u);
        
        worldMatrix[0][0] = u.x;
        worldMatrix[0][1] = u.y;
        worldMatrix[0][2] = u.z;
        worldMatrix[0][3] = 0;
        
        worldMatrix[1][0] = v.x;
        worldMatrix[1][1] = v.y;
        worldMatrix[1][2] = v.z;
        worldMatrix[1][3] = 0;
        
        worldMatrix[2][0] = w.x;
        worldMatrix[2][1] = w.y;
        worldMatrix[2][2] = w.z;
        worldMatrix[2][3] = 0;
        
        worldMatrix[3][0] = position.x;
        worldMatrix[3][1] = position.y;
        worldMatrix[3][2] = position.z;
        worldMatrix[3][3] = 1.0;
    }
    Ray generateRay(float x, float y, const Vector2 & sample) const {
        float scale = std::tan(fieldOfView * 0.5);
        
        Vector3 pixel;
        
        pixel.x = (2.0 * (x + sample.x + 0.5) / film.width - 1.0) * scale * film.aspectRatio();
        pixel.y = (1.0 - 2.0 * (y + sample.y + 0.5) / film.height) * scale;
        pixel.z = -1.0;
        
        pixel *= worldMatrix;
        
        Vector3 position(worldMatrix[3][0], worldMatrix[3][1], worldMatrix[3][2]);
        Vector3 direction = (pixel - position).normalize();
        
        return Ray(position, direction);
    }
};

struct RenderOptions {
    int width;
    int height;
    int maximumDepth;
    int cameraSamples;
    int lightSamples;
    int diffuseSamples;
    float filterWidth;
    float gamma;
    float exposure;
    
    RenderOptions() {}
    RenderOptions(int width, int height, int maximumDepth,
            int cameraSamples, int lightSamples, int diffuseSamples,
            float filterWidth, float gamma, float exposure) {
        this->width = width;
        this->height = height;
        this->maximumDepth = maximumDepth;
        this->cameraSamples = cameraSamples;
        this->lightSamples = lightSamples;
        this->diffuseSamples = diffuseSamples;
        this->filterWidth = filterWidth;
        this->gamma = gamma;
        this->exposure = exposure;
    }
};

struct Renderer {
    RenderOptions * options;
    Camera * camera;
    Scene * scene;
    
    Renderer() {}
    Renderer(RenderOptions * options, Camera * camera, Scene * scene) {
        this->options = options;
        this->camera = camera;
        this->scene = scene;
    }
    
    Color3 computeDirectIllumination(
            const BSDF & bsdf, ShaderGlobals & shaderGlobals) const {
        return Color3();
    }
    Color3 computeIndirectIllumination(
            const BSDF & bsdf, ShaderGlobals & shaderGlobals) const {
        return Color3();
    }
    
    Color3 trace(const Ray & ray, int depth) const {
        Intersection intersection;
        
       	if (scene->intersects(ray, intersection)) {
            const Shape * shape = scene->shapes[intersection.index];
            const BSDF * bsdf = shape->bsdf;
            
            ShaderGlobals shaderGlobals;
            shape->calculateShaderGlobals(ray, intersection, shaderGlobals);
            
            float cosTheta = shaderGlobals.viewDirection.dot(shaderGlobals.normal);
            
            return bsdf->color * cosTheta;
        }
        
        return Color3();
    }
    Color3 render(Image3 * image) const {
        const Vector2 half(0.5, 0.5);
        
        for (int i = 0; i < options->width; i++) {
            for (int j = 0; j < options->height; j++) {
                std::vector<Vector2> samples;
                stratifiedSample(options->cameraSamples, samples);
                
                Color3 color;
                float weight = 0;
                
                for (int k = 0; k < options->cameraSamples; k++) {
                    Vector2 sample = (samples[k] - half) * options->filterWidth;
                    Ray ray = camera->generateRay(i, j, sample);
                    
                    float w = gaussian2D(sample, options->filterWidth);
                    
                    color += trace(ray, 0) * w;
                    weight += w;
                }
                
                color /= weight;
                
                color.applyExposure(options->exposure);
                color.applyGamma(options->gamma);
                
                image->setPixel(i, j, color);
            }
        }
    }
};

int main(int argc, char **argv)
{
	RenderOptions options(500, 250, 1, 16, 1, 1, 2.0, 2.2, 0);
	
	Film film(options.width, options.height);
	
	Camera camera(radians(20.0), film, Matrix4());
	camera.lookAt(Vector3(0, 0, 35.0), Vector3(0, 0, 0), Vector3(0, 1.0, 0));
	
	BSDF * whiteDiffuse = new BSDF(BSDFType::Diffuse, Color3(1.0, 1.0, 1.0));
	BSDF * redDiffuse = new BSDF(BSDFType::Diffuse, Color3(1.0, 0, 0));	
	BSDF * greenDiffuse = new BSDF(BSDFType::Diffuse, Color3(0, 1.0, 0));		
	BSDF * lightMaterial = new BSDF(BSDFType::Diffuse, Color3(1.0, 1.0, 1.0));		
	
	Shape * left = new Sphere(Vector3(-1.0e5 - 5.0, 0, 0), 1.0e5, redDiffuse);
	Shape * right = new Sphere(Vector3(1.0e5 + 5.0, 0, 0), 1.0e5, greenDiffuse);
	Shape * bottom = new Sphere(Vector3(0, -1.0e5 -5.0, 0), 1.0e5, whiteDiffuse);	
	Shape * top = new Sphere(Vector3(0, 1.0e5 + 5.0, 0), 1.0e5, whiteDiffuse);
	Shape * back = new Sphere(Vector3(0, 0, -1.0e5 -5.0), 1.0e5, whiteDiffuse);
	Shape * frontSphere = new Sphere(Vector3(2.0, -3.0, 2.0), 2.0, whiteDiffuse);
	Shape * backSphere = new Sphere(Vector3(-2.0, -3.0, -2.0), 2.0, whiteDiffuse);
	Shape * light = new Sphere(Vector3(0, 3.0, 0), 0.5, lightMaterial);
		
	std::vector<Shape *> shapes;
	shapes.push_back(left);
	shapes.push_back(right);
	shapes.push_back(bottom);
	shapes.push_back(top);
	shapes.push_back(back);
	shapes.push_back(frontSphere);
	shapes.push_back(backSphere);
	shapes.push_back(light);
	
	Scene scene(shapes);
		
	Renderer renderer(&options, &camera, &scene);
	Image3 * image = new Image3(options.width, options.height);
	renderer.render(image);
	
	if (writeImage("output.ppm", image))
        std::cout << "Success." << std::endl;
    else
        std::cout << "Failure." << std::endl;
    
    delete whiteDiffuse;
    delete redDiffuse;
    delete greenDiffuse;
    delete lightMaterial;
    
    delete image;
    
    for (Shape * shape : shapes)
        delete shape;
    
    return 0;
}
