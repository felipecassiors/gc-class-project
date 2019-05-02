#include<iostream>
#include<Vector.cpp>

using namespace std;

class Shape
{
	Intersection* intersection;	
	intersection = new Intersection();
	 
	public:
		
		BSDF* bsdf;
		bsdf = new BSDF();
		Shape(){
		}
		Shape::Shape(BSDF bsdf){
			this->bsdf = bsdf;
		}
		Intersection Intersects(Ray ray){
			
			
			intersection.hit = false;
			intersection.distance = Infinity;
			intersection.object = null;
			for object in scene{
				Intersection temp = object.intersects(ray);
				
				if (temp.hit && temp.distance < intersection.distance){
				
					intersection = temp;
					intersection.object = object;
				}
			}
		return intersection;	

			
		}
		float calculateShaderGlobals(Ray ray, Intersection intersection ){
			
			
		}
		float surfaceArea{
			
			
		}
		
			
};

class BSDF{

	public:
		
		BSDFType type;
		Color3 color;
		BSDF(){
		}
		BSDF::BSDF(BSDFType type, Color3 color){
			this->type = type;
			this->color = color;
		}
};

class BSDFType{
	
	public:
		int Light = 0;
		int Diffuse = 1;
		int Specular = 2;
		int None = 3;
};

class Intersection{
	
	public:
		
		bool hit;
		float distance;
		int index;
		Intersection(){
		}
		Intersection::Intersection(bool hit, float distance, int index){
			this->hit = hit;
			this->distance = distance;
			this->index = index;
		}
		
};

class Scene{
	
	Shape* shape;
	shape = new Shape();
	
	public:
		
		vector<Shape> shapes;
		Scene(){
		}
		Scene::Scene(vector<Shapes> shapes){
			
			this->vector<Shapes> = shapes;
		}
		Intersection intersects(Ray ray){
			
		}
};

class Ray{
	
	public:
		
		Vector3 origin;
		Vector3 direction;
		Ray(){
			
		}
		Ray::Ray(Vector3 origin, Vector3 direction){
			
			this->origin = origin;
			this->direction = direction;
		}
		Vector3 point(float distance){
			
			return distance;
		}
};

class ShaderGlobals{
	
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
		
		ShaderGlobals(){
		}
		ShaderGlobals::ShaderGlobals(Vector3 point, Vector3 normal, Vector2 uv, Vector3 tangentU, Vector3 tangentV, Vector3 viewDirection, 
		Vector3 lightDirection, Vector3 lightPoint, Vector3 lightNormal){
			
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

class Sphere : public Shape{
	

	public:
		Vector3 position;
		float radius;
		
		Sphere(){
		}
		Sphere::Sphere(Vector3 position, float radius, BSDF bsdf){
			this->position = position;
			this->radius = radius;
			
		}
};

int main(int argc, char ** argv) {
	return 0;
}
