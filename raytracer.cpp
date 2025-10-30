#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <windows.h>
#include <chrono>

typedef unsigned char RGB[3];



inline float determinant (const parser::Vec3f &a, const parser::Vec3f &b, const parser::Vec3f &c) {
    return a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z) + a.z * (b.x * c.y - b.y * c.x);
}




void intersectsTriangle(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Triangle &t, parser::HitRecord& h, const parser::Scene& s, bool BFC){


    const parser::Vec3f &a = s.vertex_data[t.indices.v0_id-1];
    const parser::Vec3f &b = s.vertex_data[t.indices.v1_id-1];
    const parser::Vec3f &c = s.vertex_data[t.indices.v2_id-1];

    parser::Vec3f ab = b - a;
    parser::Vec3f ac = c - a;
    parser::Vec3f ao = o - a;

    parser::Vec3f normal = ab ^ ac;

    if(BFC && normal * d > 0){
        return;
    }

    float detA = determinant(ab,ac,d);

    if (fabs(detA) < 1e-8) return;

    float beta = determinant(ao,ac,d)/detA;
    float gamma = determinant(ab,ao,d)/detA;
    float tRay = -determinant(ab,ac,ao)/detA;

    if( beta >= 0 && gamma >= 0 && beta + gamma <= 1 && tRay < h.t && tRay > 1e-8f ){
        h.hit = true;
        h.normal = normal.normalized();
        h.t = tRay;
        h.point = o + d * tRay;
        h.Mid = t.material_id;

    } 
    

}

void intersectsSphere(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Sphere &sp, parser::HitRecord& h, const parser::Scene& s){

    const parser::Vec3f &c = s.vertex_data[sp.center_vertex_id - 1];
    
    parser::Vec3f oc = o - c;

    float A = d * d;
    float B = d * (oc) * 2;
    float C = (oc) * (oc) - (sp.radius * sp.radius);

    float disc = B * B - 4 * A * C;

    if(disc < 0.0f)    return;
    

    float sqrt = std::sqrt(disc);
    float t1 =  (-B + sqrt) / (2.0f * A);
    float t2 =  (-B - sqrt) / (2.0f * A);
    float t = std::min(t1,t2);
    
    if (t < 0) return;
    
    if( h.t > t){
        h.hit = true;
        h.t = t;  
        h.point = o + d * t;
        h.normal = (h.point - c).normalized();
        if(h.normal * d > 0) h.normal = h.normal * -1;
        h.Mid = sp.material_id;
    }
    

}

void intersectsPlane(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Plane &p, parser::HitRecord& h, const parser::Scene& s){

    parser::Vec3f n = p.normal.normalized();
    const parser::Vec3f &a = s.vertex_data[p.center_vertex_id - 1];
    
    float dn = d * n;

    if (fabs(dn) < 1e-8f) return;

    float t = (a - o) * n / dn ;

    if (t < 0) return;

    if(t < h.t){
        h.hit = true;
        h.t = t;
        if (d * n > 0.0f)
            n = n * -1.0f;
        h.normal = n;
        h.point = o + d * t;
        h.Mid = p.material_id;
    }

     

}

void intersectsCylinder(const parser::Vec3f &o, const parser::Vec3f &d
, const parser::Cylinder &cyl, parser::HitRecord& h, const parser::Scene& s)
{
    parser::Vec3f C = s.vertex_data[cyl.center_vertex_id - 1];
    parser::Vec3f v = cyl.axis.normalized();  

    
    parser::Vec3f delta = o - C;
    float dv = d*v;
    float delta_v = delta * v;

    parser::Vec3f d_perp = d - v * dv;
    parser::Vec3f delta_perp = delta - (v * delta_v);

    float A = (d_perp * d_perp);
    float B = 2.0f * (d_perp * delta_perp);
    float Cc = (delta_perp * delta_perp) - cyl.radius * cyl.radius;

    float disc = B * B - 4.0f * A * Cc;
    if (disc < 0.0f) return; 

    float sqrt_disc = sqrt(disc);
    float t1 = (-B - sqrt_disc) / (2.0f * A);
    float t2 = (-B + sqrt_disc) / (2.0f * A);

    float t_cyl = std::numeric_limits<float>::infinity();

    if (t1 > s.shadow_ray_epsilon) t_cyl = t1;
    else if (t2 > s.shadow_ray_epsilon) t_cyl = t2;

    if (t_cyl < std::numeric_limits<float>::infinity()) {
        parser::Vec3f P = o + d * t_cyl;
        float proj = (P - C) * v;
        float halfH = cyl.height * 0.5f;

        
        if (proj >= -halfH && proj <= halfH) {
            if (t_cyl < h.t) {
                h.hit = true;
                h.t = t_cyl;
                h.point = P;
                parser::Vec3f radial = (P - C) - (v * proj);
                h.normal = radial.normalized();
                
                h.Mid = cyl.material_id;
            }
        }
    }

    
    float halfH = cyl.height * 0.5f;
    parser::Vec3f capTopCenter = C + v * halfH;
    parser::Vec3f capBotCenter = C - v * halfH;

    
    {
        float denom = d * v;
        if (fabs(denom) > 1e-6f) {
            float t_cap = (capTopCenter - o ) *  v / denom;
            if (t_cap > s.shadow_ray_epsilon && t_cap < h.t) {
                parser::Vec3f P = o + d * t_cap;
                if ((P - capTopCenter).len() <= cyl.radius) {
                    h.hit = true;
                    h.t = t_cap;
                    h.point = P;
                    h.normal = v;
                    h.Mid = cyl.material_id;

                }
            }
        }
    }

    
    {
        float denom = d * v;
        if (fabs(denom) > 1e-6f) {
            float t_cap = (capBotCenter - o) * v / denom;
            if (t_cap > s.shadow_ray_epsilon && t_cap < h.t) {
                parser::Vec3f P = o + d * t_cap;
                if ((P - capBotCenter).len() <= cyl.radius) {
                    h.hit = true;
                    h.t = t_cap;
                    h.point = P;
                    h.normal = v * -1;
                    
                    h.Mid = cyl.material_id;
                }
            }
        }
    }
}



void intersectsMesh(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Mesh &m, parser::HitRecord& h, const parser::Scene& s,bool BFC){

    for(const parser::Face &f : m.faces){

        const parser::Vec3f &a = s.vertex_data[f.v0_id-1];
        const parser::Vec3f &b = s.vertex_data[f.v1_id-1];
        const parser::Vec3f &c = s.vertex_data[f.v2_id-1];

        parser::Vec3f ab = b - a;
        parser::Vec3f ac = c - a;
        parser::Vec3f ao = o - a;

        parser::Vec3f normal = ab ^ ac;

        if(BFC && normal * d > 0){
            continue;   
        }

        float detA = determinant(ab,ac,d);

        if (fabs(detA) < 1e-8) continue;

        float beta = determinant(ao,ac,d)/detA;
        float gamma = determinant(ab,ao,d)/detA;
        float tRay = -determinant(ab,ac,ao)/detA;

        if( beta >= 0 && gamma >= 0 && beta + gamma <= 1 && tRay < h.t && tRay > 1e-8f ){
            h.hit = true;
            h.normal = normal.normalized();
            h.t = tRay;
            h.point = o + d * tRay;
            h.Mid = m.material_id;
        } 
    
    }

}




parser::Vec3f getPixelRay(const parser::Camera &c, const int i, const int j){

    parser::Vec3f e = c.position;
    parser::Vec3f v = c.up;
    parser::Vec3f w = c.gaze * -1;
    parser::Vec3f u = v ^ w;


    float l = c.near_plane.x;
    float r = c.near_plane.y;
    float b = c.near_plane.z;
    float t = c.near_plane.w;

    parser::Vec3f m = e + w * c.near_distance * -1;
    parser::Vec3f q = m + u * l + v * t;

    float su = (i+0.5) * (r-l) / c.image_width;
    float sv = (j+0.5) * (t-b) / c.image_height;
    
    parser::Vec3f s = q + u * su - v * sv;

    return (s - e).normalized();


}


void renderRay(const parser::Scene &scene, const parser::Ray &r, parser::HitRecord & closestHit, bool BFC){
        
    for(const parser::Triangle &t : scene.triangles) {

        intersectsTriangle(r.origin,r.direction,t,closestHit,scene,BFC);

    }
    for(const parser::Sphere &s : scene.spheres) {

        intersectsSphere(r.origin,r.direction,s,closestHit,scene);
            
    }
    for(const parser::Plane &p : scene.planes) {
        
        intersectsPlane(r.origin,r.direction,p,closestHit,scene);
            
    }    

    for(const parser::Cylinder &c : scene.cylinders){

        intersectsCylinder(r.origin,r.direction,c,closestHit,scene);

    }

    for(const parser::Mesh &m : scene.meshes){

        intersectsMesh(r.origin,r.direction,m,closestHit,scene,BFC);

    }

}

bool inShadow(const parser::Scene &scene,const parser::Vec3f &offSet,const parser::Vec3f &ShadowDirection, const float t){

    parser::HitRecord hr;
    hr.t = t;

    for(const parser::Triangle &t : scene.triangles) {

        intersectsTriangle(offSet,ShadowDirection,t,hr,scene,false);
        if(hr.hit && hr.t >= 0.0f) return true;

    }
    
    for(const parser::Sphere &s : scene.spheres) {

        intersectsSphere(offSet,ShadowDirection,s,hr,scene);
        if(hr.hit) return true;
            
    }
    for(const parser::Plane &p : scene.planes) {
        
        intersectsPlane(offSet,ShadowDirection,p,hr,scene);
        if(hr.hit) return true;
            
    }    

    for(const parser::Cylinder &c : scene.cylinders){

        intersectsCylinder(offSet,ShadowDirection,c,hr,scene);
        if(hr.hit) return true;

    }

    for(const parser::Mesh &m : scene.meshes){

        intersectsMesh(offSet,ShadowDirection,m,hr,scene,false);
        if(hr.hit && hr.t >= 1e-8f) return true;

    }
    
    return false;
    
}

parser::Ray Reflect(const parser::HitRecord &hitrecord, const parser::Scene &s, const parser::Ray &r ) {
    parser::Ray incomingRay = r;
    parser::Ray refRay;

    
    parser::Vec3f wo = incomingRay.direction.normalized();
    parser::Vec3f normal = hitrecord.normal.normalized();
   
    parser::Vec3f wr = wo - normal * (normal * wo * 2.0f);

    refRay.direction = wr.normalized();
    refRay.origin = hitrecord.point + (normal * 0.001);
    
    return refRay;
}

parser::Vec3i GetColor(const parser::HitRecord &hitrecord, const parser::Scene &s, const parser::Ray &r, int depth){

    
    if ( depth > s.max_recursion_depth  ){
        return parser::Vec3i {0,0,0};
    }

    float color_r = 0;
    float color_g = 0;
    float color_b = 0;

    parser::Material material = s.materials[hitrecord.Mid -1];


    for(parser::PointLight light : s.point_lights){

        
        parser::Vec3f l = (light.position - hitrecord.point).normalized();
        parser::Vec3f v = r.direction.normalized() * -1;
        //check shadow

        parser::Vec3f antiCollision = hitrecord.point + (hitrecord.normal * s.shadow_ray_epsilon);
        parser::Vec3f sV = (light.position - antiCollision);
        float t = sV.len();
        sV = sV.normalized();
        

        
        if(inShadow(s,antiCollision,sV,t)){
            continue;
        }

        

        //diffuse

        float dSquare = (light.position - hitrecord.point).len() * (light.position - hitrecord.point).len();

        float cosT =  l * hitrecord.normal;

        if(cosT < 0) cosT = 0;
        
        color_r += material.diffuse.x * light.intensity.x * cosT / dSquare;
        color_g += material.diffuse.y * light.intensity.y * cosT / dSquare;
        color_b += material.diffuse.z * light.intensity.z * cosT / dSquare; 

        //specular

        parser::Vec3f h = (l + v) / ((l + v).len());

        float cosA = h * hitrecord.normal;

        if(cosA < 0 ) cosA = 0;

        float phong = pow(cosA, material.phong_exponent);

        color_r += material.specular.x * light.intensity.x * phong / dSquare;
        color_g += material.specular.y * light.intensity.y * phong / dSquare;
        color_b += material.specular.z * light.intensity.z * phong / dSquare;

    }

    //ambient 

        color_r += material.ambient.x * s.ambient_light.x;
        color_g += material.ambient.y * s.ambient_light.y;
        color_b += material.ambient.z * s.ambient_light.z;


    //refl

        if(material.is_mirror)
        {
            
            parser::Ray reflection = Reflect (hitrecord,s,r);
            
            parser::HitRecord reflectionRecord;

            renderRay(s,reflection,reflectionRecord,false);

            if (reflectionRecord.hit)
            {
                parser::Vec3i refcolor = GetColor(reflectionRecord,s,reflection,depth +1);
                color_r += material.mirror.x * refcolor.x;
                color_g += material.mirror.y * refcolor.y;
                color_b += material.mirror.z * refcolor.z; 
            }        


        }

    if(color_r > 255) color_r = 255;
    if(color_g > 255) color_g = 255;
    if(color_b > 255) color_b = 255;

    int icolor_r = (int) round(color_r);
    int icolor_g = (int) round(color_g);
    int icolor_b = (int) round(color_b);
    
    return parser::Vec3i{icolor_r, icolor_g, icolor_b};

}




int main(int argc, char* argv[])
{
    
    
    auto start = std::chrono::high_resolution_clock::now();

    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    int cameraId = 1;

    
    std::cout << scene.max_recursion_depth;
    
    
   


    for(const parser::Camera& c : scene.cameras){


        std::string filename = "ouroutputs/" + std::string(argv[1]) + "_camera" + std::to_string(cameraId++) + ".ppm";

        int width = c.image_width;
        int height = c.image_height;

        unsigned char* image = new unsigned char [width * height * 3];

        
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                

                parser::Ray r = parser::Ray(c.position, getPixelRay(c, x, y).normalized());

                parser::HitRecord closestHit;

                renderRay(scene,r,closestHit,true);
                
        

                int pixelIndex = (y * width + x) * 3; 

                if(closestHit.hit)
                        {                    
                            parser::Vec3i color = GetColor(closestHit, scene, r, 0);
                            image[pixelIndex] = color.x; 
                            image[pixelIndex + 1] = color.y; 
                            image[pixelIndex + 2] = color.z; 
                        }
                        else
                        {
                            
                            image[pixelIndex] = scene.background_color.x; 
                            image[pixelIndex + 1] = scene.background_color.y; 
                            image[pixelIndex + 2] = scene.background_color.z; 
                        }
                
            }
        }
        
        
        write_ppm(filename.c_str(), image, width, height);

        delete[] image;
        std::cout << "done\n";

    }

    

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;    
    std::cout << "Total execution time: " << elapsed.count() << " seconds\n";
    std::cout << "\a" << std::flush;

    return 0;

}

















