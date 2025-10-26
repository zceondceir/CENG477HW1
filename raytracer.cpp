#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <chrono>

typedef unsigned char RGB[3];



//render.h



inline float determinant (const parser::Vec3f &a, const parser::Vec3f &b, const parser::Vec3f &c) {
    return a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z) + a.z * (b.x * c.y - b.y * c.x);
}


//const parser::Vec3f &o, const parser::Vec3f &d

inline void intersectsTriangle(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Triangle &t, parser::HitRecord& h, const parser::Scene& s){


    const parser::Vec3f &a = s.vertex_data[t.indices.v0_id-1];
    const parser::Vec3f &b = s.vertex_data[t.indices.v1_id-1];
    const parser::Vec3f &c = s.vertex_data[t.indices.v2_id-1];

    parser::Vec3f ab = b - a;
    parser::Vec3f ac = c - a;
    parser::Vec3f ao = o - a;

    parser::Vec3f normal = ab ^ ac;

    if(normal * d > 0){
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

inline void intersectsSphere(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Sphere &sp, parser::HitRecord& h, const parser::Scene& s){

    const parser::Vec3f &c = s.vertex_data[sp.center_vertex_id - 1];
    
    parser::Vec3f oc = o - c;

    float B = d * (oc) * 2;
    float C = (oc) * (oc) - (sp.radius * sp.radius);

    float disc = B * B - 4 * C;

    if(disc < 0.0f)    return;
    

    float sqrt = std::sqrt(disc);
    float t1 =  (-B + sqrt) / (2.0f);
    float t2 =  (-B - sqrt) / (2.0f);
    float t = std::min(t1,t2);
    if(t < 0.0f ) return;
 
    if( h.t > t){
        h.hit = true;
        h.t = t;  
        h.point = o + d * t;
        h.normal = (h.point - c).normalized();
        if(h.normal * d > 0) h.normal = h.normal * -1;
        h.Mid = sp.material_id;
    }
    

}

inline void intersectsPlane(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Plane &p, parser::HitRecord& h, const parser::Scene& s){

    const parser::Vec3f &n = p.normal;
    const parser::Vec3f &a = s.vertex_data[p.center_vertex_id - 1];
    
    float dn = d * n;

    if (fabs(dn) < 1e-8f) return;

    float t = (a - o) * n / dn ;

    if (t < s.shadow_ray_epsilon) return;

    if(h.t > t){
        h.hit = true;
        h.t = t;
        h.normal = n;
        h.point = o + d * t;
        h.Mid = p.material_id;
    }

}

void intersectsCylinder(const parser::Ray &r, const parser::Cylinder &cyl, parser::HitRecord& h, const parser::Scene& s)
{
    parser::Vec3f o = r.origin;
    parser::Vec3f d = r.direction;
    parser::Vec3f C = s.vertex_data[cyl.center_vertex_id - 1];
    parser::Vec3f v = cyl.axis.normalized();  // eksen yönü normalize olmalı

    // --- 1. Sonsuz silindir kesişimi ---
    parser::Vec3f delta = o - C;
    float dv = d*v;
    float delta_v = delta * v;

    parser::Vec3f d_perp = d - v * dv;
    parser::Vec3f delta_perp = delta - (v * delta_v);

    float A = (d_perp * d_perp);
    float B = 2.0f * (d_perp * delta_perp);
    float Cc = (delta_perp * delta_perp) - cyl.radius * cyl.radius;

    float disc = B * B - 4.0f * A * Cc;
    if (disc < 0.0f) return; // kesişim yok

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

        // --- 2. Yükseklik kontrolü ---
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

    // --- 3. Üst ve alt kapak kontrolü ---
    float halfH = cyl.height * 0.5f;
    parser::Vec3f capTopCenter = C + v * halfH;
    parser::Vec3f capBotCenter = C - v * halfH;

    // Üst kapak
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

    // Alt kapak
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



inline void intersectsMesh(const parser::Ray &r, const parser::Mesh &m, parser::HitRecord& h, const parser::Scene& s, bool BFC){

    for(const parser::Face &f : m.faces){

        const parser::Vec3f a = s.vertex_data[f.v0_id-1];
        parser::Vec3f b = s.vertex_data[f.v1_id-1];
        parser::Vec3f c = s.vertex_data[f.v2_id-1];

        parser::Vec3f o = r.origin;
        parser::Vec3f d = r.direction;

        float detA = determinant(a-b,a-c,d);

        if (fabs(detA) < 1e-8) continue;

        float beta = determinant(a-o,a-c,d)/detA;
        float gamma = determinant(a-b,a-o,d)/detA;
        float tRay = determinant(a-b,a-c,a-o)/detA;

        parser::Vec3f normal = ((b - a) ^ (c - a)).normalized();

        if(BFC){

            
            if(normal  * d > 0 ){
                continue;
            }
        }
        
    
        if( beta >= 0 && gamma >= 0 && beta + gamma <= 1 && tRay < h.t && tRay > 1e-6f ){
            h.hit = true;
            h.t = tRay;
            h.normal = normal;
            h.point = r.at(tRay);
            h.Mid = m.material_id;
        } 
    }
    
    //          v0 = a
    //
    //     v1=b         v2 = c

}


inline void intersectsMesh(const parser::Vec3f &o, const parser::Vec3f &d, const parser::Mesh &m, parser::HitRecord& h, const parser::Scene& s){

    for(const parser::Face &f : m.faces){

        const parser::Vec3f &a = s.vertex_data[f.v0_id-1];
        const parser::Vec3f &b = s.vertex_data[f.v1_id-1];
        const parser::Vec3f &c = s.vertex_data[f.v2_id-1];

        parser::Vec3f ab = b - a;
        parser::Vec3f ac = c - a;
        parser::Vec3f ao = o - a;

        parser::Vec3f normal = ab ^ ac;

        if(normal * d > 0){
            continue;;
        }

        float detA = determinant(ab,ac,d);

        if (fabs(detA) < 1e-8) continue;;

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
    //          v0 = a
    //
    //     v1=b         v2 = c

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

    return s - e;


}


void renderRay(const parser::Scene &scene, const parser::Ray &r, parser::HitRecord & closestHit){
        
    for(const parser::Triangle &t : scene.triangles) {

        intersectsTriangle(r.origin,r.direction,t,closestHit,scene);

    }
    for(const parser::Sphere &s : scene.spheres) {

        intersectsSphere(r.origin,r.direction,s,closestHit,scene);
            
    }
    for(const parser::Plane &p : scene.planes) {
        
        intersectsPlane(r.origin,r.direction,p,closestHit,scene);
            
    }    

    for(const parser::Cylinder &c : scene.cylinders){

        intersectsCylinder(r,c,closestHit,scene);

    }

    for(const parser::Mesh &m : scene.meshes){

        intersectsMesh(r,m,closestHit,scene,true);

    }

}

bool inShadow(const parser::Scene &scene,const parser::Ray &r, const float t){

    parser::HitRecord hr;
    hr.t = t;

    for(const parser::Triangle &t : scene.triangles) {

        intersectsTriangle(r.origin,r.direction,t,hr,scene);
        if(hr.hit) return true;

    }
    
    for(const parser::Sphere &s : scene.spheres) {

        intersectsSphere(r.origin,r.direction,s,hr,scene);
        if(hr.hit) return true;
            
    }
    for(const parser::Plane &p : scene.planes) {
        
        intersectsPlane(r.origin,r.direction,p,hr,scene);
        if(hr.hit) return true;
            
    }    

    for(const parser::Cylinder &c : scene.cylinders){

        intersectsCylinder(r,c,hr,scene);
        if(hr.hit) return true;

    }

    for(const parser::Mesh &m : scene.meshes){

        intersectsMesh(r,m,hr,scene,true);
        if(hr.hit) return true;

    }
    
    return false;
    
}

parser::Ray Reflect(const parser::HitRecord &hitrecord, const parser::Scene &s, const parser::Ray &r ) {
    parser::Ray incomingRay = r;
    parser::Ray refRay;

    
    parser::Vec3f wo = incomingRay.direction.normalized() * -1;
    parser::Vec3f normal = hitrecord.normal.normalized();

   
    parser::Vec3f wr = (wo * -1) + normal * (normal * wo) * 2.0f;

    refRay.direction = wr.normalized();
    refRay.origin = hitrecord.point + normal * 1e-4f;
    
    refRay.depth = incomingRay.depth + 1;

    return refRay;
}

parser::Vec3i GetColor(const parser::HitRecord &hitrecord, const parser::Scene &s, const parser::Ray &r ){


    if ( r.depth > s.max_recursion_depth ){
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

        parser::Vec3f antiCollision = hitrecord.point + hitrecord.normal * (1e-8);
        parser::Vec3f shadowVector = (light.position - antiCollision).normalized();
        float t = (light.position - hitrecord.point).len();
        parser::Ray shadowRay(antiCollision,shadowVector,0);
        

        if(inShadow(s,shadowRay,t)){
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

            if(reflection.depth < s.max_recursion_depth){
                
                renderRay(s,reflection,reflectionRecord);

                if (reflectionRecord.hit)
                {
                    parser::Vec3i refcolor = GetColor(reflectionRecord,s,r);
                    color_r += material.mirror.x * refcolor.x;
                    color_g += material.mirror.y * refcolor.y;
                    color_b += material.mirror.z * refcolor.z; 
                }   

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
    // Sample usage for reading an XML scene file
    
    auto start = std::chrono::high_resolution_clock::now();

    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    int cameraId = 1;

    //rendering part

    



    for(const parser::Camera& c : scene.cameras){

        

        std::string filename = "ouroutputs/" + std::string(argv[1]) + "_camera" + std::to_string(cameraId++) + ".ppm";

        int width = c.image_width;
        int height = c.image_height;

        unsigned char* image = new unsigned char [width * height * 3];

        int i = 0;
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                // (sx , sy)

                parser::Ray r = parser::Ray(c.position, getPixelRay(c, x, y).normalized() , 0);

                parser::HitRecord closestHit;

                renderRay(scene,r,closestHit);
                
                int pixelIndex = (y * width + x) * 3; 

                if(closestHit.hit)
                        {                    
                            parser::Vec3i color = GetColor(closestHit, scene, r);
                            image[pixelIndex] = color.x; // R
                            image[pixelIndex + 1] = color.y; // G
                            image[pixelIndex + 2] = color.z; // B
                        }
                        else
                        {
                            
                            image[pixelIndex] = scene.background_color.x; // R
                            image[pixelIndex + 1] = scene.background_color.y; // G
                            image[pixelIndex + 2] = scene.background_color.z; // B
                        }
                
            }
        }
        
        
        write_ppm(filename.c_str(), image, width, height);

    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;    // saniye cinsinden fark
    std::cout << "Total execution time: " << elapsed.count() << " seconds\n";

}

















