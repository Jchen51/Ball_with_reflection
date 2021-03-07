//
//  main.cpp
//  Assignment 3
//
//  Created by Jason Chen on 11/26/15.
//  Copyright © 2015 Jason Chen. All rights reserved.
//
//  Part1 (Ray Casting), Part 2 (Reflection), Part 3 (Refraction)
//  When the image is rendered, Red Ball = reflection, White ball = refraction
//  Was unable to locate reason behind black grainy effect


#include <iostream>
#include <math.h>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define screenWidth 720
#define screenHeight 720

#define MAXOBJ 20

float pixelbuffert[720][720][3];//720x720x3(3 colors RGB)

int frontobj(std::vector<float> intersec);

//represents both coordinate points and vectors
class coordinate{
public:
    coordinate(){}
    coordinate(float x, float y, float z){
        xc = x;
        yc = y;
        zc = z;
    }
    float xc;
    float yc;
    float zc;
    
    float dotproduct(coordinate two){
        return ((this->xc * two.xc) + (this->yc * two.yc) + (this->zc * two .zc));
    }
    
    coordinate addvect(coordinate add){
        return coordinate((xc + add.xc),(yc + add.yc),(zc + add.zc));
    }
    
    coordinate scalevect(float f){
        return coordinate((xc*f), (yc*f), (zc*f));
        
        
    }
    
    coordinate neg(){
        return coordinate((-xc), (-yc), (-zc));
    }
    
    coordinate & operator=(const coordinate&);
    
    coordinate normalize(){
        float mag = sqrt((xc*xc) + (yc*yc)+(zc*zc));
        return coordinate ((xc/mag), (yc/mag), (zc/mag));
    }
    
    coordinate crossproduct(coordinate other){
        float totalx = (yc*other.zc)-(zc*other.yc);
        float totaly = (zc*other.xc)-(xc*other.zc);
        float totalz = (xc*other.yc)-(yc*other.xc);
        return coordinate(totalx, totaly, totalz);
    }
    
    
};

coordinate& coordinate::operator=(const coordinate& other){
    if ( this == &other)
        return *this;
    xc = other.xc;
    yc = other.yc;
    zc = other.zc;
    return *this;
}

class color{
public:
    
    color(){
        red = 0.5;
        green = 0.5;
        blue = 0.5;
        alpha = 0;
    }
    
    color(float nred, float ngreen, float nblue, float special){
        red = nred;
        green = ngreen;
        blue = nblue;
        alpha = special;
    }
    
    float red;
    float green;
    float blue;
    float alpha;
    
    color addcolor(color n){
        return color((red+n.red), (green+n.green), (blue+n.blue), alpha);
    }
    
    color scalecolor(float n){
        return color((red*n), (green*n), (blue*n), alpha);
    }
    
    color multicolor(color n){
        return color((red*n.red), (green*n.green), (blue*n.blue), alpha);
    }
    
    color avg(color n){
        return color(((red+n.red)/2), ((green+n.green)/2), ((blue+n.blue)/2), alpha);
    }
    
    //color must be less than or equal to 1 or greater than or equal to 0
    color clip(){
        float all = red + green + blue;
        float excess = all - 3;
        if (excess > 0){
            red = red + excess * (red/all);
            blue = blue + excess * (blue/all);
            green = green + excess * (green/all);
        }
        if (red > 1){red = 1;}
        if (green > 1){green = 1;}
        if (blue > 1){blue = 1;}
        if (red < 0){red = 0;}
        if (green < 0){green = 0;}
        if (blue < 0){blue = 0;}
        
        return color(red, green, blue, alpha);
    }
    
    bool eq(color c){
        if ((red == c.red) && (green == c.green) && (blue == c.blue) && (alpha == c.alpha))
            return true;
        else
            return false;
    }
    
    
};


class ray{
public:
    
    float indexofrefraction = 1;
    
    ray(){
        xc = NULL;
        yc = NULL;
        zc = NULL;
        xd = NULL;
        yd = NULL;
        zd = NULL;
    }
    ray(float startx, float starty, float startz, float nxd, float nyd, float nzd){
        xc = startx;
        yc = starty;
        zc = startz;
        xd = nxd;
        yd = nyd;
        zd = nzd;
    }
    
    ray(coordinate cent, coordinate dir){
        xc = cent.xc;
        yc = cent.yc;
        zc = cent.zc;
        xd = dir.xc;
        yd = dir.yc;
        zd = dir.zc;
    }
    
    //origin of ray
    float xc;
    float yc;
    float zc;
    
    //direction of ray
    float xd;
    float yd;
    float zd;
    
    
};

class newcam: public ray{
public:
    coordinate camright;
    coordinate camdown;
    
    newcam(){
        xc = 0;
        yc = 0;
        zc = 0;
        xd = 0;
        yd = 0;
        zd = 1;
        camright = coordinate(0,0,0);
        camdown = coordinate(0,0,0);
    }
    
    newcam(coordinate center, coordinate dir, coordinate ncamright, coordinate ncamdown){
        xc = center.xc;
        yc = center.yc;
        zc = center.zc;
        xd = dir.xc;
        yd = dir.yc;
        zd = dir.zc;
        camright = ncamright;
        camdown = ncamdown;
    }
};


class object{
public:
    object(){}
    virtual color getobjcolor(){
        return color(0, 0, 0, 0);
    }
    
    //calculate intersections between the object and R1. intersection points held in interray, roots held in root
    virtual bool intersection(ray R1, ray *interray, float &root){
        return true;
    }
    
    //get normals of different objects
    virtual coordinate getNormal(coordinate point){
        return coordinate(0,0,0);
    }
    
    float refr = 1;
};


class light{
public:
    
    //lights have centers, lights have colors
    coordinate center;
    color lcolor;
    light(){
        center = coordinate(0, 0, 0);
        lcolor = color(1, 1, 1, 0);
    }
    light(coordinate ncenter, color color1){
        center = ncenter;
        lcolor = color1;
    }
};

//function to get color at an intersection
color getColor(coordinate collisionstart, coordinate collisiondir, std::vector<object*> allobjects, int frontindex, std::vector<light*> lights, float a, float ambientlight){
    
    color frontcolor = allobjects.at(frontindex)->getobjcolor();
    coordinate frontnormal = allobjects.at(frontindex)->getNormal(collisionstart);
    
    if (frontcolor.alpha == 2){
        //checkerboard
        int square = (int) floor (collisionstart.xc) + (int)floor(collisionstart.yc) + (int)floor(collisionstart.zc);
        
        if ((square%2) == 0){
            frontcolor.red = 0;
            frontcolor.green = 0;
            frontcolor.blue = 0;
        }else{
            frontcolor.red = 1;
            frontcolor.green = 1;
            frontcolor.blue = 1;
        }
    }
    
    
    //final color to be returned
    color retcolor = frontcolor.scalecolor(ambientlight);
    
    
    //checks all lights
    for(int i = 0; i < lights.size(); i++){
        
        //directin of light to point
        coordinate lightdir = lights.at(i)->center.addvect(collisionstart.neg()).normalize();
        
        //cosine between intersection and light
        float cos = frontnormal.dotproduct(lightdir);
        
        if (cos > 0){//angle between close intersection and light is positive
            bool inshadow = false;
            
            //distance between collision point and position of first light
            coordinate dist = lights.at(i)->center.addvect(collisionstart.neg()).normalize();
            float dist_mag = sqrt((dist.xc)*(dist.xc)+(dist.yc)*(dist.yc)+(dist.zc)*(dist.zc));
            
            //starting from the collision, and heading towards the light, looking to see if intersects with anything
            ray secondray(collisionstart, lights.at(i)->center.addvect(collisionstart.neg()).normalize());
            
            std::vector<float> secondary;
            
            ray hold;
            float temproot;
            
            //holds every single intersection that exists between the second ray and the light
            for (int j = 0; j < allobjects.size() && !inshadow; j++){
                allobjects.at(j)->intersection(secondray, &hold, temproot);
                secondary.push_back(temproot);
                /*
                if (allobjects.at(j)->intersection(secondray, &hold, temproot)){
                    if (temproot > a){
                        if (temproot <= dist_mag){
                            inshadow = true;
                            break;
                        }
                    }
                }*/
            }
            
            //is it in the shadow of something?
            for (int k = 0; k < secondary.size(); k++){
                if (secondary.at(k) > a){
                    if (secondary.at(k) <= dist_mag){
                        inshadow = true;
                    }
                }
            }
            
            
            if (!inshadow){
                
                //add color of the object multiplied by the light intensity, becuase not in shadow
                retcolor = retcolor.addcolor(frontcolor.multicolor(lights.at(i)->lcolor).scalecolor(cos));
                
                if ((frontcolor.alpha > 0)&& (frontcolor.alpha <= 1)){ //object is shiny
                    
                    coordinate reflectdir = collisiondir.addvect(frontnormal.scalevect(2*(collisiondir.dotproduct(frontnormal))).neg()).normalize();
                    
                    /*
                    //math to determine direction of reflection
                    float dot1 = frontnormal.dotproduct(collisiondir.neg());
                    coordinate vect1 = frontnormal.scalevect(dot1);
                    coordinate add1 = vect1.addvect(collisiondir);
                    coordinate vect2 = add1.scalevect(2);
                    coordinate add2 = collisiondir.neg().addvect(vect2);
                    coordinate reflectdir = add2.normalize();
                    */
                    
                    float specular = reflectdir.dotproduct(lightdir);
                    if (specular > 0){
                        specular = pow(specular, 10);
                        retcolor = retcolor.addcolor(lights.at(i)->lcolor.scalecolor(specular*frontcolor.alpha));
                    }
                }
                
            }
        }
    }
    
    
    return retcolor.clip();
}

class sphere: public object{
public:
    coordinate center;
    float r; //radius
    color scolor;
    
    sphere();
    sphere(float newx, float newy, float newz, float newr, float newred, float newgreen, float newblue, float newspec, float refractive){
        
        center = coordinate(newx, newy, newz);
        r = newr;
        scolor = color(newred, newgreen, newblue, newspec);
        refr = refractive;
    }
    
    sphere(coordinate newcenter, float radius, color spherecol, float refractive){
        center = newcenter;
        r = radius;
        scolor = spherecol;
        refr = refractive;
    }
    
    //get the normal of the sphere at a point
    virtual coordinate getNormal(coordinate point){
        coordinate normalvec = point.addvect(center.neg()).normalize();
        return normalvec;
    }
    
    //must be an inherited function so can be called by object class
    color getobjcolor(){
        return scolor;
    }
    
    virtual bool intersection(ray R1, ray *interray, float &root);
    
};

//calculate intersection between a ray and the sphere
bool sphere::intersection(ray R1, ray *interray, float &root){
    float b = 2 * (R1.xd*(R1.xc - center.xc)+R1.yd*(R1.yc - center.yc)+R1.zd*(R1.zc - center.zc)); //sphere is second?
    float c = (R1.xc-center.xc)*(R1.xc-center.xc) + (R1.yc-center.yc)*(R1.yc-center.yc) + (R1.zc-center.zc)*(R1.zc-center.zc) - r*r;
    
    float disc = b*b-4*c;
    if (disc < 0){
        //no intersections
        root = -1;
        return false;
    }else{
        float t = (-b-sqrt(disc))/2;
        ray newray;
        if (t >= 0){
            //first root is closest or only one root
            newray.xc = R1.xc + t * R1.xd;
            newray.xd = R1.xd;
            newray.yc = R1.yc + t * R1.yd;
            newray.yd = R1.yd;
            newray.zc = R1.zc + t * R1.zd;
            newray.zd = R1.zd;
            interray = &newray;
            root = t;
        }
        else if ((t + sqrt(disc)) >= 0){
            //second root is closest
            t += disc;
            newray.xc = R1.xc + t * R1.xd;
            newray.xd = R1.xd;
            newray.yc = R1.yc + t * R1.yd;
            newray.yd = R1.yd;
            newray.zc = R1.zc + t * R1.zd;
            newray.zd = R1.zd;
            interray = &newray;
            root = t;
        }
        else{
            //both roots are futher than camera
            root = -1;
            return false;
        }
        return true;
    }
}

//plane class
class plane: public object{
public:
    coordinate normal;
    float distance;//to center of scene
    color color1;
    
    plane(){
        normal = coordinate(1, 0, 0);
        
        distance = 0.0;
        color1 = color(0.5, 0.5, 0.5, 0);
        
        refr = 1;
    }
    
    plane(coordinate newnormal, float newdist, color newcol){
        normal = newnormal;
        distance = newdist;
        color1 = newcol;
        
        refr = 1;
    }
    
    virtual bool intersection(ray R1, ray *interray, float &root);
    
    
    color getobjcolor(){
        return color1;
    }
    
    //must be inherited class due to needing to call on object class
    virtual coordinate getNormal(coordinate point){
        return normal;
    }
    
    
};

//calculate intersection between plane and ray R1, return roots to &root
bool plane::intersection(ray R1, ray *interray, float &root){
    coordinate raystart(R1.xc, R1.yc, R1.zc);
    coordinate raydir(R1.xd, R1.yd, R1.zd);
    ray newray;
    
    float d = raydir.dotproduct(normal);
    if (d == 0){
        //ray is parallel to plane, no intersect
        root = -1;
        return false;
    }
    else{
        float b = normal.dotproduct(raystart.addvect(normal.scalevect(distance).neg()));
        
        float t = -1*(b/d);
        
        newray.xc = R1.xc + t*R1.xd;
        newray.xd = R1.xd;
        newray.yc = R1.yc + t*R1.yd;
        newray.yd = R1.yd;
        newray.zc = R1.zc + t*R1.zd;
        newray.zd = R1.zd;
        
        interray = &newray;
        
        root = t;
        
        return true;
    }
}

//unfortunately, not enough time
class ellipse: public object{
    bool intersection(ray R1, coordinate*interray);
};

//calculate the frontmost object, closest to the camera
int frontobj(std::vector<float> intersec){
    int min;
    
    if (intersec.size() == 0){
        //no intersections
        return -1;
    }else if(intersec.size() == 1){
        //one intersec
        if (intersec.at(0) > 0){
            //first is in front
            return 0;
        }
        else {
            //only intersec is neg
            return -1;
        }
    }else{
        //more than one intersection
        //find maxvalue of vector
        
        float max = 0;
        for (int i = 0; i < intersec.size(); i++){
            if (max < intersec.at(i)){
                max = intersec.at(i);
            }
        }
        
        //now need to find minimum positive value
        if (max > 0){
            for (int i = 0; i < intersec.size(); i++){
                if (intersec.at(i) > 0 && intersec.at(i) <= max){
                    max = intersec.at(i);
                    min = i;
                }
            }
            return min;
        }else{
            //all intersec neg
            return -1;
        }
        
    }
    
}


//our main raytrace function, called recursively to calculate refraction and reflection
void raytrace(ray R1, int recurdepth, color &col, std::vector<object*> allobj, std::vector<light*> lights, float a, float ambientlight){
    
    //stop calling self
    if(recurdepth == 0)
        return;
    
    float temproot;
    
    std::vector<float> intersections;
    ray intersecpoint;
    std::vector<coordinate> normals;
    
    //origin of incoming ray
    coordinate start(R1.xc, R1.yc, R1.zc);
    
    //direction of incoming ray
    coordinate dir(R1.xd, R1.yd, R1.zd);
    
    //get the roots
    for (int i = 0; i < allobj.size(); i++){
        allobj.at(i)->intersection(R1, &intersecpoint, temproot);
        intersections.push_back(temproot);
    }
    
    //frontmost object
    int newfront = frontobj(intersections);
    
    if (newfront == -1) //front object either behind camera or doesn't exist
        return;
    else {
        
    //point of collision
    coordinate colstart = start.addvect(dir.scalevect(intersections.at(newfront)));

    for (int i = 0; i < allobj.size();i++)
        normals.push_back(allobj.at(i)->getNormal(colstart));
    
    
    color final;
    color defaultcol;
    
    
    //get color of intersection
        final = getColor(colstart, dir, allobj, newfront, lights, a, ambientlight);
    
    
    //calculate reflection
    
    if (final.alpha > 0 && final.alpha <= 1){
    

        coordinate refdir;
        
        //normal of reflection
        coordinate refnormal = normals.at(newfront).normalize();
        
        //new direction of reflected ray
        refdir = dir.addvect(refnormal.scalevect(2*(dir.dotproduct(refnormal))).neg()).normalize();
        
        //ray being reflected out
        ray reflectout(colstart, refdir);
        
        color reflectcol;
        
        raytrace(reflectout, recurdepth-1, reflectcol, allobj, lights, a, ambientlight);
        
        final = final.addcolor(reflectcol.scalecolor(final.alpha));
       
        final.clip();
    }
    
    //calculate refraction n1sin1 = n2sin2
    if (allobj.at(newfront)->refr != 1.000){
        
        coordinate refrdir;
        
        //incoming vector = dir
        //collision = colstart
        coordinate refrnormal = normals.at(newfront);
        
        //ray coming in
        
        float n1 = R1.indexofrefraction;
        float n2 = allobj.at(newfront)->refr;
        
        float n = n1/n2;
        
        float newindex = allobj.at(newfront)->refr;
        
        float cosI = -1*refrnormal.dotproduct(dir);
        float sinI2 = n*n*(1.0-(cosI*cosI));
        if (sinI2 > 1.0)
            return;
        
        //coordinate refractdirect = dir.addvect(newnormal);
        
        color refractcol;
        
        coordinate refractdirect = dir.scalevect(n).addvect(refrnormal.scalevect(n+sqrt(1.0-sinI2)).neg());
        
        //outgoing refracting ray
        ray refractout(colstart, refractdirect);
        
        refractout.indexofrefraction = newindex;
        
        raytrace(refractout, recurdepth-1, refractcol, allobj, lights, a, ambientlight);
        
        final = final.addcolor(refractcol);
        
        final.clip();
    }
    
    col = final;
    }
}

//the base logic behind creating the shapes, camera, and lights
void logic(void){
    float xamount;
    float yamount;
    
    coordinate O(0,0,0); //origin
    coordinate X(1,0,0); //default X
    coordinate Y(0,1,0); //default Y
    coordinate Z(0,0,1); //default Z
    
    coordinate center(3, 1.5, -4);
    
    coordinate lookat(0,0,0);
    coordinate diff((center.xc-lookat.xc),(center.yc - lookat.yc), (center.zc - lookat.zc));
    
    //direction that the camera is looking at
    coordinate camdir = diff.neg().normalize();
    
    coordinate camright = Y.crossproduct(camdir).normalize();
    coordinate camdown = camright.crossproduct(camdir);
    
    newcam cam(center, camdir, camright, camdown);
    
    coordinate dirvec;
    
    //property of ambient light
    float ambientlight = 0.2;
    
    coordinate sphere2cent(1.75, -0.5, 0);
    coordinate sphere3cent(1.75, -0.5, 1.75);
    
    //some test colors
    color white(1, 1, 1, 0);
    color black(0, 0, 0, 0);
    color blue(0,0,1,0);
    color red(1, 0, 0, 0);
    color grey(0.5, 0.5, 0.5, 0);
    color maroon(0.5, 0.25, 0.25, 0);
    color floorcolor(1, 1, 1, 2);
    
    //objects
    sphere sphere1(0 , 0 , 0 , 1, 0.5, 0.25, 0.25, 0.3 /*shininess*/, 1);
    sphere sphere2(sphere2cent, 0.5, grey, 1.3);
    sphere sphere3(sphere3cent, 0.5, blue, 1);
    plane plane1(Y, -1 /* underneath the sphere */, floorcolor);
    
    //light
    coordinate lightloc(-7, 10, -10);
    light light1(lightloc, white);
    
    std::vector<object*> allobjects;
    allobjects.push_back(dynamic_cast<object*>(&sphere1));
    allobjects.push_back(dynamic_cast<object*>(&plane1));
    allobjects.push_back(dynamic_cast<object*>(&sphere2));
    allobjects.push_back(dynamic_cast<object*>(&sphere3));
    
    
    std::vector<light*> lights;
    lights.push_back(dynamic_cast<light*>(&light1));
    
    ray drawray;
    
    for (int i = 0; i < screenWidth; i++){
        for (int j = 0; j < screenHeight; j++){
            
            
         
            // the image is square
            xamount = (i + 0.5)/screenWidth;
            yamount = ((screenHeight - j) + 0.5)/screenHeight;
            
            coordinate camorigin = center;
            coordinate camraydir = camdir.addvect(camright.scalevect(xamount-0.5).addvect(camdown.scalevect(yamount-0.5))).normalize();
            
            
            //ray containing the origin and direction of our camera
            ray newray(camorigin, camraydir);
            
            std::vector<float> intersections; //array of intersections
            ray intersecpoint;
            std::vector<ray> intersecpoints; //array of intersection points
            
            float temproot;
            
            for (int k = 0; k < allobjects.size(); k++){
                
                
                allobjects.at(k)->intersection(newray, &intersecpoint, temproot);
                intersecpoints.push_back(intersecpoint);
                intersections.push_back(temproot);
            }
            
            
            //object closest to the camera, need to return color
            int frontindex = frontobj(intersections);
            
            if (frontindex == -1){//missed, background = black
                
                
                pixelbuffert[j][i][0] = 0.0;
                pixelbuffert[j][i][1] = 0.0;
                pixelbuffert[j][i][2] = 0.0;
                
            }else{//hit, found
                
                if (intersections.at(frontindex) > 0.00000000001){
                    
                    //collision is at the origin + t * direction, which is the front intersection which we are looking at
                    
                    //depth of tracing
                    int tracer = 4;
                    
                    color final;
                    
                    raytrace(newray,  tracer, final, allobjects, lights, 0.0000000001, ambientlight);
                    
                    //color each pixel one by one
                    pixelbuffert[j][i][0] = final.red;
                    pixelbuffert[j][i][1] = final.green;
                    pixelbuffert[j][i][2] = final.blue;
                    
                }
            }
            
        }
    }

}

void init(void){
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    logic();
}

//draw whole thing using glDrawPixels
void display(){
    
    //endif
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(screenHeight, screenWidth, GL_RGB, GL_FLOAT, pixelbuffert);
    glFlush();
    
}

int main(int argc, char * argv[]) {
    glutInit(&argc, argv);
    
    //Set Display Mode
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
    
    //Set the window size
    glutInitWindowSize(screenWidth, screenHeight);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("HW3");
    
    init();
    
    //Call "display" function
    glutDisplayFunc(display);
    
    glutMainLoop();
}
