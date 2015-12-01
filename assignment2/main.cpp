#ifdef __APPLE__
#include <GLUT/glut.h>

#include <OpenGL/gl.h>
#else
#include <GL/glut.h>
#endif

#include <cmath>
#include <vector>
#include "neheTGATexture.h"

#include"Vector.h"
#include "Camera.h"
#include "SpringMesh.h"
#include "gauss.h"
#include "Particle.h"



#include <Magick++.h>
#include "FrameGrabber.h"

using namespace Magick;

static string MYPATH = "/Users/yuyaolong/Desktop/animation";
static string MYFILENAME = "flagFlock";

#define HDWIDTH	      1920		// HD image dimensions
#define HDHEIGHT      1080
#define WIDTH	      (HDWIDTH/2)	// window dimensions = 1/2 HD
#define HEIGHT	      (HDHEIGHT/2)
#define STARTFRAME    0

FrameGrabber framegrabber(MYPATH, MYFILENAME, HDWIDTH, HDHEIGHT, STARTFRAME);

int Recording = false;


#define automaticSpeed 20

//#define OBSTACLENO 1

#define ERRORTHRESH 0.001

#define Cd 0.5
#define Cl 0.5

#define Ksita 0.4
#define Dsita 0.1

#define SITA0 PI
#define EPSILON  0.000001

double hStep=0.05;

int SCREENWIDTH = 1000;
int SCREENHEIGHT = 800;

Camera *camera;

bool  resetSign = false;
bool showGrid = false;
bool clothShow = true;
bool explode = false;
SpringMesh springMesh;


unsigned int windCount = 0;
Vector3d gravity(0,-1,0);

//particle thing
int ParticleNum = 150;
float particleSize = 0.2;
std::vector<Particle>particles;
Vector3d circleCenterPosition(0,6,0);

std::vector<Particle>fireWorks;


imageTGA texture;

enum calcuState
{
    SPRINGCAL=0,
    PARTICLECAL,
    FIREWORKCAL
};


Vector4d color[10]={ Vector4d(1.00f, 0.00f, 0.00f, 1.0f),
                    Vector4d(1.00f, 0.50f, 0.00f, 1.0f),
                    Vector4d(0.50f, 1.00f, 0.00f, 1.0f),
                    Vector4d(0.20f, 0.20f, 1.00f, 1.0f),
                    Vector4d(0.10f, 1.00f, 0.20f, 1.0f),
                    Vector4d(0.80f, 0.90f, 0.20f, 1.0f),
                    Vector4d(0.90f, 0.10f, 0.40f, 1.0f),
                    Vector4d(0.10f, 0.10f, 0.90f, 1.0f),
                    Vector4d(0.10f, 0.90f, 0.90f, 1.0f),
                    Vector4d(1.00f, 0.20f, 0.10f, 1.0f)
};


// draws a simple grid
void makeGrid()
{
    //Setting  material
    GLfloat ball_mat_ambient[]  = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat ball_mat_diffuse[]  = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat ball_mat_specular[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat ball_mat_emission[] = {0.0f, 0.0f, 0.0f, 1.0f};
    GLfloat ball_mat_shininess  = 10.0f;
    
    glMaterialfv(GL_FRONT, GL_AMBIENT,   ball_mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   ball_mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  ball_mat_specular);
    glMaterialfv(GL_FRONT, GL_EMISSION,  ball_mat_emission);
    glMaterialf (GL_FRONT, GL_SHININESS, ball_mat_shininess);
    
    glLineWidth(1.0);
    
    for (float i=-12; i<12; i++) {
        for (float j=-12; j<12; j++) {
            glBegin(GL_LINES);
            glVertex3f(i, 0, j);
            glVertex3f(i, 0, j+1);
            glEnd();
            glBegin(GL_LINES);
            glVertex3f(i, 0, j);
            glVertex3f(i+1, 0, j);
            glEnd();
            
            if (j == 11){
                glBegin(GL_LINES);
                glVertex3f(i, 0, j+1);
                glVertex3f(i+1, 0, j+1);
                glEnd();
            }
            if (i == 11){
                glBegin(GL_LINES);
                glVertex3f(i+1, 0, j);
                glVertex3f(i+1, 0, j+1);
                glEnd();
            }
        }
    }
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glVertex3f(-12, 0, 0);
    glVertex3f(12, 0, 0);
    glEnd();
    glBegin(GL_LINES);
    glVertex3f(0, 0, -12);
    glVertex3f(0, 0, 12);
    glEnd();
    glLineWidth(1.0);
}

void mouseEventHandler(int button, int state, int x, int y) {
    // let the camera handle some specific mouse events (similar to maya)
    camera->HandleMouseEvent(button, state, x, y);
    glutPostRedisplay();
}

void motionEventHandler(int x, int y) {
    // let the camera handle some mouse motions if the camera is to be moved
    camera->HandleMouseMotion(x, y);
    glutPostRedisplay();
}

double calAngle(Vector3d& a, Vector3d& b)
{
    return acos(a.normalize() * b.normalize())*180/PI;
}


void particlesGenerator(int number)
{
    for (int i=0; i<number; i++) {
        Particle tmp(Vector3d(-10,gauss(7, 1, 1),gauss(0, 2, 1)), Vector3d(gauss(2.5, 0.5, 1),2,0), Vector3d(0,0,0), Vector4d(0,0.5,0.2,1), 1, 0, particleSize, false, "chaos");//position,velocity,color,mass,lifespan,pointsize,stopSign
        particles.push_back( tmp );
        
        
    }
    
}

void fireWorkGenerator(Vector3d center, unsigned int number)
{
    unsigned int index = rand()%10;
    for (int i = 0; i<number; i++) {
        fireWorks.push_back(Particle(Vector3d(center.x, center.y, center.z), Vector3d(gauss(0, 5, 1),gauss(0, 5, 1),gauss(0, 5, 1)), Vector3d(0,0,0), color[index], 1, 0, 2, false,"firework"));
        
    }
}


void traceGenerator(Vector3d center, unsigned int number)
{
    for (int i = 0; i<number; i++) {
        fireWorks.push_back(Particle(Vector3d(center.x, center.y, center.z), Vector3d(gauss(0, 0.2, 1),gauss(0, 0.2, 1),gauss(0, 0.2, 1)), Vector3d(0,0,0), Vector4d(0.3,0.3,0.3,0), 1, 0, 2, false,"traces"));
        
    }
}




void myDisplay(void)
{
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // draw the camera created in perspective
    camera->PerspectiveDisplay(WIDTH, HEIGHT);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glPushMatrix();
    glTranslatef(0,-EDGELENGTH*LINES/2,0);
    if (showGrid)
        makeGrid();
    glPopMatrix();
     
    //spring mesh
    GLfloat ball_mat_ambient[]  = {0.0f, 0.0f, 0.0f, 1.0f};
    GLfloat ball_mat_diffuse[]  = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat ball_mat_specular[] = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat ball_mat_emission[] = {0.7f, 0.7f, 0.7f, 1.0f};
    GLfloat ball_mat_shininess  = 8.0f;
    
    glMaterialfv(GL_FRONT, GL_AMBIENT,   ball_mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   ball_mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  ball_mat_specular);
    glMaterialfv(GL_FRONT, GL_EMISSION,  ball_mat_emission);
    glMaterialf (GL_FRONT, GL_SHININESS, ball_mat_shininess);
    

    for (int i=0; i<FACENUMBER; i++) {
 
        int T1 = springMesh.faces[i].pointIndices[0];
        int T2 = springMesh.faces[i].pointIndices[1];
        int T3 = springMesh.faces[i].pointIndices[2];
        
        if (clothShow) {
            
            
            glBindTexture(GL_TEXTURE_2D, texture.getId());        //使用贴图纹理
            glBegin(GL_TRIANGLES);
            glTexCoord2f(springMesh.points[T1].textureCoordinate.x, springMesh.points[T1].textureCoordinate.y);
            glVertex3f(springMesh.points[T1].xposition.x, springMesh.points[T1].xposition.y, springMesh.points[T1].xposition.z);
            glTexCoord2f(springMesh.points[T2].textureCoordinate.x, springMesh.points[T2].textureCoordinate.y);
            glVertex3f(springMesh.points[T2].xposition.x, springMesh.points[T2].xposition.y, springMesh.points[T2].xposition.z);
            glTexCoord2f(springMesh.points[T3].textureCoordinate.x, springMesh.points[T3].textureCoordinate.y);
            glVertex3f(springMesh.points[T3].xposition.x, springMesh.points[T3].xposition.y, springMesh.points[T3].xposition.z);
            
            
            glEnd();
            
        }
        else
        {
            glLineWidth(1);
        
            glBegin(GL_LINES);
            glVertex3f(springMesh.points[T1].xposition.x, springMesh.points[T1].xposition.y, springMesh.points[T1].xposition.z);
            glVertex3f(springMesh.points[T2].xposition.x, springMesh.points[T2].xposition.y, springMesh.points[T2].xposition.z);
            glEnd();

            glBegin(GL_LINES);
            glVertex3f(springMesh.points[T2].xposition.x, springMesh.points[T2].xposition.y, springMesh.points[T2].xposition.z);
            glVertex3f(springMesh.points[T3].xposition.x, springMesh.points[T3].xposition.y, springMesh.points[T3].xposition.z);
            glEnd();
            
            glBegin(GL_LINES);
            glVertex3f(springMesh.points[T3].xposition.x, springMesh.points[T3].xposition.y, springMesh.points[T3].xposition.z);
            glVertex3f(springMesh.points[T1].xposition.x, springMesh.points[T1].xposition.y, springMesh.points[T1].xposition.z);
            glEnd();
        }
        
        
    }
    
    
    
    //particle
    for (int i = 0; i < particles.size(); ++i) {
        if (particles[i].getStopSign() == false) {
            GLfloat mat_ambient[]  = {0, 0.2, 0.5, 1.0f};
            GLfloat mat_diffuse[]  = {0, 0.2, 0.5, 1.0f};
            GLfloat mat_specular[] = {0, 0.2, 0.5, 1.0f};
            GLfloat mat_emission[] = {0, 0.2, 0.5, 1.0f};
            GLfloat mat_shininess  = 8.0f;
            glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
            glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
            glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
            glMaterialfv(GL_FRONT, GL_EMISSION,  mat_emission);
            glMaterialf (GL_FRONT, GL_SHININESS, mat_shininess);
            
        }
        else
        {
            float r=particles[i].getColor().x;
            float g=particles[i].getColor().y;
            float b=particles[i].getColor().z;
            float al=particles[i].getColor().w;
            
            GLfloat mat_ambient[]  = {r, g, b, al};
            GLfloat mat_diffuse[]  = {r, g, b, al};
            GLfloat mat_specular[] = {r, g, b, al};
            GLfloat mat_emission[] = {r, g, b, al};
            GLfloat mat_shininess  = 8.0f;
            glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
            glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
            glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
            glMaterialfv(GL_FRONT, GL_EMISSION,  mat_emission);
            glMaterialf (GL_FRONT, GL_SHININESS, mat_shininess);
        }
        

        glPushMatrix();
        glTranslatef(particles[i].getPosition().x, particles[i].getPosition().y, particles[i].getPosition().z);
        glutSolidSphere(particleSize, 10, 10);
        glPopMatrix();
    }
    
    
    //firework

    if (fireWorks.size() > 0) {
        for (int i = 0; i < fireWorks.size(); ++i) {
            float r=fireWorks[i].getColor().x;
            float g=fireWorks[i].getColor().y;
            float b=fireWorks[i].getColor().z;
            float al=fireWorks[i].getColor().w;
            
            GLfloat mat_ambient[]  = {r, g, b, al};
            GLfloat mat_diffuse[]  = {r, g, b, al};
            GLfloat mat_specular[] = {r, g, b, al};
            GLfloat mat_emission[] = {r, g, b, al};
            GLfloat mat_shininess  = 8.0f;
            glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
            glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
            glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
            glMaterialfv(GL_FRONT, GL_EMISSION,  mat_emission);
            glMaterialf (GL_FRONT, GL_SHININESS, mat_shininess);
            
            glPointSize(fireWorks[i].getPointSize());
            glBegin(GL_POINTS);
            glVertex3f(fireWorks[i].getPosition().x, fireWorks[i].getPosition().y, fireWorks[i].getPosition().z);
            glEnd();
        }
    }
    
}






//*************floacking acceleration addd ************************//

float mini(float a, float b)
{
    if (a>b) {
        return b;
    }
    else
    {
        return a;
    }
}

void flockAcceleration(std::vector<Vector3d>& states, std::vector<Vector3d>& statesA, unsigned int BOIDNUMBER)
{
    double KA = gauss(0.05, 0.01, 1);
    double KV = gauss(0.003, 0.001, 1);
    double KC = gauss(0.02, 0.007, 1);
    double R1 = 8;
    double R2 = 20;
    double AMAX = 0.1;
    
//    //Computer environmental acceleration on all particals
//    for (int i = 0; i<BOIDNUMBER; i++) {
//        statesA[i+BOIDNUMBER] = 0;
//    }
    
    //Add interactive acceleration to all particals
    for (int i = 0; i<BOIDNUMBER; i++) {
        if (particles[i].getStopSign() == true) {
            
            for (int j=0; j<BOIDNUMBER; j++) {
                if (i!=j) {
                    
                    Vector3d Xij = states[j] - states[i];
                    double distance = Xij.norm();
                    //avoidance acceleration
                    Vector3d Aija = (Xij.normalize()*(-1*KA/distance)) / (particles[i].getMass());
                    //velocity matching acceleration
                    Vector3d Aijv = (KV * (statesA[j] - statesA[i])) / (particles[i].getMass());
                    //centering acceleration
                    Vector3d Aijc = (KC * Xij) / (particles[i].getMass());
                    
                    Vector3d acceleration(0,0,0);
                    
                    float ar = AMAX;
                    if (0 != Aija.norm()) {
                        acceleration = mini(ar, Aija.norm())*Aija.normalize();
                    }
                    else
                    {
                        acceleration = Aija;
                    }
                    
                    
                    ar = AMAX - acceleration.norm();
                    if (0 != Aijv.norm()) {
                        acceleration = acceleration + Aijv.normalize()*mini(ar, Aijv.norm());
                    }
                    else
                    {
                        acceleration = acceleration + Aijv;
                    }
                    
                    
                    ar = AMAX - acceleration.norm();
                    if (0 != Aijc.norm()) {
                        acceleration = acceleration + mini(ar, Aijc.norm())*Aijc.normalize();
                    }
                    else
                    {
                        acceleration = acceleration +Aijc;
                    }
                    
                    double KD = 0;
                    if (Xij.norm()<R1) {
                        KD = 1;
                    }else if(Xij.norm()>R2)
                    {
                        KD = 0;
                    }else
                    {
                        KD = (R2-Xij.norm())*1.0/(R2 - R1);
                    }
                    statesA[i+BOIDNUMBER] =statesA[i+BOIDNUMBER] + KD*acceleration;
                    
                }
            }
            
        }
        
        
        
    }
}


//*******************  calculate states *******************//

void springForceUpdate( std::vector<Vector3d>& states, unsigned int pNum)
{
    
    //calculate springForce and dampForce
    for (int k=0; k<STRUTNUMBER; k++) {
        int i = springMesh.struts[k].vertexIndices[0];
        int j = springMesh.struts[k].vertexIndices[1];
        
        //std::cout<<"i: "<<i<<"j: "<<j<<"\n";
        
        
        double Kij = springMesh.struts[k].k;
        double Dij = springMesh.struts[k].d;
        double L0ij = springMesh.struts[k].l;
        
        //std::cout<<"k: "<<Kij<<"d: "<<Dij<<"L0 "<<L0ij<<"\n";
        
        Vector3d Xij = states[j] - states[i];
        //std::cout<<Xij<<"\n";
        
        if (Xij.norm() > ERRORTHRESH) {
            Vector3d XijN = Xij.normalize();
            double Lij = Xij.norm();
            
            Vector3d Fis(0,0,0);
            Vector3d Fjs(0,0,0);
            Vector3d Fid(0,0,0);
            Vector3d Fjd(0,0,0);
            
            Fis = Kij*(Lij - L0ij)* XijN;
            Fjs = -1*Fis;
           
            Fid = Dij*((states[j+pNum] - states[i+pNum])*XijN)*XijN;
            Fjd = -1*Fid;
            
            springMesh.points[i].force = springMesh.points[i].force + Fis + Fid;
            springMesh.points[j].force = springMesh.points[j].force + Fjs + Fjd;
        }
    }
    
}

//system Dynamical Function
std::vector<Vector3d> sysDynaFunc(std::vector<Vector3d>& states, const unsigned int pNum , const calcuState state)
{
    std::vector<Vector3d> statesA(pNum*2);
    
    if (state == SPRINGCAL) {
        springForceUpdate(states,pNum);
        for (int i =0; i<pNum; i++) {
            statesA[i] = states[i+pNum];
            if (i%ROWS == 0 || i<ROWS  || i>(pNum-ROWS-1) || i%ROWS == (ROWS-1)) {
                
                statesA[i+pNum] = Vector3d(0,0,0);
                statesA[i] = Vector3d(0,0,0);
                states[i+pNum] = Vector3d(0,0,0);
            }
            
            else
            {
                statesA[i+pNum] = (springMesh.points[i].force)*(1.0/springMesh.points[i].mass);
            }
            //clean the force every loop !!!!
            springMesh.points[i].force = Vector3d(0,0,0);
        }
    }
    
    if (state == PARTICLECAL) {
        for (int i=0; i<pNum; i++) {
            statesA[i] = states[i+pNum];
 
            Vector3d Xai = states[i] - circleCenterPosition;
            double distance  = Xai.norm();
            Vector3d particleAcceleration((-0.001 * (distance*distance))*Xai);
            
            //Vector3d particleAcceleration1((-10 * 1.0/(distance*distance))*Xai);

            statesA[i+pNum] = particleAcceleration;
            
        }
        flockAcceleration(states, statesA, pNum);
        
    }
    
    if (state == FIREWORKCAL) {
        for (int i=0; i<pNum; i++) {
            statesA[i] = states[i+pNum];
            
            
            Vector3d Xai = states[i] - circleCenterPosition;
            double distance  = Xai.norm();
            Vector3d particleAcceleration((-0.001 * (distance*distance))*Xai);

            statesA[i+pNum] = states[i+pNum]*-0.1 + particleAcceleration + gravity;
        }
    }
    
    
    return statesA;
}

//numerical integration function
void statesNumInt(std::vector<Vector3d>& states, std::vector<Vector3d>& statesNew, std::vector<Vector3d>& statesA, const unsigned int pNum, const calcuState state)
{
    statesA = sysDynaFunc(states,pNum, state);
    for (int i = 0; i<pNum*2; i++) {
        statesNew[i] = states[i] + statesA[i]*hStep;
    }
    
    
}

void statesNumIntRK4(std::vector<Vector3d>& states, std::vector<Vector3d>& statesNew, std::vector<Vector3d> &statesA, const unsigned int pNum, const calcuState state)
{
    std::vector<Vector3d> statesK1(pNum*2);
    std::vector<Vector3d> statesK2(pNum*2);
    std::vector<Vector3d> statesK3(pNum*2);
    std::vector<Vector3d> statesK4(pNum*2);
    
    std::vector<Vector3d> statesTMP(pNum*2);
    
    for (int i = 0; i<pNum*2; i++) {
        statesTMP[i] = states[i];
    }
    
    
    //1st
    statesK1 = sysDynaFunc(statesTMP, pNum, state);
    for (int i = 0; i<pNum*2; i++) {
        statesTMP[i] = states[i] + statesK1[i]*hStep;
    }

    //2ed
    statesK2 = sysDynaFunc(statesTMP, pNum, state);
    for (int i = 0; i<pNum*2; i++) {
        statesTMP[i] = states[i] + statesK2[i]*0.5*hStep;
    }
    
    //3rd
    statesK3 = sysDynaFunc(statesTMP, pNum, state);
    for (int i = 0; i<pNum*2; i++) {
        statesTMP[i] = states[i] + statesK3[i]*0.5*hStep;
    }
    
    //4st
    statesK4 = sysDynaFunc(statesTMP, pNum, state);
    
    for (int i = 0; i<pNum*2; i++) {
        statesA[i] = (statesK1[i] + 2*statesK2[i] + 2*statesK3[i] + statesK4[i])*(1.0/6);
        statesNew[i] = states[i] + statesA[i]*hStep;
    }

}

//************************************dectect collision part*********************************************************//
bool detectTriangleCollision(const Vector3d& particlePosition, const Vector3d& particleVelocity,  Vector3d& particlePositionNew  , Vector3d&  particleVelocityNew, const Vector3d& P0, const Vector3d& P1, const Vector3d& P2, int pIndex)
{
    

    Vector3d vertexPosition =particlePosition;
    vertexPosition.y = vertexPosition.y-particleSize;
    
    Vector3d vertexPositionNew =particlePositionNew;
    vertexPositionNew.y = vertexPositionNew.y-particleSize;
    
    
    
        Vector3d vn =  (P2 - P1)%(P1 - P0);
        double vnnorm = vn.norm();
        Vector3d n = vn.normalize();
        double a = (vertexPosition - P0)*n;
        double b = (vertexPositionNew - P0)*n;
        
        Vector3d D = (vertexPositionNew - vertexPosition).normalize();
        
        Vector3d hitPosition;
        if (a*b<0) {
            double s = a/(a-b);
            
            hitPosition = vertexPosition + (vertexPositionNew - vertexPosition)*(s*hStep);
            
            Vector3d e1 = P1-P0;
            Vector3d e2 = P2-P0;
            
            Vector3d P = D % e2;
            double det = e1 * P;
            
            if(det > -EPSILON && det < EPSILON) return false;
            double inv_det = 1.f / det;
            
            Vector3d T = particlePosition - P0;
            double u = (T*P)*inv_det;
            
            if(u < 0.f || u > 1.f) return false;
            
            
            Vector3d Q = T % e1;
            double v = (D * Q) * inv_det;
            
            if(v < 0.f || u + v  > 1.f) return false;
            
            double t = (e2 * Q) * inv_det;
            
            if(t > EPSILON) { //ray intersection
                double d1 = (vertexPosition - hitPosition)*n;
                double d2 = (vertexPositionNew - hitPosition)*n;
                vertexPositionNew = vertexPositionNew - 1.7*d2*n;
                vertexPositionNew.y = vertexPositionNew.y+particleSize;
                particlePositionNew = vertexPositionNew;
                
                Vector3d vn = (particleVelocity*n)*n;
                Vector3d vt = particleVelocity - vn;
                
                particleVelocityNew = -0.2*vn + 0.6*vt;
                
                
                
                return true;
            }
            
        }
    
    
    return false;
}


//*********************************************************************************************//

void springUpdate(const unsigned int pNum, const calcuState state)
{
    
    std::vector<Vector3d>vertexStates(pNum*2);
    std::vector<Vector3d>vertexStatesNew(pNum*2);
    std::vector<Vector3d>vertexStatesA(pNum*2);
    
    
    for (int i=0; i<pNum; i++) {
        vertexStates[i] = springMesh.points[i].xposition;
        vertexStates[i+pNum] = springMesh.points[i].velocity;
    }
    
    //statesNumInt(vertexStates, vertexStatesNew, vertexStatesA, pNum);
    statesNumIntRK4(vertexStates, vertexStatesNew, vertexStatesA, pNum, state);
    
    for (int i=0; i<pNum; i++) {
        springMesh.points[i].xposition = vertexStatesNew[i];
        springMesh.points[i].velocity = vertexStatesNew[i+pNum];
    }
    
}





bool checkInTriangleFace(const Vector3d& p1, const Vector3d& v)
{
    float thresh = 3;
    if (v.x>(p1.x+thresh) || v.x<(p1.x-thresh)) return false;
    if (v.y>(p1.y+thresh) || v.y<(p1.y-thresh)) return false;
    if (v.z>(p1.z+thresh) || v.z<(p1.z-thresh)) return false;
    return true;
}

void particlesUpdate(const unsigned int pNum, const calcuState state)
{
    std::vector<Vector3d>vertexStates(pNum*2);
    std::vector<Vector3d>vertexStatesNew(pNum*2);
    std::vector<Vector3d>vertexStatesA(pNum*2);
    for (int i=0; i<pNum; i++) {
        vertexStates[i] = particles[i].getPosition();
        vertexStates[i+pNum] = particles[i].getVelocity();
    }
    statesNumInt(vertexStates, vertexStatesNew, vertexStatesA, pNum, state);
    
    for (unsigned int f=0; f<FACENUMBER; f++) {
        Vector3d P0 = springMesh.points[springMesh.faces[f].pointIndices[0]].xposition;
        Vector3d P1 = springMesh.points[springMesh.faces[f].pointIndices[1]].xposition;
        Vector3d P2 = springMesh.points[springMesh.faces[f].pointIndices[2]].xposition;

        for (int i=0; i<pNum; i++) {
            if (vertexStatesNew[i].y < -1  && checkInTriangleFace(P1, particles[i].getPosition())) {
                if(detectTriangleCollision(vertexStates[i], vertexStates[i+pNum], vertexStatesNew[i], vertexStatesNew[i+pNum],P0,P1,P2,i))
                {
                    int V0 = springMesh.faces[f].pointIndices[0];
                    int V1 = springMesh.faces[f].pointIndices[1];
                    int V2 = springMesh.faces[f].pointIndices[2];
                    
                    if (i%ROWS == 0 || i<ROWS  || i>(pNum-ROWS-1) || i%ROWS == (ROWS-1)) {
                        springMesh.points[V0].velocity =  Vector3d(0,0,0);
                        springMesh.points[V1].velocity =  Vector3d(0,0,0);
                        springMesh.points[V2].velocity =  Vector3d(0,0,0);
                    }
                    else
                    {
                        springMesh.points[V0].velocity =  vertexStates[i+pNum]*0.5;
                        springMesh.points[V1].velocity =  vertexStates[i+pNum]*0.5;
                        springMesh.points[V2].velocity =  vertexStates[i+pNum]*0.5;
                    }
                    
                    //particles[i].setStopSign(true);
                    
                }
            }
        }
    }
    
    
    unsigned int flockCounter = 0;
    Vector3d flockCenter(0,0,0);
    for (int i=0; i<pNum; i++)
    {
        if (particles[i].getStopSign() == true) {
            flockCounter++;
            flockCenter = flockCenter + particles[i].getPosition();
        }

    }
    flockCenter = flockCenter / flockCounter;
    
    
    for (int i=0; i<pNum; i++) {
        particles[i].setPosition(vertexStatesNew[i]);
        particles[i].setVelocity(vertexStatesNew[i+pNum]);
        if (particles[i].getPosition().y<-10) {
            particles[i].setPosition( Vector3d(-10, gauss(7, 1, 1), gauss(0, 2, 1) ) );
            particles[i].setVelocity( Vector3d(gauss(2.5, 0.5, 1), 0, 0) );
        }
        
        
        
        if ( particles[i].getVelocity().norm() <1  || particles[i].getVelocity().norm() >8) {
             particles[i].setStopSign(true);
        }
        
        
        
        if (particles[i].getName() == "trace") {
            
            traceGenerator(particles[i].getPosition(), 20);
            if ((particles[i].getPosition() - flockCenter).norm() < 1 ) {
                explode = true;
                particles[i].setName("chaos");
                fireWorkGenerator(particles[i].getPosition(), 2000);
            }
            
        }
        
    }
    
    
    if (explode) {
        
        for (int i=0; i<pNum; i++) {
            if ( particles[i].getStopSign()) {
                
                if (particles[i].getName() != "trace") {
                    particles[i].setStopSign(false);
                }
                Vector4d color(rand()*0.8/RAND_MAX, rand()*0.9/RAND_MAX, rand()*0.7/RAND_MAX, rand()*1.0/RAND_MAX);
                Vector3d tmp = particles[i].getVelocity();
                particles[i].setVelocity(Vector3d(tmp.x+gauss(0, 2, 1), tmp.y+gauss(0, 2, 1), tmp.z+gauss(0, 2, 1)));
                particles[i].setColor(color);
                
            }
            
            
        }
        explode = false;
    }
    
    
}


void firewordUpdate(const unsigned int pNum, const calcuState state)
{
    
    std::vector<Vector3d>vertexStates(pNum*2);
    std::vector<Vector3d>vertexStatesNew(pNum*2);
    std::vector<Vector3d>vertexStatesA(pNum*2);
    for (int i=0; i<pNum; i++) {
        vertexStates[i] = fireWorks[i].getPosition();
        vertexStates[i+pNum] = fireWorks[i].getVelocity();
    }
    
   statesNumInt(vertexStates, vertexStatesNew, vertexStatesA, pNum, state);
    
   for (int i=0; i<pNum; i++) {
       fireWorks[i].setPosition(vertexStatesNew[i]);
       fireWorks[i].setVelocity(vertexStatesNew[i+pNum]);
   }
    
    
    std::vector<Particle>::iterator ptr = fireWorks.begin();
    while (ptr != fireWorks.end()) {
        
        float speedLimit = 0.2;
        if (ptr->getName() == "firework") {
            speedLimit = 2;
        }
        else
        {
            speedLimit = 0.5;
        }

        if ((ptr->getVelocity()).norm() < speedLimit ) {
            ptr = fireWorks.erase(ptr);
        }
        else
            ptr++;
        
    }
    
}



//*******************************************************//


void timeProc(int id)
{
    if (id == 1) {
        
        springUpdate(VERTEXNUMBER, SPRINGCAL);
        
        if (particles.size()<ParticleNum) {
            particlesGenerator(1);
        }
        particlesUpdate(particles.size(), PARTICLECAL);
    
        
        if (fireWorks.size() > 0) {
            firewordUpdate(fireWorks.size(), FIREWORKCAL);
        }
        

        
        glutPostRedisplay();
        if (resetSign == false) {
            glutTimerFunc(automaticSpeed, timeProc, 1);
        }
    }
    
}

void handleKey(unsigned char key, int x, int y){
    unsigned int t = gauss(VERTEXNUMBER/2, VERTEXNUMBER/6, 1);
    unsigned int n = rand()%VERTEXNUMBER;
    
    switch(key){
        case 'a':
        case 'A':
            
            resetSign = false;
            glutTimerFunc(automaticSpeed, timeProc, 1);
            break;
        case 's':
        case 'S':
            resetSign = true;
            break;
            
        case 'c':
        case 'C':
            clothShow = !clothShow;
            glutPostRedisplay();
            break;
            
        case 'w':
        case 'W':
            windCount++;
            break;
            
        case 'r':			// R -- toggle between recording or not
        case 'R':
            Recording = !Recording;
            glutPostRedisplay();
            break;
        
        case 'k':
        case 'K':
            if (t<VERTEXNUMBER && t>0 ) {
                if (!(t%ROWS == 0 || t<ROWS  || t>(VERTEXNUMBER-ROWS-1) || t%ROWS == (ROWS-1))) {
                    springMesh.points[t].velocity.y = 20;
                    
                }
            }
            break;
        case 'g':
        case 'G':
            fireWorkGenerator(circleCenterPosition, 100);
            break;
            
        case 'z':
        case 'Z':
            
            particles.push_back( Particle(springMesh.points[n].xposition, Vector3d(gauss(2.5, 0.5, 1),2,0), Vector3d(0,0,0), Vector4d(0,0.5,0.2,1), 1, 0, particleSize, true, "trace") );
            break;
            
        case 'm':
        case 'M':
            showGrid = !showGrid;
            break;
            
        case 'q':       // q - quit
        case 'Q':
        case 27:        // esc - quit
            exit(0);
            
        default:        // not a valid key -- just ignore it
            return;
    }
}


void doDisplay(){
    myDisplay();
    glutSwapBuffers();
    
    // recordimage() method needs a pointer to the routine that does all of the drawing
    if(Recording)
        framegrabber.recordimage(myDisplay);
}


void init() {
    //LoadParameters(ParamFilename);
    // set up camera
    // parameters are eye point, aim point, up vector
    camera = new Camera(Vector3d(0, 0, 1), Vector3d(0, 0, 0), Vector3d(0, 1, 0));
    
    // black background for window
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glShadeModel(GL_FLAT);
    //glDepthRange(0, 1);
    
    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    //Setting lights
    
    GLfloat light0_position[] = {15.0f, 15.0f, 10.0f, 1.0f}; //spot light
    GLfloat light1_position[] = {-15.0f, 15.0f,  10.0f, 1.0f};
    GLfloat light2_position[] = {0.0f, 15.0f, -15.0f, 1.0f};
    GLfloat light_ambient[]  = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat light_diffuse[]  = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat light_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
     
    
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    
    
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    
    glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
    glLightfv(GL_LIGHT2, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT2, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular);
    
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    

    
    particles.reserve(ParticleNum);
    fireWorks.reserve(1000000);
    srand (time(NULL));
    
    texture.loadTGA("/Users/yuyaolong/Desktop/finalProj/clemson.tga");
    glEnable(GL_CULL_FACE);
    
    
    
}


int main(int argc, char *argv[])
{
    
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutCreateWindow("Particle");
    
    init();
    
    glutDisplayFunc(doDisplay);
    glutMouseFunc(mouseEventHandler);
    glutMotionFunc(motionEventHandler);
    glutKeyboardFunc(handleKey);
    
    glutMainLoop();
    
    
    return 0;
}