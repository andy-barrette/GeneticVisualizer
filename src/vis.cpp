//g++ -o vis vis.cpp glut32.lib -lopengl32 -lglu32  //use this format to compile
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <windows.h>
#include "GL/glut.h"
#include <math.h>
#include <vector>
using namespace std;

/*
short PAT_SIZE=2;
short BLOCK_SIZE=82;
short BLOCK_NUM=5;
float DECAY_FACTOR=.7;
float START_STEP=.1;
short DOT_SIZE=1;
*/
#define PAT_SIZE 2
#define BLOCK_SIZE 82
#define BLOCK_NUM 5
#define DECAY_FACTOR .7
#define START_STEP .1
#define DOT_SIZE 1


GLint winw=750,winh=750;
char current_key;

class arraytype{
    public:
        
    float a[BLOCK_NUM][BLOCK_NUM][BLOCK_SIZE][BLOCK_SIZE];
    
    void set(int i,int j,int k,int l,float val){
        a[i][j][k][l]=val;
    }
    
    float get(int i,int j,int k,int l){
        return a[i][j][k][l];
    }
    
    void copyto(arraytype* v){
        int i,j,k,l;
        for(i=0;i<BLOCK_NUM;i++)
            for(j=0;j<BLOCK_NUM;j++)
                for(k=0;k<BLOCK_SIZE;k++)
                    for(l=0;l<BLOCK_SIZE;l++)
                        v->a[i][j][k][l]=a[i][j][k][l];
    }
};

class algorithmtype{
    vector<float> genes;
    arraytype arrays;
    short curgene[2];
    float curstep[2];
    
    public:
    
    void resetgenes(){//reset genes
        srand(time(NULL));
        genes.clear();
        float genenum=(2*PAT_SIZE+1)*(2*PAT_SIZE+1);
        int i;
        for(i=0;i<genenum;i++){
            genes.push_back(round(pow(fmod(rand(),1000)/1000,2))*fmod(rand(),1000)/(genenum*1000.0));
        }
    }
    /*
    void genegain(float factor){
        int i;
        for(i=0;i<genes.size();i++)genes[i]*=factor;
    }
    */
    void resetarrays(){//reset ararys
        int i,j,k,l;
        for(i=0;i<BLOCK_NUM;i++)
            for(j=0;j<BLOCK_NUM;j++)
                for(k=0;k<BLOCK_SIZE;k++)
                    for(l=0;l<BLOCK_SIZE;l++){
                        if(k<BLOCK_SIZE*pow(1.02,DOT_SIZE)/2.0 && k>BLOCK_SIZE*pow(.98,DOT_SIZE)/2.0 && l<BLOCK_SIZE*pow(1.02,DOT_SIZE)/2.0 && l>BLOCK_SIZE*pow(.98,DOT_SIZE)/2.0)
                            arrays.a[i][j][k][l]=1;
                        else
                            arrays.a[i][j][k][l]=0;
                    }
    }
    
    algorithmtype(){
        curgene[0]=0;
        curgene[1]=1;
        curstep[0]=START_STEP;
        curstep[1]=START_STEP;
        resetgenes();
        resetarrays();
    }
    
    void init(){
        curgene[0]=0;
        curgene[1]=1;
        curstep[0]=START_STEP;
        curstep[1]=START_STEP;
        resetgenes();
        resetarrays();
    }

    void timeevolve(){
        int i,j,k,l;
        arraytype temparrays;
        arrays.copyto(&temparrays);
        for(i=0;i<BLOCK_NUM;i++)
            for(j=0;j<BLOCK_NUM;j++)
                for(k=PAT_SIZE;k<BLOCK_SIZE-PAT_SIZE;k++)
                    for(l=PAT_SIZE;l<BLOCK_SIZE-PAT_SIZE;l++){
                        int a,b;
                        for(a=-PAT_SIZE;a<=PAT_SIZE;a++)
                            for(b=-PAT_SIZE;b<=PAT_SIZE;b++){
                                if((a+PAT_SIZE)+(b+PAT_SIZE)*(2*PAT_SIZE+1)==curgene[0])
                                    temparrays.a[i][j][k][l]+=(genes[curgene[0]]+(curstep[0]*(i-floor(BLOCK_NUM/2.0))))*arrays.a[i][j][k+a][l+b];
                                else if((a+PAT_SIZE)+(b+PAT_SIZE)*(2*PAT_SIZE+1)==curgene[1])
                                    temparrays.a[i][j][k][l]+=(genes[curgene[1]]+(curstep[1]*(j-floor(BLOCK_NUM/2.0))))*arrays.a[i][j][k+a][l+b];
                                else 
                                    temparrays.a[i][j][k][l]+=genes[(a+PAT_SIZE)+(b+PAT_SIZE)*(2*PAT_SIZE+1)]*arrays.a[i][j][k+a][l+b];
                            }
                        temparrays.a[i][j][k][l]*=DECAY_FACTOR;
                    }
        temparrays.copyto(&arrays);
    }
    
    void disp(){
        int i,j,k,l;
        float pixelsize=(float)winw/(BLOCK_NUM*(BLOCK_SIZE-(2*PAT_SIZE)));
        for(i=0;i<BLOCK_NUM;i++)
            for(j=0;j<BLOCK_NUM;j++)
                for(k=PAT_SIZE;k<BLOCK_SIZE-PAT_SIZE;k++)
                    for(l=PAT_SIZE;l<BLOCK_SIZE-PAT_SIZE;l++){
                        glColor3f(arrays.a[i][j][k][l],arrays.a[i][j][k][l],arrays.a[i][j][k][l]);
                        glBegin(GL_QUADS);
                            glVertex2f((winw*i/(float)BLOCK_NUM)+k*pixelsize,(winh*j/(float)BLOCK_NUM)+l*pixelsize);
                            glVertex2f((winw*i/(float)BLOCK_NUM)+(k+1)*pixelsize,(winh*j/(float)BLOCK_NUM)+l*pixelsize);
                            glVertex2f((winw*i/(float)BLOCK_NUM)+(k+1)*pixelsize,(winh*j/(float)BLOCK_NUM)+(l+1)*pixelsize);
                            glVertex2f((winw*i/(float)BLOCK_NUM)+k*pixelsize,(winh*j/(float)BLOCK_NUM)+(l+1)*pixelsize);
                        glEnd();
                    }
    }
    
    void applyselection(short xblock,short yblock){
        if(xblock==0){
            curgene[0]=(curgene[0]+2)%genes.size();
            curstep[0]=START_STEP;
        }
        else{
            genes[curgene[0]]+=xblock*curstep[0];
            curstep[0]*=abs(xblock)*.75;
        }
        if(yblock==0){
            curgene[1]=(curgene[1]+2)%genes.size();
            curstep[1]=START_STEP;
        }
        else{
            genes[curgene[1]]+=yblock*curstep[1];
            curstep[1]*=abs(yblock)*.75;
        }
    }
}p;

void resize(int x,int y){
    if(x>0)winw=x;
    if(y>0)winh=y;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,winw,winh,0);
}

void disp() 
{
	glClear(GL_COLOR_BUFFER_BIT);		     // Clear Screen and Depth Buffer
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	glViewport(0,0,winw,winh);
	p.disp();
	p.timeevolve();
	glutSwapBuffers();
}

void init() 
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,winw,winh,0);	
	p.init();
}

/*
void specialkeyfunc(int key, int x, int y){
    if(key==GLUT_KEY_UP)
        switch(current_key){
            case 'r': BLOCK_SIZE+=10; p.resetarrays(); break;
            case 'b': BLOCK_NUM+=2; p.resetarrays(); break;
            case 'a': genegain(1.05); p.resetarrays(); break;
            case 'd': DECAY_FACTOR*=1.05; p.resetarrays(); break;
            case 'p': DOT_SIZE*=1.1; p.resetarrays(); break;
            case 'n': PAT_SIZE++; p.resetgenes(); p.resetarrays.(); break;
        }
    else if(key==GLUT_KEY_DOWN)
        switch(current_key){
            case 'r': if(BLOCK_SIZE>11)BLOCK_SIZE-=10; p.resetarrays(); break;
            case 'b': if(BLOCK_NUM>=5)BLOCK_NUM-=2; p.resetarrays(); break;
            case 'a': genegain(.95); p.resetarrays(); break;
            case 'd': DECAY_FACTOR*=.95; p.resetarrays(); break;
            case 'p': DOT_SIZE*=.9; p.resetarrays(); break;
            case 'n': if(PAT_SIZE>1)PAT_SIZE--; p.resetgenes(); p.resetarrays.(); break;
        }
}
            
void keyupfunc(unsigned char key, int x, int y){
    if(key==current_key)current_key=NULL;
}

void keyfunc(unsigned char key, int x, int y){
    if(key=='0')p.resetarrays();
    else current_key=key;
}
*/
void mouse(int but,int state,int x,int y){
    if(but==GLUT_LEFT_BUTTON && state==GLUT_DOWN){
        p.applyselection((short)(x*BLOCK_NUM/winw)-floor(BLOCK_NUM/2.0),(short)(y*BLOCK_NUM/winh)-floor(BLOCK_NUM/2.0));
        p.resetarrays();
    }
}

int main(int argc, char **argv) 
{
	glutInit(&argc, argv);                                      // GLUT initialization
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);  // Display Mode
	glutInitWindowSize(winw,winh);					// set window size
	glutCreateWindow("evolving visualizer");								// create Window
	glutDisplayFunc(disp);									// register Display Function
	glutIdleFunc(disp);						
    glutMouseFunc(mouse);
	glutReshapeFunc(resize);
	//glutKeyboardFunc(keyfunc);
	//glutSpecialFunc(specialkeyfunc);
	//glutKeyUpFunc(keyupfunc);
	init();
	glutMainLoop();												// run GLUT mainloop
	return 0;
}
