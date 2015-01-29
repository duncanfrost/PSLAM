
#include "SimuWindow.h"


namespace SimuWindowStuff
{
//pose from which map is displayed
HomogeneousMatrix22 currentCamPose;
void keyboard_process(unsigned char key, int x, int y);
void mouse_process(int button, int state, int x, int y);
void motion_process(int x, int y);

float z_near_map_win=0.2;

int mouse_old_x, mouse_old_y;
float velocity_rotation=0.05;
float velocity_translation=0.03;
}
  
using namespace SimuWindowStuff;


SimuWindow::SimuWindow(std::string _title,std::string _objBaseName,std::string _objFile,int _width,int _height,Camera *_myCamera)
{
	title=_title;
	objBaseName=_objBaseName;
	objFile=_objFile;
	width=_width;
	height=_height;
	
	myCamera=_myCamera;
	show_texture=true;
	current_img.create(height,width, CV_8UC3);
	current_img_BW.create(height,width, CV_8UC1);
	current_img_Z.create(height, width, CV_32FC1);
}


void SimuWindow::prepare_draw()
{
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
 	set3DGLProjection(myCamera,currentCamPose);
        glCallList(objId);       //and display it
		
	displayedPose=currentCamPose;
	current_img = Capture_Image();
	cvtColor(current_img,current_img_BW,CV_RGB2GRAY);
	current_img_Z= Capture_realZ();
	
	glFlush();
	glutSwapBuffers();
	
	//get gl ready for add draw function from main
	glDisable(GL_TEXTURE_2D);
	set2DGLProjection();
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);

	
	/*glColor3f(1,1,1);
	unset2DGLProjection();	
	glEnable(GL_TEXTURE_2D);*/
}


GLfloat light_ambient[] = {1.0, 1.0, 1.0, 0.3};
GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 0.7};  /* Red diffuse light. */
GLfloat light_position[] = {10.0, 10, 10.0, 0.0};  /* Infinite light location. */

void SimuWindow::setLight()
{
	if(show_texture)
	{
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHTING);
		GLfloat global_ambient[] = {1.0, 1.0, 1.0, 1.0f };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	}
	else
	{
		GLfloat light_ambient[] = {1.0, 1.0, 1.0, 0.3};
		
		GLfloat light_diffuse[] = {1.0, 0.0, 0.5, 0.7};  
		GLfloat light_position[] = {-10.0, 10, 10.0, 0.0};  
		glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient );
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
		glLightfv(GL_LIGHT0, GL_POSITION, light_position);
		glEnable(GL_LIGHT0);
		GLfloat light_diffuse1[] = {0.5, 1.0, 0.0, 0.7}; 
		GLfloat light_position1[] = {10.0, -10, 10.0, 0.0};  
		//glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1 );
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
		glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
		GLfloat light_diffuse2[] = {0.9, 1.0, 0.2, 0.7}; 
		GLfloat light_position2[] = {10.0, 10, -10.0, 0.0};  
		glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
		glLightfv(GL_LIGHT2, GL_POSITION, light_position2);
		
		
		glEnable(GL_LIGHT1);
		glEnable(GL_LIGHT2);
		glEnable(GL_LIGHTING);
		
		GLfloat mat_specular[] = { 0.0, 0.0, 0.0, 0.0 };
		GLfloat mat_shininess[] = { 100.0 };
		GLfloat mat_surface[] = { 0.2, 0.0, 0.0, 0.9 };
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_surface);
	}
}
void SimuWindow::CreateWindow()
{
#ifdef AMOVERBOSE
	std::cout<<"SimuWindow::CreateWindow()"<<std::endl;
#endif
	glutInitWindowSize(width, height);	
	window = glutCreateWindow(title.c_str());
	
	 glClearColor(0.5,0.5,0.5,1.0);
	
	//position camera in interesting place in map (on ground)
	 currentCamPose.SetIdentity();
	currentCamPose.TranslateX(-3);
	currentCamPose.TranslateY(-2.1);
	
	currentCamPose.TranslateZ(0.1);
	currentCamPose.RotateX(-90);
	currentCamPose.RotateY(90);
	
	//config planar scene
	//currentCamPose.RotateY(90);
	
	std::cout<<currentCamPose<<std::endl;
       	set3DGLProjection(myCamera,currentCamPose);
 	glEnable (GL_BLEND);
        glEnable(GL_DEPTH_TEST);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        objId=obj.load(objBaseName,objFile.c_str());      //no transparent texture
       // objId=obj.load(WIMU_WORLD_DIR,"city.obj");      //no transparent texture
	setLight();
}

void SimuWindow::setEvents()
{
	//get id in static array
	glutKeyboardFunc(keyboard_process);
	glutMouseFunc(mouse_process);
	glutMotionFunc(motion_process);         
}



void SimuWindowStuff::keyboard_process(unsigned char key, int x, int y)
{
	int window=glutGetWindow();
	int idS=VisuWindowNs::windowIDtoStatic(window);
	
	 if(VisuWindowNs::mImageWinStatic[idS].processKeys_fct)
		(*VisuWindowNs::mImageWinStatic[idS].processKeys_fct)(key,x,y);
	
	switch(key) {
		case 52://left
			currentCamPose.TranslateX(-velocity_translation);
			break;
		case 54://right
			currentCamPose.TranslateX(velocity_translation);
			break;
		case 56:
			currentCamPose.TranslateZ(velocity_translation);
			break;
		case 53:
			currentCamPose.TranslateZ(-velocity_translation);
			break;
		case 57:
			currentCamPose.TranslateY(-velocity_translation);
			break;
		case 51:
			currentCamPose.TranslateY(velocity_translation);
			break;
	}


	glutPostRedisplay();
}

void SimuWindow::set_velocity_translation(float _f){velocity_translation=_f;}

void SimuWindowStuff::mouse_process(int button, int state, int x, int y)
{
	if(state == GLUT_DOWN)
	{
	    mouse_old_x=x;
	    mouse_old_y=y;
	}
	
	int window=glutGetWindow();
	int idS=VisuWindowNs::windowIDtoStatic(window);
	    
	//do normal thing to do in window => check button action
	for(int b=0;b<VisuWindowNs::mButtons[idS].size();b++)
	      VisuWindowNs:: mButtons[idS][b].checkForAtivity(button,state,x,y);
	
	 if(VisuWindowNs::mImageWinStatic[idS].processMouse_fct)
		(*VisuWindowNs::mImageWinStatic[idS].processMouse_fct)(button,state,x,y);
}

void SimuWindowStuff::motion_process(int x, int y)
{
	currentCamPose.RotateY(velocity_rotation*(x-mouse_old_x));
	currentCamPose.RotateX(-velocity_rotation*(y-mouse_old_y));

	mouse_old_x = x;
	mouse_old_y = y;
	
 	int window=glutGetWindow();
	int idS=VisuWindowNs::windowIDtoStatic(window);
	    
	//do normal thing to do in window => check button action
	for(int b=0;b<VisuWindowNs::mButtons[idS].size();b++)
	      VisuWindowNs:: mButtons[idS][b].checkForAtivityMotion(x,y);
	
	 if(VisuWindowNs::mImageWinStatic[idS].processMouseActiveMotion_fct)
		(*VisuWindowNs::mImageWinStatic[idS].processMouseActiveMotion_fct)(x,y);   
    
}
void SimuWindow::moveCamera(HomogeneousMatrix22 _cHc2)
{
	currentCamPose=_cHc2*currentCamPose;
}
void SimuWindow::setCameraPose(HomogeneousMatrix22 _cHc2)
{
	currentCamPose=_cHc2;	
}

	