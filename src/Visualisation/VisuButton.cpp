
#include "VisuButton.h"

VisuButton::VisuButton(std::string _caption,void (*_onClick)(),int _x,int _y, int _w,int _h)
{
	caption=_caption;
	onClick=_onClick;
	x=_x;y=_y;
	h=_h;w=_w;
	
	clicked=false;
	mouse_in=false;
}


void VisuButton::checkForAtivity(int button, int state, int _x, int _y)
{
	mouse_in=(_x>x && _y>y && _x<x+w && _y<y+h );
	if(state == GLUT_DOWN)
	{
		if(mouse_in)
			clicked=true;
	}
	if(state == GLUT_UP)
	{
		if(mouse_in && clicked)
		{
			(*onClick)();
		}
		clicked=false;
	}
}
void VisuButton::checkForAtivityMotion(int _x,int _y)
{
	mouse_in=(_x>x && _y>y && _x<x+w && _y<y+h );

}

void VisuButton::draw()
{
	set2DGLProjection();
	
	if(clicked && mouse_in)
		glColor3f(0.8,0.5,0);
	else
		glColor3f(1.0,0.7,0);
	glBegin(GL_QUADS);
	glVertex2f(x, y); glVertex2f(x+w, y); glVertex2f(x+w, y+h); glVertex2f(x, y+h);
	glEnd()	;
	glColor3f(0,0.,0);
	
	glBegin(GL_LINE_LOOP);
	glVertex2f(x, y); glVertex2f(x+w, y); glVertex2f(x+w, y+h); glVertex2f(x, y+h);
	glEnd()	;
	
	drawText(x+7, y+7, caption.c_str(), 1,1,1);

	glColor3f(1,1,1);
	unset2DGLProjection();	
}