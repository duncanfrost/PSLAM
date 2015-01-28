
#include "PlottingWindow.h"

PlottingWindow::PlottingWindow(std::string _title,int _width,int _height,int _nb_plots,int _timeLine,bool _incTimeWhenDisplay)
{
	title=_title;
	width=_width;
	height=_height;
	nb_plots=_nb_plots;
	timeLine=_timeLine;
	incTimeWhenDisplay=_incTimeWhenDisplay;
	current_time=0;
	
	//init plots with one value to display by default
	mPlotStructs.resize(nb_plots);
	for(int i=0;i<nb_plots;i++)
		mPlotStructs[i].init(1,timeLine);
	
	//init random colors
	for(int i=0;i<MAX_VAL_PER_PLOTS;i++)
	{
		Vector3f randColCamera = Vector3f::Random(3);randColCamera=randColCamera/randColCamera.maxCoeff();
		colPlots.push_back(randColCamera);
	}

}
void PlottingWindow::setValuesPerPlot(int _p,int _n_val)
{
	mPlotStructs[_p].init(_n_val,timeLine);
}
void PlottingWindow::setVal(int _p,int _ival,float _val)
{
	mPlotStructs[_p].setVal(current_time,_ival,_val);
}

void PlottingWindow::prepare_draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	set2DGLProjection();
	
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPointSize(3.0);

	//glColor3f(0.5,1,0.7);	
	//glBegin(GL_POINTS);	glVertex2f(10,10);	glEnd();
	
	int margin=5;
	int height_plot=height/nb_plots-2*margin;//put margin on top and bottom of each plot
	int width_plot=width-2*margin;//put margin on top and bottom of each plot
	for(int p=0;p<nb_plots;p++)
	{
		plotStruct &currentPlot=mPlotStructs[p];
		
		int y1=p*(height/nb_plots)+margin;//y top coord of plot
		int x1=margin;//x left coord of plot
		int y2=y1+height_plot;//y bottom coord of plot
		int x2=width-margin;//x right coord of plot
		
		//plot value of min and max
		char max_txt[200];	sprintf (max_txt, "%f",currentPlot.max);
		drawText(x1, y1, max_txt, 1.,1.,1.);
		char min_txt[200];	sprintf (min_txt, "%f",currentPlot.min);
		drawText(x1, y2-20, min_txt, 1.,1.,1.);
		
		//and total
		/*float total=0;for(int iv=0;iv<currentPlot.nb_val;iv++)total+=currentPlot.vals[current_time][iv];
		char tot_txt[200];	sprintf (tot_txt, "%f",total);
		drawText(x2-50, y2-20, tot_txt, 1.,1.,1.);*/
		
		//draw limits
		glColor3f(1.,1.,1.);	
		glDisable(GL_LINE_SMOOTH);
		glLineWidth(2.0);
		glBegin(GL_LINES);	glVertex2f(x1,y1);glVertex2f(x2,y1);	glEnd();	//horizontal
		glBegin(GL_LINES);	glVertex2f(x1,y2);glVertex2f(x2,y2);	glEnd();	
		glBegin(GL_LINES);	glVertex2f(x1,y1);glVertex2f(x1,y2);	glEnd();	//vert
		glBegin(GL_LINES);	glVertex2f(x2,y1);glVertex2f(x2,y2);	glEnd();	
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2.0);
				
		//draw 0
		if(currentPlot.max-currentPlot.min!=0)
		{
			float y_zero_rel=-currentPlot.min/(currentPlot.max-currentPlot.min)*height_plot;
			glBegin(GL_LINES);	glVertex2f(x1,y2-y_zero_rel);glVertex2f(x2,y2-y_zero_rel);	glEnd();
			
			//draw values
			for(int iv=0;iv<currentPlot.nb_val;iv++)
			{
				glColor3f(colPlots[iv][0],colPlots[iv][1],colPlots[iv][2]);
				for(int t=0;t<timeLine;t++)
					if(t!=current_time)
				{
					//want to have current time going from left to right
					/*int t0=t;
					int t1=t+1;
					float x_t0=x1+width_plot*t0/timeLine;
					float x_t1=x1+width_plot*t1/timeLine;
					float y_t0_rel=(currentPlot.vals[t][iv]-currentPlot.min)/(currentPlot.max-currentPlot.min)*height_plot;
					float y_t1_rel=(currentPlot.vals[(t+1) % timeLine][iv]-currentPlot.min)/(currentPlot.max-currentPlot.min)*height_plot;
					glBegin(GL_LINES);	glVertex2f(x_t0,y2-y_t0_rel);glVertex2f(x_t1,y2-y_t1_rel);	glEnd();*/
				
					//want to have current time plotted on the right end
					int t0=t-current_time;if(t0>0)t0=t0-timeLine;
					int t1=t+1-current_time;if(t1>0)t1=t1-timeLine;
					float x_t0=x2+width_plot*t0/timeLine;
					float x_t1=x2+width_plot*t1/timeLine;
					float y_t0_rel=(currentPlot.vals[t][iv]-currentPlot.min)/(currentPlot.max-currentPlot.min)*height_plot;
					float y_t1_rel=(currentPlot.vals[(t+1) % timeLine][iv]-currentPlot.min)/(currentPlot.max-currentPlot.min)*height_plot;
					glBegin(GL_LINES);	glVertex2f(x_t0,y2-y_t0_rel);glVertex2f(x_t1,y2-y_t1_rel);	glEnd();
				}
			}
		}
		
		
		
		
	}

	//glBegin(GL_LINES);	glVertex2f(10,20);glVertex2f(30,40);	glEnd();

	glColor3f(1,1,1);	
	if(incTimeWhenDisplay)
		current_time=(current_time+1)%timeLine;
	

}


void PlottingWindow::CreateWindow()
{
#ifdef AMOVERBOSE
	std::cout<<"PlottingWindow::CreateWindow()"<<std::endl;
#endif
	glutInitWindowSize(width, height);	
	window = glutCreateWindow(title.c_str());

}


	