#pragma once

/*#include <cvd/image.h>
#include <cvd/byte.h>
#include <cvd/rgba.h>
*/
#include <objectr3d/MotherWindow.h>
#define MAX_VAL_PER_PLOTS 10

//TextureSet imgTexture;
struct plotStruct
{
	plotStruct(){nb_val=0;timeLine=0;};
	plotStruct(int _nv,int _t){init(_nv,_t);};
	void init(int _nv,int _t)
	{
		desallocate();
		nb_val=_nv;timeLine=_t;
		allocate();	
		min=0;
		max=0;
	}
	void allocate()
	{
		vals=new float*[timeLine];
		for(int i=0;i<timeLine;i++)
		{
			vals[i]=new float[nb_val];
			for(int j=0;j<nb_val;j++)
				vals[i][j]=0;
		}
	}
	void desallocate()
	{
		if(timeLine!=0 && vals!=0)
		{
			for(int i=0;i<timeLine;i++)
				delete[] vals[i];
			delete[] vals;	
			
			timeLine=0;
			nb_val=0;
		}
	};
	~plotStruct(){desallocate();};
	
	void setVal(int _t,int _iv,float _v)
	{
		vals[_t][_iv]=_v;
		if(_v<min)min=_v;
		if(_v>max)max=_v;
	}
	
	//basis storage properties
	int nb_val;
	int timeLine;
	float **vals;//use as vals[timeLine][nb_val]
	
	//dislaying props
	float min,max;
};
	
class PlottingWindow:public MotherWindow
{
 public:
	//creation of drawing function
	PlottingWindow(std::string _title,int _width,int _height,int _nb_plots,int _timeLine,bool _incTimeWhenDisplay=true);
	~PlottingWindow();
	void CreateWindow();
	void prepare_draw();
	
	//setting of each graph
	void setValuesPerPlot(int _p,int _n_val);
	void setVal(int _p,int _ival,float _val);
	
	//increment timeLine manually
	void incrementTimeLine(){current_time=(current_time+1)%timeLine;};
	
	
 protected: 
	int nb_plots;//number of plots
	int timeLine;//number of value plotted per plot
	std::vector<plotStruct> mPlotStructs;
	
	int current_time;
	std::vector<Vector3f> colPlots;
	
	bool incTimeWhenDisplay;
};
