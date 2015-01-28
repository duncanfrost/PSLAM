#include "MapPoint.h"

void MapPoint::removeView(int idkf)
{
	if(idKFmeas.size()!=0)
	{
		std::vector<int>::iterator it=std::find(idKFmeas.begin(), idKFmeas.end(), idkf);
		if(it!=idKFmeas.end())  
		{
			int index=it-idKFmeas.begin();
			//remove correponding i1p
			i1p.erase(i1p.begin()+index);
			idKFmeas.erase(it);
		}
		else
			std::cout<<"Point["<<id<<"] : try to remove non existing view "<<std::endl;
	}
	else
		std::cout<<"Point["<<id<<"] : try to remove non existing view "<<std::endl;
		
}
void MapPoint::removeViewId(int idv)
{
	if(idKFmeas.size()>idv)
	{
		i1p.erase(i1p.begin()+idv);
		idKFmeas.erase(idKFmeas.begin()+idv);		  
	}
}

void MapPoint::saveToStream(std::ofstream &fout)
{
	char isused=used;
	fout.write((const char*)&isused,sizeof(char));
	
	for(int i=0;i<3;i++)
	      fout.write((const char*)&position[i],sizeof(float));

	fout.write((const char*)&id,sizeof(int));
	fout.write((const char*)&outlierCount,sizeof(short));
	fout.write((const char*)&inlierCount,sizeof(short));

	fout.write((const char*)&weight,sizeof(float));

	int nb_view=idKFmeas.size();
	fout.write((const char*)&nb_view,sizeof(int));
	for(int v=0;v<nb_view;v++)
	{
		fout.write((const char*)&idKFmeas[v],sizeof(int));
		fout.write((const char*)&i1p[v],sizeof(int));
	}
#ifdef SAVE_POINT_COLOR
	//fout.write((const char*)&grayVal,sizeof(char));
	for(int i=0;i<3;i++)
		fout.write((const char*)&col[i],sizeof(char));
	#endif
}

void MapPoint::loadFromStream(std::ifstream &fout)
{
	char isused;	fout.read((char*)&isused,sizeof(char));
	used=(isused==1);

	for(int i=0;i<3;i++)
	      fout.read((char*)&position[i],sizeof(float));
	
	fout.read((char*)&id,sizeof(int));
	fout.read((char*)&outlierCount,sizeof(short));
	fout.read((char*)&inlierCount,sizeof(short));
	fout.read((char*)&weight,sizeof(float));
	
	int nb_view;
	fout.read((char*)&nb_view,sizeof(int));
	idKFmeas.resize(nb_view);
	i1p.resize(nb_view);
	for(int v=0;v<nb_view;v++)
	{
		fout.read((char*)&idKFmeas[v],sizeof(int));
		fout.read((char*)&i1p[v],sizeof(int));
	}
#ifdef SAVE_POINT_COLOR
	//fout.read((char*)&grayVal,sizeof(char));
	for(int i=0;i<3;i++)
		fout.read((char*)&col[i],sizeof(char));
#endif


}