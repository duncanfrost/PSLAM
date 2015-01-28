#include "GraphFunctions.h"


std::vector<kf_edgeNew> getEdgesInInnerWin(obMap *myMap, std::vector<int> _innerWindowKFs)
{
	std::vector<kf_edgeNew> kf_edges;
	for(int i=0;i<_innerWindowKFs.size();i++)
	{
		KeyFrame *ptKF=myMap->getKF(_innerWindowKFs[i]);
		for(int j=0;j<ptKF->getNbNeigbours();j++)
		{
			int &id_neigbour=ptKF->getPtNeigbour(j)->neighboring_kf;
			if(std::find(_innerWindowKFs.begin(), _innerWindowKFs.end(), id_neigbour)!=_innerWindowKFs.end())
			if(_innerWindowKFs[i]<id_neigbour)//neigboring links are defined in two direction, => this if avoid duplicating edges
			{;
				kf_edgeNew edge;
				edge.kf1=_innerWindowKFs[i];
				edge.kf2=id_neigbour;
				edge.weight=ptKF->getPtNeigbour(j)->edgeScore;
				kf_edges.push_back(edge);
			}
		}
	}
	return kf_edges;  
}


