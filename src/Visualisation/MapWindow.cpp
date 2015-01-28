
#include "MapWindow.h"


namespace MapWindowStuff
{
//pose from which map is displayed
HomogeneousMatrix currentCamPose;
void keyboard_process(unsigned char key, int x, int y);
void mouse_process(int button, int state, int x, int y);
void motion_process(int x, int y);

float z_near_map_win=0.2;

int mouse_old_x, mouse_old_y;
float velocity_rotation=0.05;
float velocity_translation=0.03;
}
  
using namespace MapWindowStuff;


MapWindow::MapWindow(std::string _title,int _width,int _height,Camera *_myCamera)
{
	title=_title;
	width=_width;
	height=_height;
	
	myCamera=_myCamera;
	
	for(int i=0;i<100;i++)
	{
		Vector3f randCol = Vector3f::Random(3);randCol=randCol/randCol.maxCoeff();
		lotsOfRandColors.push_back(randCol);
	}
	lotsOfRandColors[0]=Vector3f(0,1,0.5);
	lotsOfRandColors[1]=Vector3f(0,0.5,1);
	
	closestKF=-1;
	b_showTextures=false;
	b_showFeatureConnections=true;
	b_showLocalFeature=true;
}

void MapWindow::addCamera(HomogeneousMatrix *_camPose,Vector3f _col,float _lineSize)
{
	if(cameraPoses.size()==0 && MapsNew.size()==0)
		currentCamPose=*_camPose;
	cameraPoses.push_back(_camPose);
	colCamera.push_back(_col);
	lineSizeCamera.push_back(_lineSize);
}

void MapWindow::addPointCloud(std::vector<Vector3f> *_map)
{
	pointClouds.push_back(_map);
}

void MapWindow::addMap(obMap *_map)
{
	if(MapsNew.size()==0 && _map->getNbKeyFrames()>0)
		currentCamPose=_map->getKF(0)->getPose();
	
	//generate random color to draw associated keyFrames
	Vector3f randColCamera = Vector3f::Random(3);randColCamera=randColCamera/randColCamera.maxCoeff();
	colMapsNew.push_back(randColCamera);
	
	MapsNew.push_back(_map);	
}

float camera_drawn_size=0.05;
void MapWindow::set_camera_drawn_size(float _f)
{
	camera_drawn_size=_f;
}


void drawCamera(HomogeneousMatrix cameraPose,Camera *myCamera,HomogeneousMatrix &)
{
	HomogeneousMatrix PoseInverse =cameraPose.inverse();
	Vector3f pt_center_cam=currentCamPose*PoseInverse.get_translation();
	Vector2f center_cam=myCamera->Project(pt_center_cam);
		
	Vector2f pt_cam_corners[4];
	pt_cam_corners[0]<< 0 , 0;pt_cam_corners[1]<< myCamera->get_width() , 0;
	pt_cam_corners[3]<< 0 , myCamera->get_height();pt_cam_corners[2]<< myCamera->get_width() , myCamera->get_height();
	
	float depth_corner=camera_drawn_size;
	if(pt_center_cam[2]>z_near_map_win)
	for(int i=0;i<4;i++)
	{
		//draw center to corner
		Vector3f pt_corner=currentCamPose*PoseInverse*(depth_corner*myCamera->UnProjectZ1(pt_cam_corners[i]));
		Vector3f pt_corner2=currentCamPose*PoseInverse*(depth_corner*myCamera->UnProjectZ1(pt_cam_corners[(i+1)%4]));
		if(pt_corner[2]>z_near_map_win && pt_corner2[2]>z_near_map_win)
		{
			Vector2f corner=myCamera->Project(pt_corner);
			glBegin(GL_LINES);	glVertex2f(center_cam[0],center_cam[1]);glVertex2f(corner[0],corner[1]);	glEnd();
			
			//draw corner to next corner
			Vector2f corner2=myCamera->Project(pt_corner2);
			glBegin(GL_LINES);	glVertex2f(corner2[0],corner2[1]);glVertex2f(corner[0],corner[1]);	glEnd();
		}
	}	
}

void MapWindow::prepare_draw()
{
	//glClearColor(1.0f, 1.0f, 1.0f, 0.0f);//set background as white
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	set2DGLProjection();
	
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glPointSize(3.0);

	GLint factor = 1;    // Stippling factor
	//GLushort pattern = 0x5555;//alternate
	GLushort pattern = 0x0C0F;//normal
	{//new system with entire map given as parameter

		//draw neigbouring links
		glLineWidth(2);
		float maxScoreEdge=0;
		for(int m=0;m<MapsNew.size();m++)
		{
			std::vector<int> kfToShow;
			if(closestKF==-1)
				for(int c=0;c<MapsNew[m]->getNbKeyFrames();c++)kfToShow.push_back(c);
			else
			{
				std::vector<int> closestVect;closestVect.push_back(closestKF);
				kfToShow=MapsNew[m]->getConnectedKeyframes(closestKF);
			}
			std::sort(kfToShow.begin(),kfToShow.end());
		  
		      
			//get max score edge
			for(int c=0;c<kfToShow.size();c++)
			{
				glColor3f(1.-colMapsNew[m][0],0.5*colMapsNew[m][1],0.5*colMapsNew[m][2]);
				//center cam current		
				KeyFrame *tKF=MapsNew[m]->getKF(kfToShow[c]);
				HomogeneousMatrix PoseInverse=tKF->getPose().inverse();
				
				Vector3f pt_center_cam=currentCamPose*PoseInverse.get_translation();
				Vector2f center_cam=myCamera->Project(pt_center_cam);
				if(pt_center_cam[2]>z_near_map_win)
				for(int i=0;i<tKF->getNbNeigbours();i++)
				{
					float scoreEdge=tKF->getPtNeigbour(i)->edgeScore;
					if(scoreEdge>maxScoreEdge)
					  maxScoreEdge=scoreEdge;
				}
			}

			for(int c=0;c<kfToShow.size();c++)
			{
				glColor3f(1.-colMapsNew[m][0],0.5*colMapsNew[m][1],0.5*colMapsNew[m][2]);
				//center cam current		
				KeyFrame *tKF=MapsNew[m]->getKF(kfToShow[c]);
				HomogeneousMatrix PoseInverse=tKF->getPose().inverse();
				
				Vector3f pt_center_cam=currentCamPose*PoseInverse.get_translation();
				Vector2f center_cam=myCamera->Project(pt_center_cam);
				if(pt_center_cam[2]>z_near_map_win)
				for(int i=0;i<tKF->getNbNeigbours();i++)
				{
					int neigb=tKF->getPtNeigbour(i)->neighboring_kf;
					float scoreEdge=tKF->getPtNeigbour(i)->edgeScore;
					float lineWidth=4*scoreEdge/maxScoreEdge;
					if(lineWidth<0.5)lineWidth=0.5;
					glLineWidth(lineWidth);
					HomogeneousMatrix PoseInverse2 =MapsNew[m]->getKF(neigb)->getPose().inverse();
					Vector3f pt_center_cam2=currentCamPose*PoseInverse2.get_translation();
					if(pt_center_cam2[2]>z_near_map_win)
					{
						Vector2f center_cam2=myCamera->Project(pt_center_cam2);
						glBegin(GL_LINES);	glVertex2f(center_cam[0],center_cam[1]);glVertex2f(center_cam2[0],center_cam2[1]);	glEnd();
					}
				}
			}
			//draw temp point miniBA
			/*glColor3f(1.0,0.0,0.0);
			glPointSize(5.0);
			for(int c=0;c<kfToShow.size();c++)
			{
				KeyFrame *tKF=MapsNew[m]->getKF(kfToShow[c]);			
				for(int i=0;i<tKF->nb_feat_ref;i++)
				{
					if(tKF->invDepthFeat[i]!=-1)
					{
						float idepth=tKF->invDepthFeat[i];
						if(idepth<1e-6)idepth=1e-6;
						
						Vector3f mapPointsc=tKF->getPose().inverse()*(toHomogeneous(tKF->MeanFeatPos[i])/idepth);
						Vector3f pt_c2;pt_c2=currentCamPose* mapPointsc;
						if(pt_c2[2]>0)
						{
							Vector2f pt_pc2=myCamera->Project(pt_c2);
							glBegin(GL_POINTS);	
							glVertex2f(pt_pc2[0],pt_pc2[1]);	
							glEnd();
						}
					}
				}
			}*/
		
			//draw points and features 
			glEnable(GL_LINE_STIPPLE);
			glLineWidth(1);
			glLineStipple(factor,pattern);	//factor= 			

			//cameras
			for(int c=0;c<kfToShow.size();c++)
			{
				KeyFrame *tKF=MapsNew[m]->getKF(kfToShow[c]);
				
				//get max score of points for this KF
				float score_max=0;
				for(int i=0;i<tKF->getNbMapPoint();i++)
				{
					MapPoint *point=tKF->getPtMapPoint(i);
					if(point->isUsed())					
					{
						float score_p=point->getWeight();
						if(score_p>score_max)score_max=score_p;

					}
				}
				
				
				//draw map points
				//std::cout<<"Kf "<<tKF->getId()<<"nb pt = "<<tKF->getNbMapPoint()<<std::endl;
				for(int i=0;i<tKF->getNbMapPoint();i++)
				{
					MapPoint *point=tKF->getPtMapPoint(i);
#ifndef SHOW_BAD_POINTS
					if(point->isUsed())					
#endif
					{
						float score_p=point->getWeight();

						//draw point
#ifdef  SAVE_POINT_COLOR						
						if(b_showTextures)
						{
							//float grayLvlGl=(float)point->getGrayVal()/255.;
							
							float col[3];
							for(int ch=0;ch<3;ch++)
								col[ch] =(float)point->getCol(ch)/255.;
							glColor3f(col[0],col[1],col[2]);							
						}
						else
						  	glColor3f(1.0,0.2,0.2);//in purple
#else
						glColor3f(1.0,0.2,0.2);//in purple
#endif						
						if(!point->isUsed())
						{
							glColor3f(1.0,0.0,0.0);
							glPointSize(5.0);
						}
						else
							glPointSize(7.0*score_p/score_max);
						//glPointSize(score_p*5.0*10);
							

						Vector3f pt_c1;pt_c1=currentCamPose* point->getPosition();
						if(pt_c1[2]>0)
						{
							Vector2f pt_pc1=myCamera->Project(pt_c1);
							glBegin(GL_POINTS);	glVertex2f(pt_pc1[0],pt_pc1[1]);	glEnd();
							
							//draw views
							if(b_showFeatureConnections)
							{
								glColor3f(1.0,0.0,0.0);
								//std::cout<<"Kf "<<tKF->getId()<<" pt = "<<point->getId()<<" nb views = "<<point->nbViews()<<std::endl;
								for(int v=0;v<point->nbViews();v++)
								{
									int kfView=point->getView(v);
									int i1pView=point->getI1p(v);
									//get corresponding local feature:
									int id_feat_v=MapsNew[m]->getKF(kfView)->indexCandidateFeatureFromVisoId(i1pView);
									if(id_feat_v!=-1)
									{
										uptoscaleFeature &feat_v=*MapsNew[m]->getKF(kfView)->getPtLocalBestFeatures(id_feat_v);
										HomogeneousMatrix invHv=MapsNew[m]->getKF(kfView)->getPose().inverse();
										Vector3f pt_c2;pt_c2=currentCamPose* invHv*(toHomogeneous(feat_v.posRef)*feat_v.depthInRef);
										if(pt_c2[2]>0)
										{
											Vector2f pt_pc2=myCamera->Project(pt_c2);
											glBegin(GL_LINES);	
											glVertex2f(pt_pc2[0],pt_pc2[1]);	
											glVertex2f(pt_pc1[0],pt_pc1[1]);	
											glEnd();
										}
									}
								}
							}
							
						}
					}
				}
				
				//draw kf features
				if(b_showLocalFeature)
				{
					HomogeneousMatrix invH=tKF->getPose().inverse();
					
					//get rec Angle max
					float recAngl_max=0;
					for(int i=0;i<tKF->getNbLocalBestFeatures();i++)
					{
						uptoscaleFeature &tFeat=*tKF->getPtLocalBestFeatures(i);
						
						float recAngl=tFeat.recAngle;
						if(recAngl>recAngl_max)recAngl_max=recAngl;
					}		
					
					//glPointSize(3.0);
					for(int i=0;i<tKF->getNbLocalBestFeatures();i++)
					{
						uptoscaleFeature &tFeat=*tKF->getPtLocalBestFeatures(i);
						
						float score_f=tFeat.recAngle;
						//glPointSize(score_f*5.0/1.5);
						
						Vector3f pt_c1;pt_c1=currentCamPose* invH*(toHomogeneous(tFeat.posRef)*tFeat.depthInRef);
						if(pt_c1[2]>0)
						{

							glPointSize(3.0*score_f/recAngl_max);

							
	#if  defined (SAVE_POINT_COLOR)						
							if(b_showTextures)
							{
								//float grayLvlGl=(float)tFeat.grayVal/255.;
								//glColor3f(grayLvlGl,grayLvlGl,grayLvlGl);//in purple
								
								float col[3];
								for(int ch=0;ch<3;ch++)
									col[ch] =(float)tFeat.col[ch]/255.;
								glColor3f(col[0],col[1],col[2]);//in purple
							}
							else
							{
								Vector3f col_cam=lotsOfRandColors[c % lotsOfRandColors.size()];
								glColor3f(col_cam[0],col_cam[1],col_cam[2]);
							}
	#else
							Vector3f col_cam=lotsOfRandColors[c % lotsOfRandColors.size()];
							glColor3f(col_cam[0],col_cam[1],col_cam[2]);
	#endif	
							
							Vector2f pt_pc1=myCamera->Project(pt_c1);
							glBegin(GL_POINTS);	glVertex2f(pt_pc1[0],pt_pc1[1]);	glEnd();
							
							//draw link with point (use one color per point)
							/*if(tFeat.matched)
							{
								//int id_col=100*tFeat.ptKForigin->getId()+tFeat.idPoint;
								//Vector3f col_p=lotsOfRandColors[id_col%lotsOfRandColors.size()];
								//Vector3f col_p=lotsOfRandColors[tFeat.ptKForigin->getId()%lotsOfRandColors.size()];
								//glColor3f(col_p[0],col_p[1],col_p[2]);//in purple
								int id_ref=tFeat.ptKForigin->getId();
								glColor3f(col_cam[id_ref],col_cam[id_ref],col_cam[id_ref]);//in purple
								
								MapPoint *point=tFeat.ptKForigin->getPtMapPoint(tFeat.idPoint);
								Vector3f pt_p;pt_p=currentCamPose* point->getPosition();
								if(pt_c1[2]>0)
								{
									Vector2f pt_=myCamera->Project(pt_p);
									glBegin(GL_LINES);	
									glVertex2f(pt_pc1[0],pt_pc1[1]);	
									glVertex2f(pt_[0],pt_[1]);	
									glEnd();
								}
								
							}*/
						}
						
						
					}
				}
				
				

			}
			
		
			glLineStipple(0,pattern);	
			glDisable(GL_LINE_STIPPLE);
			
			//draw cameras 
			glDisable(GL_LINE_STIPPLE);
			glLineWidth(2);

			//cameras
			for(int c=0;c<kfToShow.size();c++)
			{
				KeyFrame *tKF=MapsNew[m]->getKF(kfToShow[c]);
				
				//cehck if closestKF or if in active window
				if(tKF->getId()==closestKF)//closest
					glColor3f(1.,0.6,0);
				else if(std::find(activeKF.begin(), activeKF.end(), tKF->getId())!=activeKF.end())
					glColor3f(1.,1,0);
				else
					glColor3f(colMapsNew[m][0],colMapsNew[m][1],colMapsNew[m][2]);
				drawCamera(tKF->getPose(),myCamera,currentCamPose);
			}

		
			//draw cameras best local pair
			
	#ifdef SHOW_BEST_KF_PAIR		
			glLineWidth(1);

			//cameras
			for(int c=0;c<kfToShow.size();c++)
			{
				KeyFrame *tKF=MapsNew[m]->getKF(kfToShow[c]);
				
				glEnable(GL_LINE_STIPPLE);
				glColor3f(1-colMapsNew[m][0],colMapsNew[m][1],1-colMapsNew[m][2]);
				HomogeneousMatrix camPair=tKF->getBestRelPose()*tKF->getPose();
				drawCamera(camPair,myCamera,currentCamPose);
				
				//draw connection with camera
				//glDisable(GL_LINE_STIPPLE);
				HomogeneousMatrix CamInv =tKF->getPose().inverse();
				HomogeneousMatrix CamPairInv =camPair.inverse();
				Vector3f pt_center_cam=currentCamPose*CamInv.get_translation();
				Vector3f pt_center_camInv=currentCamPose*CamPairInv.get_translation();
				if(pt_center_camInv[2]>z_near_map_win && pt_center_cam[2]>z_near_map_win)
				{
					Vector2f center_cam=myCamera->Project(pt_center_cam);
					Vector2f center_cam2=myCamera->Project(pt_center_camInv);
					glBegin(GL_LINES);	glVertex2f(center_cam[0],center_cam[1]);glVertex2f(center_cam2[0],center_cam2[1]);	glEnd();
				}
				
			}
		}
		glDisable(GL_LINE_STIPPLE);
#endif
		

	}

	{//for old system when camera and point where added independently
		//Draw maps
		glColor3f(0.5,1,0.7);	
		for(int m=0;m<pointClouds.size();m++)
		{
			for(int i=0;i<(*pointClouds[m]).size();i++)
			{
				Vector3f pt_c1;pt_c1=currentCamPose* (*pointClouds[m])[i];			
				Vector2f pt_pc1=myCamera->Project(pt_c1);
				glBegin(GL_POINTS);	glVertex2f(pt_pc1[0],pt_pc1[1]);	glEnd();
			}
		}
		
		for(int m=0;m<pointInvClouds.size();m++)
		{
			std::vector<PointInvDepth> &refPtCloud=*pointInvClouds[m];
			for(int i=0;i<refPtCloud.size();i++)
			{
				float idepth;
				if(refPtCloud[i].invDepth<1e-6)idepth=1e-6;
				else idepth=refPtCloud[i].invDepth;
				//draw incertitude
				//glColor3f(1.0,0.5,0.7);
				glColor3f(lotsOfRandColors[m][0],lotsOfRandColors[m][1],lotsOfRandColors[m][2]);
				
				Vector2f pt_pcc1;Vector2f pt_pcc2;
				Vector3f PointWorld1,PointWorld2;
				{
				  Vector3f mapPointsCam1=toHomogeneous(refPtCloud[i].meterCoord)/(refPtCloud[i].invDepth+refPtCloud[i].invDepthCovar);
				  PointWorld1=cameraPoses[refPtCloud[i].srcKf]->inverse()*mapPointsCam1;
				  
				  Vector3f pt_c1;pt_c1=currentCamPose* PointWorld1;			
				  pt_pcc1=myCamera->Project(pt_c1);
				}
				{
				  Vector3f mapPointsCam1=toHomogeneous(refPtCloud[i].meterCoord)/(refPtCloud[i].invDepth-refPtCloud[i].invDepthCovar);
				  PointWorld2=cameraPoses[refPtCloud[i].srcKf]->inverse()*mapPointsCam1;
				  
				  Vector3f pt_c1;pt_c1=currentCamPose* PointWorld2;			
				  pt_pcc2=myCamera->Project(pt_c1);
				}
				if(PointWorld1[2]>0 && PointWorld2[2]>0)
				{
					glBegin(GL_LINES);	glVertex2f(pt_pcc1[0],pt_pcc1[1]);glVertex2f(pt_pcc2[0],pt_pcc2[1]);	glEnd();
				}
				
				//glColor3f(0.5,1,0.7);	
				//glColor3f(1.,0.,0.);	
			        Vector3f mapPointsCam1;
				if(refPtCloud[i].invDepth==0)
					mapPointsCam1=toHomogeneous(refPtCloud[i].meterCoord)*1e10;
				else
					mapPointsCam1=toHomogeneous(refPtCloud[i].meterCoord)/refPtCloud[i].invDepth;
				Vector3f PointWorld=cameraPoses[refPtCloud[i].srcKf]->inverse()*mapPointsCam1;
				
				Vector3f pt_c1;pt_c1=currentCamPose* PointWorld;			
				Vector2f pt_pc1=myCamera->Project(pt_c1);
				if(PointWorld[2]>0)
				{
					glBegin(GL_POINTS);	glVertex2f(pt_pc1[0],pt_pc1[1]);	glEnd();
				}
				
			}

		}
		
		//draw error between map1 and map2 if there is more than one map
		if(pointClouds.size()>1)
		{
			glColor3f(1.,0.,0.);	
			for(int i=0;i<(*pointClouds[0]).size();i++)
			{
				Vector3f pt_c1;pt_c1=currentCamPose* (*pointClouds[0])[i];			
				Vector2f pt_pc1=myCamera->Project(pt_c1);
				Vector3f pt_c2;pt_c2=currentCamPose* (*pointClouds[1])[i];			
				Vector2f pt_pc2=myCamera->Project(pt_c2);
				glBegin(GL_LINES);	glVertex2f(pt_pc1[0],pt_pc1[1]);glVertex2f(pt_pc2[0],pt_pc2[1]);	glEnd();
			}
			
		}
		
		//draw cameras
		for(int c=0;c<cameraPoses.size();c++)
		{
			glLineWidth(lineSizeCamera[c]);
			glColor3f(colCamera[c][0],colCamera[c][1],colCamera[c][2]);
			drawCamera(*cameraPoses[c],myCamera,currentCamPose);
		}
		glLineWidth(1.);
	}	
	
	glColor3f(1,1,1);	
	//glutSwapBuffers();	
}


void MapWindow::CreateWindow()
{
#ifdef AMOVERBOSE
	std::cout<<"MapWindow::CreateWindow()"<<std::endl;
#endif
	glutInitWindowSize(width, height);	
	window = glutCreateWindow(title.c_str());
}

void MapWindow::setEvents()
{
	//get id in static array
	glutKeyboardFunc(keyboard_process);
	glutMouseFunc(mouse_process);
	glutMotionFunc(motion_process);         
}



void MapWindowStuff::keyboard_process(unsigned char key, int x, int y)
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

void MapWindow::set_velocity_translation(float _f){velocity_translation=_f;}

void MapWindowStuff::mouse_process(int button, int state, int x, int y)
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

void MapWindowStuff::motion_process(int x, int y)
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
void MapWindow::moveCamera(HomogeneousMatrix _cHc2)
{
	currentCamPose=_cHc2*currentCamPose;
}
void MapWindow::setCameraPose(HomogeneousMatrix _cHc2)
{
	currentCamPose=_cHc2;	
}
HomogeneousMatrix MapWindow::getCameraPose()
{
	return currentCamPose;
}
	