/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libviso2.
Authors: Andreas Geiger

libviso2 is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libviso2 is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libviso2; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

#include "motionEstim.h"




using namespace std;
using namespace Viso;

MotionEstimation::MotionEstimation (float _fu,float _fv,float _cu,float _cv){
      pitch            = 0.0;
      ransac_iters     = 2000;
      //inlier_threshold = 0.00001;
      inlier_threshold = 0.000001;
      motion_threshold = 100.0;
      
      fu=_fu;
      fv=_fv;
      cu=_cu;
      cv=_cv;
}

MotionEstimation::~MotionEstimation () {
	
}


bool MotionEstimation::estimateMotion (std::vector<p_match> p_matched) {

  // get number of matches
  int32_t N = p_matched.size();
  if (N<10)
    return false;
   
  // create calibration matrix
  double K_data[9] = {fu,0,cu,0,fv,cv,0,0,1};
  Matrix K(3,3,K_data);
    
  // normalize feature points and return on errors
  Matrix Tp,Tc;
  vector<p_match> p_matched_normalized = p_matched;
  if (!normalizeFeaturePoints(p_matched_normalized,Tp,Tc))
    return false;

  // initial RANSAC estimate of F
  Matrix E,F;
  inliers.clear();
  for (int32_t k=0;k<ransac_iters;k++) {

    // draw random sample set
    vector<int32_t> active = getRandomSample(N,8);

    // estimate fundamental matrix and get inliers
    fundamentalMatrix(p_matched_normalized,active,F);
    vector<int32_t> inliers_curr = getInlier(p_matched_normalized,F);

    // update model if we are better
    if (inliers_curr.size()>inliers.size())
      inliers = inliers_curr;
  }
  
  // are there enough inliers?
  if (inliers.size()<10)
    return false;
  
  // refine F using all inliers
  fundamentalMatrix(p_matched_normalized,inliers,F); 
  
  // denormalise and extract essential matrix
  F = ~Tc*F*Tp;
  E = ~K*F*K;
  
  // re-enforce rank 2 constraint on essential matrix
  Matrix U,W,V;
  E.svd(U,W,V);
  W.val[2][0] = 0;
  E = U*Matrix::diag(W)*~V;
  
  // compute 3d points X and R|t up to scale
  //Matrix X,R,t;
  Matrix X;
  EtoRt(E,K,p_matched,X,R,t);
  
  if(R.m==0)
	 return false;
  /*std::cout<<"R = "<<std::endl;
  std::cout<<R<<std::endl;
  std::cout<<"t = "<<std::endl;
  std::cout<<t<<std::endl;*/
  
  return true;
}

bool MotionEstimation::estimateMotionAndRemoveOutliers (std::vector<p_match> &p_matched) {

  // get number of matches
  int32_t N = p_matched.size();
  if (N<10)
    return false;
   
  // create calibration matrix
  double K_data[9] = {fu,0,cu,0,fv,cv,0,0,1};
  Matrix K(3,3,K_data);
    
  // normalize feature points and return on errors
  Matrix Tp,Tc;
  vector<p_match> p_matched_normalized = p_matched;
  if (!normalizeFeaturePoints(p_matched_normalized,Tp,Tc))
    return false;

  // initial RANSAC estimate of F
  Matrix E,F;
  inliers.clear();
  for (int32_t k=0;k<ransac_iters;k++) {

    // draw random sample set
    vector<int32_t> active = getRandomSample(N,8);

    // estimate fundamental matrix and get inliers
    fundamentalMatrix(p_matched_normalized,active,F);
    vector<int32_t> inliers_curr = getInlier(p_matched_normalized,F);

    // update model if we are better
    if (inliers_curr.size()>inliers.size())
      inliers = inliers_curr;
  }
  
  // are there enough inliers?
  if (inliers.size()<10)
    return false;
  
  // refine F using all inliers
  fundamentalMatrix(p_matched_normalized,inliers,F); 
  
  // denormalise and extract essential matrix
  F = ~Tc*F*Tp;
  E = ~K*F*K;
  
  // re-enforce rank 2 constraint on essential matrix
  Matrix U,W,V;
  E.svd(U,W,V);
  W.val[2][0] = 0;
  E = U*Matrix::diag(W)*~V;
  
  // compute 3d points X and R|t up to scale
  //Matrix X,R,t;
  Matrix X;
  EtoRtOutlierReject(E,K,p_matched,X,R,t);
  
  /*std::cout<<"R = "<<std::endl;
  std::cout<<R<<std::endl;
  std::cout<<"t = "<<std::endl;
  std::cout<<t<<std::endl;*/
  return true;
}
void MotionEstimation::getEstimatedTransfo(float *rot,float *trans)
{
	//if(R.m!=0)
	{
		int k=0;
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			{
				rot[k]=R.val[j][i];
				k++;
			}
		for(int i=0;i<3;i++)
			trans[i]=t.val[i][0];
	}
		
}


bool MotionEstimation::normalizeFeaturePoints(vector<p_match> &p_matched,Matrix &Tp,Matrix &Tc) {
  
  // shift origins to centroids
  double cpu=0,cpv=0,ccu=0,ccv=0;
  for (vector<p_match>::iterator it = p_matched.begin(); it!=p_matched.end(); it++) {
    cpu += it->u1p;
    cpv += it->v1p;
    ccu += it->u1c;
    ccv += it->v1c;
  }
  cpu /= (double)p_matched.size();
  cpv /= (double)p_matched.size();
  ccu /= (double)p_matched.size();
  ccv /= (double)p_matched.size();
  for (vector<p_match>::iterator it = p_matched.begin(); it!=p_matched.end(); it++) {
    it->u1p -= cpu;
    it->v1p -= cpv;
    it->u1c -= ccu;
    it->v1c -= ccv;
  }
  
  // scale features such that mean distance from origin is sqrt(2)
  double sp=0,sc=0;
  for (vector<p_match>::iterator it = p_matched.begin(); it!=p_matched.end(); it++) {
    sp += sqrt(it->u1p*it->u1p+it->v1p*it->v1p);
    sc += sqrt(it->u1c*it->u1c+it->v1c*it->v1c);
  }
  if (fabs(sp)<1e-10 || fabs(sc)<1e-10)
    return false;
  sp = sqrt(2.0)*(double)p_matched.size()/sp;
  sc = sqrt(2.0)*(double)p_matched.size()/sc;
  for (vector<p_match>::iterator it = p_matched.begin(); it!=p_matched.end(); it++) {
    it->u1p *= sp;
    it->v1p *= sp;
    it->u1c *= sc;
    it->v1c *= sc;
  }
  
  // compute corresponding transformation matrices
  double Tp_data[9] = {sp,0,-sp*cpu,0,sp,-sp*cpv,0,0,1};
  double Tc_data[9] = {sc,0,-sc*ccu,0,sc,-sc*ccv,0,0,1};
  Tp = Matrix(3,3,Tp_data);
  Tc = Matrix(3,3,Tc_data);
  
  // return true on success
  return true;
}

void MotionEstimation::fundamentalMatrix (const vector<p_match> &p_matched,const vector<int32_t> &active,Matrix &F) {
  
  // number of active p_matched
  int32_t N = active.size();
  
  // create constraint matrix A
  Matrix A(N,9);
  for (int32_t i=0; i<N; i++) {
    p_match m = p_matched[active[i]];
    A.val[i][0] = m.u1c*m.u1p;
    A.val[i][1] = m.u1c*m.v1p;
    A.val[i][2] = m.u1c;
    A.val[i][3] = m.v1c*m.u1p;
    A.val[i][4] = m.v1c*m.v1p;
    A.val[i][5] = m.v1c;
    A.val[i][6] = m.u1p;
    A.val[i][7] = m.v1p;
    A.val[i][8] = 1;
  }
   
  // compute singular value decomposition of A
  Matrix U,W,V;
  A.svd(U,W,V);
   
  // extract fundamental matrix from the column of V corresponding to the smallest singular value
  F = Matrix::reshape(V.getMat(0,8,8,8),3,3);
  
  // enforce rank 2
  F.svd(U,W,V);
  W.val[2][0] = 0;
  F = U*Matrix::diag(W)*~V;
}

vector<int32_t> MotionEstimation::getInlier (vector<p_match> &p_matched,Matrix &F) {

  // extract fundamental matrix
  double f00 = F.val[0][0]; double f01 = F.val[0][1]; double f02 = F.val[0][2];
  double f10 = F.val[1][0]; double f11 = F.val[1][1]; double f12 = F.val[1][2];
  double f20 = F.val[2][0]; double f21 = F.val[2][1]; double f22 = F.val[2][2];
  
  // loop variables
  double u1,v1,u2,v2;
  double x2tFx1;
  double Fx1u,Fx1v,Fx1w;
  double Ftx2u,Ftx2v;
  
  // vector with inliers
  vector<int32_t> inliers;
  
  // for all matches do
  for (int32_t i=0; i<(int32_t)p_matched.size(); i++) {

    // extract matches
    u1 = p_matched[i].u1p;
    v1 = p_matched[i].v1p;
    u2 = p_matched[i].u1c;
    v2 = p_matched[i].v1c;
    
    // F*x1
    Fx1u = f00*u1+f01*v1+f02;
    Fx1v = f10*u1+f11*v1+f12;
    Fx1w = f20*u1+f21*v1+f22;
    
    // F'*x2
    Ftx2u = f00*u2+f10*v2+f20;
    Ftx2v = f01*u2+f11*v2+f21;
    
    // x2'*F*x1
    x2tFx1 = u2*Fx1u+v2*Fx1v+Fx1w;
    
    // sampson distance
    double d = x2tFx1*x2tFx1 / (Fx1u*Fx1u+Fx1v*Fx1v+Ftx2u*Ftx2u+Ftx2v*Ftx2v);
    
    // check threshold
    if (fabs(d)<inlier_threshold)
      inliers.push_back(i);
  }

  // return set of all inliers
  return inliers;
}

void MotionEstimation::EtoRt(Matrix &E,Matrix &K,vector<p_match> &p_matched,Matrix &X,Matrix &R,Matrix &t) {

  // hartley matrices
  double W_data[9] = {0,-1,0,+1,0,0,0,0,1};
  double Z_data[9] = {0,+1,0,-1,0,0,0,0,0};
  Matrix W(3,3,W_data);
  Matrix Z(3,3,Z_data); 
  
  // extract T,R1,R2 (8 solutions)
  Matrix U,S,V;
  E.svd(U,S,V);
  Matrix T  = U*Z*~U;
  Matrix Ra = U*W*(~V);
  Matrix Rb = U*(~W)*(~V);
  
  // convert T to t
  t = Matrix(3,1);
  t.val[0][0] = T.val[2][1];
  t.val[1][0] = T.val[0][2];
  t.val[2][0] = T.val[1][0];
  
  // assure determinant to be positive
  if (Ra.det()<0) Ra = -Ra;
  if (Rb.det()<0) Rb = -Rb;
  
  // create vector containing all 4 solutions
  vector<Matrix> R_vec;
  vector<Matrix> t_vec;
  R_vec.push_back(Ra); t_vec.push_back( t);
  R_vec.push_back(Ra); t_vec.push_back(-t);
  R_vec.push_back(Rb); t_vec.push_back( t);
  R_vec.push_back(Rb); t_vec.push_back(-t);
  
  // try all 4 solutions
  Matrix X_curr;
  int32_t max_inliers = 0;
  for (int32_t i=0; i<4; i++) {
    int32_t num_inliers = triangulateChieral(p_matched,K,R_vec[i],t_vec[i],X_curr);
    if (num_inliers>max_inliers) {
      max_inliers = num_inliers;
      X = X_curr;
      R = R_vec[i];
      t = t_vec[i];
    }
  }
}

int32_t MotionEstimation::triangulateChieral (vector<p_match> &p_matched,Matrix &K,Matrix &R,Matrix &t,Matrix &X) {
  
  // init 3d point matrix
  X = Matrix(4,p_matched.size());
  
  // projection matrices
  Matrix P1(3,4);
  Matrix P2(3,4);
  P1.setMat(K,0,0);
  P2.setMat(R,0,0);
  P2.setMat(t,0,3);
  P2 = K*P2;
  
  // triangulation via orthogonal regression
  Matrix J(4,4);
  Matrix U,S,V;
  for (int32_t i=0; i<(int)p_matched.size(); i++) {
    for (int32_t j=0; j<4; j++) {
      J.val[0][j] = P1.val[2][j]*p_matched[i].u1p - P1.val[0][j];
      J.val[1][j] = P1.val[2][j]*p_matched[i].v1p - P1.val[1][j];
      J.val[2][j] = P2.val[2][j]*p_matched[i].u1c - P2.val[0][j];
      J.val[3][j] = P2.val[2][j]*p_matched[i].v1c - P2.val[1][j];
    }
    J.svd(U,S,V);
    X.setMat(V.getMat(0,3,3,3),0,i);
  }
  
  // compute inliers
  Matrix  AX1 = P1*X;
  Matrix  BX1 = P2*X;
  int32_t num = 0;
  for (int32_t i=0; i<X.n; i++)
    if (AX1.val[2][i]*X.val[3][i]>0 && BX1.val[2][i]*X.val[3][i]>0)
      num++;
  
  // return number of inliers
  return num;
}
void MotionEstimation::EtoRtOutlierReject(Matrix &E,Matrix &K,vector<p_match> &p_matched,Matrix &X,Matrix &R,Matrix &t) {

  // hartley matrices
  double W_data[9] = {0,-1,0,+1,0,0,0,0,1};
  double Z_data[9] = {0,+1,0,-1,0,0,0,0,0};
  Matrix W(3,3,W_data);
  Matrix Z(3,3,Z_data); 
  
  // extract T,R1,R2 (8 solutions)
  Matrix U,S,V;
  E.svd(U,S,V);
  Matrix T  = U*Z*~U;
  Matrix Ra = U*W*(~V);
  Matrix Rb = U*(~W)*(~V);
  
  // convert T to t
  t = Matrix(3,1);
  t.val[0][0] = T.val[2][1];
  t.val[1][0] = T.val[0][2];
  t.val[2][0] = T.val[1][0];
  
  // assure determinant to be positive
  if (Ra.det()<0) Ra = -Ra;
  if (Rb.det()<0) Rb = -Rb;
  
  // create vector containing all 4 solutions
  vector<Matrix> R_vec;
  vector<Matrix> t_vec;
  R_vec.push_back(Ra); t_vec.push_back( t);
  R_vec.push_back(Ra); t_vec.push_back(-t);
  R_vec.push_back(Rb); t_vec.push_back( t);
  R_vec.push_back(Rb); t_vec.push_back(-t);
  
  // try all 4 solutions
  Matrix X_curr;
  Matrix isXoulier;
  Matrix isXoulier_curr;
  int32_t max_inliers = 0;
  for (int32_t i=0; i<4; i++) {
    int32_t num_inliers = triangulateChieralOutlierReject(p_matched,K,R_vec[i],t_vec[i],X_curr,isXoulier_curr);
    if (num_inliers>max_inliers) {
      max_inliers = num_inliers;
      isXoulier=isXoulier_curr;
      X = X_curr;
      R = R_vec[i];
      t = t_vec[i];
    }
  }
  
  for (int32_t i=isXoulier.n-1; i>=0; i--)
  {
	  if(isXoulier.val[0][i]==1)
		  p_matched.erase(p_matched.begin()+i);
  }
  
  
}

int32_t MotionEstimation::triangulateChieralOutlierReject (vector<p_match> &p_matched,Matrix &K,Matrix &R,Matrix &t,Matrix &X,Matrix &isXoulier) {
  
  // init 3d point matrix
  X = Matrix(4,p_matched.size());
  isXoulier = Matrix(1,p_matched.size());
  
  // projection matrices
  Matrix P1(3,4);
  Matrix P2(3,4);
  P1.setMat(K,0,0);
  P2.setMat(R,0,0);
  P2.setMat(t,0,3);
  P2 = K*P2;
  
  // triangulation via orthogonal regression
  Matrix J(4,4);
  Matrix U,S,V;
  for (int32_t i=0; i<(int)p_matched.size(); i++) {
    for (int32_t j=0; j<4; j++) {
      J.val[0][j] = P1.val[2][j]*p_matched[i].u1p - P1.val[0][j];
      J.val[1][j] = P1.val[2][j]*p_matched[i].v1p - P1.val[1][j];
      J.val[2][j] = P2.val[2][j]*p_matched[i].u1c - P2.val[0][j];
      J.val[3][j] = P2.val[2][j]*p_matched[i].v1c - P2.val[1][j];
    }
    J.svd(U,S,V);
    X.setMat(V.getMat(0,3,3,3),0,i);
  }
  
  // compute inliers
  Matrix  AX1 = P1*X;
  Matrix  BX1 = P2*X;
  int32_t num = 0;
  for (int32_t i=0; i<X.n; i++)
  {
    if (AX1.val[2][i]*X.val[3][i]>0 && BX1.val[2][i]*X.val[3][i]>0)
    {
	    num++;
	    isXoulier.val[0][i]=0;
    }
    else
	    isXoulier.val[0][i]=1;
  }
  
  // return number of inliers
  return num;
}
vector<int32_t> MotionEstimation::getRandomSample(int32_t N,int32_t num) {

  // init sample and totalset
  vector<int32_t> sample;
  vector<int32_t> totalset;
  
  // create vector containing all indices
  for (int32_t i=0; i<N; i++)
    totalset.push_back(i);

  // add num indices to current sample
  sample.clear();
  for (int32_t i=0; i<num; i++) {
    int32_t j = rand()%totalset.size();
    sample.push_back(totalset[j]);
    totalset.erase(totalset.begin()+j);
  }
  
  // return sample
  return sample;
}