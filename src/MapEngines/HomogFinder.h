//Homography estimation file


#pragma once

#include <objectr3d/AmoDefines.h>
#include <objectr3d/HomogeneousMatrix.h>
#include <objectr3d/Camera.h>
#include <objectr3d/matcher.h>
#include <vector>

//algebraic estimation
Matrix3f HomographyFromMatches(std::vector<p_match> &_hmatches);
//non linear estimation using RANSAC
Matrix3f HomographyFromMatchesRANSAC(std::vector<p_match> &_hmatches,int nb_trial=500);
//non linear estimation using gauss newton (use parameter as initial solution)
void HomographyFromMatchesRobustLM(std::vector<p_match> &_hmatches,Matrix3f &Homography);

float HomogScore(Matrix3f &Homog,p_match &_hmatch);
float HomogScoreTukey(Matrix3f &Homog,p_match &_hmatch);

VectorXf getVectorFromHomog(Matrix3f _homog);
Matrix3f getHomogFromVector(VectorXf _v);
Vector2f homogWarp(Vector2f &x,Matrix3f &Homog);
MatrixXf JacHomogWarp(Vector2f &x,Vector2f &xw,Matrix3f &Homog);




