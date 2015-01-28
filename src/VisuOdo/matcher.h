/*
Copyright 2012. All rights reserved.
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

#ifndef __MATCHER3_H__
#define __MATCHER3_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <emmintrin.h>
#include <algorithm>
#include <vector>

#include <Eigen/Core>

#include "matrix.h"

  // structure for storing matches
  struct p_match {
    float   u1p,v1p; // u,v-coordinates in previous left  image
    int32_t i1p;     // feature index (for tracking)
    float   u1c,v1c; // u,v-coordinates in current  left  image
    int32_t i1c;     // feature index (for tracking)
    p_match(){}
    p_match(float u1p,float v1p,int32_t i1p,
            float u1c,float v1c,int32_t i1c):
            u1p(u1p),v1p(v1p),i1p(i1p),
            u1c(u1c),v1c(v1c),i1c(i1c) {}
    p_match reverse(){return p_match(u1c,v1c,i1c,u1p,v1p,i1p);} ;   

  };

  
  struct p_feat {
    float   u1c,v1c; // u,v-coordinates in current  left  image
    int32_t i1c;     // feature index (for tracking)
    
    int32_t c;   // class
    p_feat(float u1p,float v1p,int32_t i1p,int32_t c):
    u1c(u1p),v1c(v1p),i1c(i1p),c(c){}    ;
  };
  
class Matcher {

public:
  // parameter settings
  struct parameters {
  
    int32_t nms_n;                  // non-max-suppression: min. distance between maxima (in pixels)
    int32_t nms_tau;                // non-max-suppression: interest point peakiness threshold
    int32_t match_binsize;          // matching bin width/height (affects efficiency only)
    int32_t match_radius;           // matching radius (du/dv in pixels)
    int32_t match_disp_tolerance;   // dv tolerance for stereo matches (in pixels)
    int32_t outlier_disp_tolerance; // outlier removal: disparity tolerance (in pixels)
    int32_t outlier_flow_tolerance; // outlier removal: flow tolerance (in pixels)
    int32_t multi_stage;            // 0=disabled,1=multistage matching (denser and faster)
    int32_t half_resolution;        // 0=disabled,1=match at half resolution, refine at full resolution
    int32_t refinement;             // refinement (0=none,1=pixel,2=subpixel)
    double  f,cu,cv,base;           // calibration (only for match prediction)
    bool refIsPrec;		    //which image we match is it first to current or previous to current
    
    // default settings
    parameters () {
      nms_n                  = 3;
      nms_tau                = 50;
      match_binsize          = 50;
      match_radius           = 50;
      match_disp_tolerance   = 2;
      outlier_disp_tolerance = 5;
      outlier_flow_tolerance = 5;
      multi_stage            = 1;
      half_resolution        = 1;
      refinement             = 2;
      refIsPrec             = 0;
    }
  };

  // constructor (with default parameters)
  Matcher(parameters param=parameters());

  // deconstructor
  ~Matcher();
  


  // computes features from left/right images and pushes them back to a ringbuffer,
  // which interally stores the features of the current and previous image pair
  // use this function for stereo or quad matching
  // input: I1,I2 .......... pointers to left and right image (row-aligned), range [0..255]
  //        dims[0,1] ...... image width and height (both images must be rectified and of same size)
  //        dims[2] ........ bytes per line (often equals width)
  //        replace ........ if this flag is set, the current image is overwritten with
  //                         the input images, otherwise the current image is first copied
  //                         to the previous image (ring buffer functionality, descriptors need
  //                         to be computed only once)    
  //void pushBack (uint8_t *I1,uint8_t* I2,int32_t* dims,const bool replace);
  
  // computes features from a single image and pushes it back to a ringbuffer,
  // which interally stores the features of the current and previous image pair
  // use this function for flow computation
  // parameter description see above
  std::vector<p_feat> getFeatures(uint8_t *I1,int32_t* dims);
  std::vector<p_feat> getFeaturesp();
  std::vector<p_feat> getFeaturespClass(int c);
  std::vector<p_feat> getFeaturesc();
  void pushBack (uint8_t *I1,int32_t* dims,const bool init=false);//init used if refIsPrec=0 to change first image
  void pushBackRef (uint8_t *I1,int32_t* dims,std::vector<p_feat> &featureKlt);//init used if refIsPrec=0 to change first image
  void pushBack (uint8_t *I1,int32_t* dims,std::vector<p_feat> &featureI1);//init used if refIsPrec=0 to change first image
  void toWarpedId(std::vector<p_match> &_matches){for(int i=0;i<_matches.size();i++)_matches[i].i1c=LUTtoIDwarp[_matches[i].i1c];};
  void computeDescriptorViso(uint8_t *I1,int32_t* dims,std::vector<Eigen::Vector2f> &maxCorners,int32_t *&_VisoDescriptors,int32_t &_nVisoDescriptor);
  void setTwoFrames(int32_t* dims,uint8_t *I1,int32_t *_VisoDescriptors1,int32_t _nVisoDescriptor1,uint8_t *I2,int32_t *_VisoDescriptors2,int32_t _nVisoDescriptor2);
  
  //get second pass features
  std::vector<Eigen::Vector2f> GetListFeatureCoord(int pass=2);

  // match features currently stored in ring buffer (current and previous frame)
  // input: method ... 0 = flow, 1 = stereo, 2 = quad matching
  //        Tr_delta: uses motion from previous frame to better search for
  //                  matches, if specified
  void matchFeatures(Viso::Matrix *Tr_delta = 0);

  // feature bucketing: keeps only max_features per bucket, where the domain
  // is split into buckets of size (bucket_width,bucket_height)
  void bucketFeatures(int32_t max_features,float bucket_width,float bucket_height);

  //return number of feature of pass 2 in previous image
  int getNbFeatureImg1() { return n1p2; }
  int getNbFeaturesRef() { return n1p2; }
 void printFeatPosFromI1p(int i){int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
				      if(i<n1p2)
					std::cout<<"feat ["<<i<<"]"<<*(m1p2+step_size*i+0)<<" "<<*(m1p2+step_size*i+1)<<std::endl;else std::cout<<"i1p not exist"<<std::endl;
    
  };
 void printFeatPosFromI1c(int i){int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
				      if(i<n1c2)
					std::cout<<"feat ["<<i<<"]"<<*(m1c2+step_size*i+0)<<" "<<*(m1c2+step_size*i+1)<<std::endl;else std::cout<<"i1c not exist"<<std::endl;
    
  };
  
  Eigen::Vector2f getCoordMeasurePix(int i) 
	    { int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);  
	      return Eigen::Vector2f(*(m1p2+step_size*i+0),*(m1p2+step_size*i+1)); }
	      
  std::vector< Eigen::Vector2f > getFeaturesImgp() 
	    { 
	      std::vector< Eigen::Vector2f > res;
	      int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);  
	      for(int i=0;i<n1p2;i++)
		res.push_back(Eigen::Vector2f(*(m1p2+step_size*i+0),*(m1p2+step_size*i+1))); 
	      return res;
	    }
	      
   // return vector with matched feature points and indices
  std::vector<p_match> getMatches() { return p_matched_2; }
  std::vector<p_match> getMatchesFromClass(int c)
  {
	  std::vector<p_match> matchesClassC;
	  for(int i=0;i<p_matched_2.size();i++)
	  {
		  int32_t step_size = sizeof(Matcher::maximum)/sizeof(int32_t);
		  if(*(m1p2+step_size*p_matched_2[i].i1p+3)==c)
			  matchesClassC.push_back(p_matched_2[i]);
	  }
	  return matchesClassC;
	  
  }
  ;
  void updateMatchesWithOffsetMatches(std::vector<p_match> &updatedMatch){p_matched_2=updatedMatch;};
  void updateMotionPriorWithNewMatches(){computePriorStatistics(p_matched_2);};
  
  
  // u/v ranges for matching stage 0-3
  struct range {
    float u_min[4];
    float u_max[4];
    float v_min[4];
    float v_max[4];
  };
  
  //update motion prior by providing ranges directly
  int getBinsize(){return param.match_binsize;};
  void setRanges(std::vector<range> &_ranges){ranges=_ranges;};

  // given a vector of inliers computes gain factor between the current and
  // the previous frame. this function is useful if you want to reconstruct 3d
  // and you want to cancel the change of (unknown) camera gain.
  float getGain (std::vector<int32_t> inliers);

private:

  // structure for storing interest points
  struct maximum {
    int32_t u;   // u-coordinate
    int32_t v;   // v-coordinate
    int32_t val; // value
    int32_t c;   // class
    int32_t d1,d2,d3,d4,d5,d6,d7,d8; // descriptor //actually it is 32 gray intensity around point
    maximum() {}
    maximum(int32_t u,int32_t v,int32_t val,int32_t c):u(u),v(v),val(val),c(c) {}
  };

  
  struct delta {
    float val[8];
    delta () {}
    delta (float v) {
      for (int32_t i=0; i<8; i++)
        val[i] = v;
    }
  };
  
  // computes the address offset for coordinates u,v of an image of given width
  inline int32_t getAddressOffsetImage (const int32_t& u,const int32_t& v,const int32_t& width) {
    return v*width+u;
  }

  // Alexander Neubeck and Luc Van Gool: Efficient Non-Maximum Suppression, ICPR'06, algorithm 4
  void nonMaximumSuppression (int16_t* I_f1,int16_t* I_f2,const int32_t* dims,std::vector<Matcher::maximum> &maxima,int32_t nms_n);
  void extractFastFeatures(uint8_t *I,const int32_t* dims,const int32_t* dims_dest,std::vector<Matcher::maximum> &maxima);

  // descriptor functions
  inline uint8_t saturate(int16_t in);
  void filterImageAll (uint8_t* I,uint8_t* I_du,uint8_t* I_dv,int16_t* I_f1,int16_t* I_f2,const int* dims);
  void filterImageSobel (uint8_t* I,uint8_t* I_du,uint8_t* I_dv,const int* dims);
  inline void computeDescriptor (const uint8_t* I_du,const uint8_t* I_dv,const int32_t &bpl,const int32_t &u,const int32_t &v,uint8_t *desc_addr);
  inline void computeSmallDescriptor (const uint8_t* I_du,const uint8_t* I_dv,const int32_t &bpl,const int32_t &u,const int32_t &v,uint8_t *desc_addr);
  void computeDescriptors (uint8_t* I_du,uint8_t* I_dv,const int32_t bpl,std::vector<Matcher::maximum> &maxima);
  
  void getHalfResolutionDimensions(const int32_t *dims,int32_t *dims_half);
  uint8_t* createHalfResolutionImage(uint8_t *I,const int32_t* dims);

  // compute sparse set of features from image
  // inputs:  I ........ image
  //          dims ..... image dimensions [width,height]
  //          n ........ non-max neighborhood
  //          tau ...... non-max threshold
  // outputs: max ...... vector with maxima [u,v,value,class,descriptor (128 bits)]
  //          I_du ..... gradient in horizontal direction
  //          I_dv ..... gradient in vertical direction
  // WARNING: max,I_du,I_dv has to be freed by yourself!
  void computeFeatures (uint8_t *I,const int32_t* dims,int32_t* &max1,int32_t &num1,int32_t* &max2,int32_t &num2,uint8_t* &I_du,uint8_t* &I_dv,uint8_t* &I_du_full,uint8_t* &I_dv_full);
  void computeFeaturesPlusKlt (uint8_t *I,const int32_t* dims,int32_t* &max1,int32_t &num1,int32_t* &max2,int32_t &num2,uint8_t* &I_du,uint8_t* &I_dv,uint8_t* &I_du_full,uint8_t* &I_dv_full,std::vector<p_feat> &featureKlt);
  void computeFeaturesForceFeats (uint8_t *I,const int32_t* dims,int32_t* &max1,int32_t &num1,int32_t* &max2,int32_t &num2,uint8_t* &I_du,uint8_t* &I_dv,uint8_t* &I_du_full,uint8_t* &I_dv_full,std::vector<p_feat> &featureI);

  // matching functions
  void computePriorStatistics (std::vector<p_match> &p_matched);
  void createIndexVector (int32_t* m,int32_t n,std::vector<int32_t> *k,const int32_t &u_bin_num,const int32_t &v_bin_num);
  inline void findMatch (int32_t* m1,const int32_t &i1,int32_t* m2,const int32_t &step_size,
                         std::vector<int32_t> *k2,const int32_t &u_bin_num,const int32_t &v_bin_num,const int32_t &stat_bin,
                         int32_t& min_ind,int32_t stage,bool flow,bool use_prior,double u_=-1,double v_=-1);
  void matching (int32_t *m1p,int32_t *m1c,
                 int32_t n1p,int32_t n1c,
                 std::vector<p_match> &p_matched,bool use_prior,Viso::Matrix *Tr_delta = 0);

  // outlier removal
  void removeOutliers (std::vector<p_match> &p_matched);

  // parabolic fitting
  bool parabolicFitting(const uint8_t* I1_du,const uint8_t* I1_dv,const int32_t* dims1,
                        const uint8_t* I2_du,const uint8_t* I2_dv,const int32_t* dims2,
                        const float &u1,const float &v1,
                        float       &u2,float       &v2,
                        Viso::Matrix At,Viso::Matrix AtA,
                        uint8_t* desc_buffer);
  void relocateMinimum(const uint8_t* I1_du,const uint8_t* I1_dv,const int32_t* dims1,
                       const uint8_t* I2_du,const uint8_t* I2_dv,const int32_t* dims2,
                       const float &u1,const float &v1,
                       float       &u2,float       &v2,
                       uint8_t* desc_buffer);
  void refinement (std::vector<p_match> &p_matched);

  // mean for gain computation
  inline float mean(const uint8_t* I,const int32_t &bpl,const int32_t &u_min,const int32_t &u_max,const int32_t &v_min,const int32_t &v_max);

  // parameters
  parameters param;
  int32_t    margin;
  
  int32_t *m1p1,*m1c1;//m1p : features in previous frame, m1p1: previous frame first pass (bigger nms), m1p2  
  int32_t *m1p2,*m1c2;
  int32_t n1p1,n1c1;//n1p1 : number of features in m1p1 array
  int32_t n1p2,n1c2;
  uint8_t *I1p,*I1c;//image precedent left(1) and right(2)
  uint8_t *I1p_du,*I1c_du;
  uint8_t *I1p_dv,*I1c_dv;
  uint8_t *I1p_du_full,*I1c_du_full; // only needed for
  uint8_t *I1p_dv_full,*I1c_dv_full; // half-res matching
  int32_t dims_p[3],dims_c[3];
  
  //for extraction of feature in unwarped image
   int32_t *m1f1;//m1p : features in previous frame, m1p1: previous frame first pass (bigger nms), m1p2  
  int32_t *m1f2;
  int32_t n1f1;//n1p1 : number of features in m1p1 array
  int32_t n1f2;
  uint8_t *I1f;//image precedent left(1) and right(2)
  uint8_t *I1f_du;
  uint8_t *I1f_dv;
  uint8_t *I1f_du_full; // only needed for
  uint8_t *I1f_dv_full; // half-res matching
  int32_t dims_f[3];
 
   std::vector<int> LUTtoIDwarp;

   int i1pFirstKltFeat;
   std::vector<p_feat> KltFeatRef;//to keep precise location of feature
   
   int i1cFirstKltFeat;
   std::vector<p_feat> KltFeatCurr;//to keep precise location of feature
   
  std::vector<p_match> p_matched_1;
  std::vector<p_match> p_matched_2;
  std::vector<Matcher::range>   ranges;
};

#endif

