#pragma once


#define NB_LEVELS 4 // nb of level in pyramid

#define lvl_Viso 1 //level of pyramid used to do continuous tracking
//check in matcher.h if viso uses half resolution of provided image, 
//ie if lvl_viso= 1 and Matcher::parameters::half_resolution == 1 => uses lvl 2 or pyramid, using lvl 2 seems to work well for continuous
//tracking, with enough precision for fundamental matrix estimation and small enough displacement of features for robust matching

//minimum number of matches to define an edge
#define MIN_MATCHES_EDGE 30

//default gain in non linear optimisations
#define DEFAULT_GAIN_GN 0.3

//maximum number of KF used in active window
#define MAX_KF_ACTIVE 2
