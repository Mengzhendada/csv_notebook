#include "FragmentationFunctions.h"
//#include "insane/tmds/FragmentationFunctions.h"

namespace insane {
namespace physics {

//______________________________________________________________________________
FragmentationFunctions::FragmentationFunctions(){
}
//______________________________________________________________________________
FragmentationFunctions::~FragmentationFunctions(){
}
//______________________________________________________________________________
double FragmentationFunctions::D_u_pi0( double z, double Q2)   const {
   double Dplus  = D_u_piplus( z,Q2);
   double Dminus = D_u_piminus(z,Q2);
   return( (Dplus+Dminus)/2.0 );
}
//______________________________________________________________________________
double FragmentationFunctions::D_ubar_pi0( double z, double Q2) const {
   double Dplus  = D_ubar_piplus( z,Q2);
   double Dminus = D_ubar_piminus(z,Q2);
   return( (Dplus+Dminus)/2.0 );
}
//______________________________________________________________________________
double FragmentationFunctions::D_d_pi0( double z, double Q2)    const {
   double Dplus  = D_d_piplus( z,Q2);
   double Dminus = D_d_piminus(z,Q2);
   return( (Dplus+Dminus)/2.0 );
}
//______________________________________________________________________________
double FragmentationFunctions::D_dbar_pi0( double z, double Q2)  const{
   double Dplus  = D_dbar_piplus( z,Q2);
   double Dminus = D_dbar_piminus(z,Q2);
   return( (Dplus+Dminus)/2.0 );
}
//______________________________________________________________________________
double FragmentationFunctions::D_s_pi0( double z, double Q2)   const  {
   double Dplus  = D_s_piplus( z,Q2);
   double Dminus = D_s_piminus(z,Q2);
   return( (Dplus+Dminus)/2.0 );
}
//______________________________________________________________________________
double FragmentationFunctions::D_sbar_pi0( double z, double Q2) const {
   double Dplus  = D_sbar_piplus( z,Q2);
   double Dminus = D_sbar_piminus(z,Q2);
   return( (Dplus+Dminus)/2.0 );
}
//______________________________________________________________________________
}}
