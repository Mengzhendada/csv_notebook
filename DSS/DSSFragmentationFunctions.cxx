#include "DSSFragmentationFunctions.h"
//#include "insane/tmds/DSSFragmentationFunctions.h"
#include "TMath.h"


namespace insane {
namespace physics {

TMutex g_DSS_mutex;

//______________________________________________________________________________
DSSFragmentationFunctions::DSSFragmentationFunctions(){

   // -----------------------------------------------
   // pi+ NLO, from Table I
   fPar_piplus_u_plus_ubar[0] =  0.345;
   fPar_piplus_u_plus_ubar[1] = -0.015;
   fPar_piplus_u_plus_ubar[2] =  1.2;
   fPar_piplus_u_plus_ubar[3] = 11.06;
   fPar_piplus_u_plus_ubar[4] =  4.23;

   fPar_piplus_d_plus_dbar[0] =  0.380;
   fPar_piplus_d_plus_dbar[1] = -0.015;
   fPar_piplus_d_plus_dbar[2] =  1.2;
   fPar_piplus_d_plus_dbar[3] = 11.06;
   fPar_piplus_d_plus_dbar[4] =  4.23;

   fPar_piplus_ubar_eq_d[0] =  0.115;
   fPar_piplus_ubar_eq_d[1] =  0.520;
   fPar_piplus_ubar_eq_d[2] =  3.27;
   fPar_piplus_ubar_eq_d[3] = 16.26;
   fPar_piplus_ubar_eq_d[4] =  8.46;

   fPar_piplus_s_plus_sbar[0] =  0.190;
   fPar_piplus_s_plus_sbar[1] =  0.520;
   fPar_piplus_s_plus_sbar[2] =  3.27;
   fPar_piplus_s_plus_sbar[3] = 16.26;
   fPar_piplus_s_plus_sbar[4] =  8.46;

   fPar_piplus_c_plus_cbar[0] =  0.271;
   fPar_piplus_c_plus_cbar[1] = -0.905;
   fPar_piplus_c_plus_cbar[2] =  3.23;
   fPar_piplus_c_plus_cbar[3] =  0.0;
   fPar_piplus_c_plus_cbar[4] =  0.0;

   fPar_piplus_b_plus_bbar[0] =  0.501;
   fPar_piplus_b_plus_bbar[1] = -1.305;
   fPar_piplus_b_plus_bbar[2] =  5.67;
   fPar_piplus_b_plus_bbar[3] =  0.0;
   fPar_piplus_b_plus_bbar[4] =  0.0;

   fPar_piplus_g[0] =  0.279;
   fPar_piplus_g[1] =  0.899;
   fPar_piplus_g[2] =  1.57;
   fPar_piplus_g[3] = 20.0;
   fPar_piplus_g[4] =  4.91;

   // -----------------------------------------------
   // Kplus NLO, from Table VI
   fPar_Kplus_u_plus_ubar[0] =  0.058;
   fPar_Kplus_u_plus_ubar[1] =  0.705;
   fPar_Kplus_u_plus_ubar[2] =  1.2;
   fPar_Kplus_u_plus_ubar[3] = 15.00;
   fPar_Kplus_u_plus_ubar[4] =  6.02;

   fPar_Kplus_s_plus_sbar[0] =  0.343;
   fPar_Kplus_s_plus_sbar[1] = -0.065;
   fPar_Kplus_s_plus_sbar[2] =  1.2;
   fPar_Kplus_s_plus_sbar[3] =  4.36;
   fPar_Kplus_s_plus_sbar[4] =  3.73;

   fPar_Kplus_d_plus_dbar[0] =  0.016;
   fPar_Kplus_d_plus_dbar[1] =  1.108;
   fPar_Kplus_d_plus_dbar[2] = 10.00;
   fPar_Kplus_d_plus_dbar[3] = 10.00;
   fPar_Kplus_d_plus_dbar[4] =  3.28;

   fPar_Kplus_ubar_eq_s[0] =  0.008;
   fPar_Kplus_ubar_eq_s[1] =  1.108;
   fPar_Kplus_ubar_eq_s[2] = 10.00;
   fPar_Kplus_ubar_eq_s[3] = 10.00;
   fPar_Kplus_ubar_eq_s[4] =  3.28;

   fPar_Kplus_c_plus_cbar[0] =  0.196;
   fPar_Kplus_c_plus_cbar[1] =  0.102;
   fPar_Kplus_c_plus_cbar[2] =  4.56;
   fPar_Kplus_c_plus_cbar[3] =  0.0;
   fPar_Kplus_c_plus_cbar[4] =  0.0;

   fPar_Kplus_b_plus_bbar[0] =  0.139;
   fPar_Kplus_b_plus_bbar[1] = -0.584;
   fPar_Kplus_b_plus_bbar[2] =  7.42;
   fPar_Kplus_b_plus_bbar[3] =  0.0;
   fPar_Kplus_b_plus_bbar[4] =  0.0;

   fPar_Kplus_g[0] =  0.017;
   fPar_Kplus_g[1] =  5.055;
   fPar_Kplus_g[2] =  1.20;
   fPar_Kplus_g[3] =  0.0;
   fPar_Kplus_g[4] =  0.0;

   fNprime = 0.83; // 1.0;

}
//______________________________________________________________________________
DSSFragmentationFunctions::~DSSFragmentationFunctions(){
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_model(double z,double *pars) const {
   // Model function (at input scale)
   double Ni     = pars[0];
   double alphai = pars[1];
   double betai  = pars[2];
   double gammai = pars[3];
   double deltai = pars[4];
   double num = Ni*TMath::Power(z,alphai)*TMath::Power(1.0-z,betai)*(1.0+gammai*TMath::Power(1.0-z,deltai));
   double den = TMath::Beta(2.0+alphai,betai+1.0) + gammai*TMath::Beta(2.0+alphai,betai+deltai+1.0);
   return(num/den);
}
//______________________________________________________________________________

double DSSFragmentationFunctions::D_u_piplus( double z, double Q2)  const   {
  TLockGuard guard(&g_DSS_mutex);
   //double Du_ubar = D_model(z,fPar_piplus_u_plus_ubar);
   //double Dubar   = D_model(z,fPar_piplus_ubar_eq_d);
   //return(Du_ubar - Dubar);
   int ih = 1; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return u;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_u_piminus(double z, double Q2)   const  {
  TLockGuard guard(&g_DSS_mutex);
   // D_u^{pi-} = D_ubar^{pi+}
   //double Du = D_model(z,fPar_piplus_ubar_eq_d);
   //return(Du);
   int ih = 1; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return u;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_ubar_piplus( double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   // D_ubar^{pi+}
   //double Dubar   = D_model(z,fPar_piplus_ubar_eq_d);
   //return(Dubar);
   int ih = 1; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return ub;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_ubar_piminus(double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   // D_ubar^{pi-} = D_u^{pi+}
   //double Du_ubar = D_model(z,fPar_piplus_u_plus_ubar);
   //double Dubar   = D_model(z,fPar_piplus_ubar_eq_d);
   //return(Du_ubar - Dubar);
   int ih = 1; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return ub;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_d_piplus( double z, double Q2)   const  {
  TLockGuard guard(&g_DSS_mutex);
   // D_d^{pi+} 
   //double Dubar   = D_model(z,fPar_piplus_ubar_eq_d);   // D_d^{pi+}
   //return(Dubar);
   int ih = 1; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return d;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_d_piminus(double z, double Q2)   const  {
  TLockGuard guard(&g_DSS_mutex);
   // D_d^{pi-}  = D_dbar^{pi+} = D_{d+dbar} - D_d
   //double Dd_dbar = D_model(z,fPar_piplus_d_plus_dbar);
   //double Dd      = D_model(z,fPar_piplus_ubar_eq_d);
   //return(Dd_dbar - Dd);
   int ih = 1; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return d;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_dbar_piplus( double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   // D_dbar^{pi+} = D_{d+dbar} - D_d
   //double Dd_dbar = D_model(z,fPar_piplus_d_plus_dbar);
   //double Dd      = D_model(z,fPar_piplus_ubar_eq_d);
   //return(Dd_dbar - Dd);
   int ih = 1; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return db;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_dbar_piminus(double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   // D_dbar^{pi-} = D_d^{pi+}
   //double Dd   = D_model(z,fPar_piplus_ubar_eq_d);   // D_d^{pi+}
   //return(Dd);
   int ih = 1; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return db;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_s_piplus( double z, double Q2)   const  {
  TLockGuard guard(&g_DSS_mutex);
   // Eq. 17
   // D_s^{pi+} = D_sbar^{pi+} = N' D_ubar^{pi+}
   //double Ds = fNprime*(this->D_ubar_piplus(z,Q2));
   //return(Ds);
   int ih = 1; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return s;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_s_piminus(double z, double Q2)    const {
  TLockGuard guard(&g_DSS_mutex);
   // Eq. 17
   // D_s^{pi+} = D_sbar^{pi+} = N' D_ubar^{pi+}
   //double Ds = fNprime*(this->D_ubar_piplus(z,Q2));
   //return(Ds);
   int ih = 1; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return s;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_sbar_piplus( double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   // Eq. 17
   // D_s^{pi+} = D_sbar^{pi+} = N' D_ubar^{pi+}
   //double Ds = fNprime*(this->D_ubar_piplus(z,Q2));
   //return(Ds);
   int ih = 1; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return sb;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_sbar_piminus(double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   // Eq. 17
   // D_s^{pi+} = D_sbar^{pi+} = N' D_ubar^{pi+}
   //double Ds = fNprime*(this->D_ubar_piplus(z,Q2));
   //return(Ds);
   int ih = 1; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return sb;
}
//______________________________________________________________________________
double DSSFragmentationFunctions::D_u_Kplus( double z, double Q2)   const  {
  TLockGuard guard(&g_DSS_mutex);
  int ih = 2; //1 = pion, kaon, proton
  int ic = 1; // charge
  int io = 1; // 1 = nlo
  double u, ub, d, db, s, sb, c,b,gl;
  fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
  return u;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_u_Kminus(double z, double Q2)   const  {
  TLockGuard guard(&g_DSS_mutex);
  int ih = 2; //1 = pion, kaon, proton
  int ic = -1; // charge
  int io = 1; // 1 = nlo
  double u, ub, d, db, s, sb, c,b,gl;
  fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
  return u;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_ubar_Kplus( double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return ub;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_ubar_Kminus(double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return ub;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_d_Kplus( double z, double Q2)    const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return d;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_d_Kminus(double z, double Q2)    const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return d;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_dbar_Kplus( double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return db;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_dbar_Kminus(double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return db;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_s_Kplus( double z, double Q2)    const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return s;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_s_Kminus(double z, double Q2)    const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return s;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_sbar_Kplus( double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = 1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return sb;
} 
//______________________________________________________________________________
double DSSFragmentationFunctions::D_sbar_Kminus(double z, double Q2) const {
  TLockGuard guard(&g_DSS_mutex);
   int ih = 2; //1 = pion, kaon, proton
   int ic = -1; // charge
   int io = 1; // 1 = nlo
   double u, ub, d, db, s, sb, c,b,gl;
   fdss_(&ih, &ic, &io, &z, &Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
   return sb;
} 
//______________________________________________________________________________
}
}
