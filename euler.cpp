#include "euler.h"
#include "cell.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <algorithm>
//#include <armadillo>
#include <stdlib.h>
#include <eigen3/Eigen/Dense>
Etat Euler::F(const Etat& W) {

    Etat Fx;

    Fx.rho = W.u;
    Fx.u =  pow(W.u, 2)/W.rho + W.p;
    Fx.v =W.u/W.rho * W.v;
    Fx.E = (W.u /W.rho)*(W.E + W.p );

    return Fx;
}


Etat Euler::G(const  Etat& W) {

    Etat Gy;

    Gy.rho = W.v;
    Gy.u =  W.u/W.rho * W.v;
    Gy.v = pow(W.v, 2)/W.rho + W.p;
    Gy.E = ( W.v/W.rho)* (W.E + W.p );

    return Gy;
}

//---------------------------STAR FLUXES-------------------------------------------------------------

Etat Euler::F_star(const Etat& W,const Etat& GW){
    //GW.p=0.4 * (GW.E -0.5 *GW.rho *(pow(GW.u /GW.rho,2)+ pow(GW.v /GW.rho, 2)));
  //  std::cout<<" rho "<<W.rho<<" .u "<<W.u<<" p  "<<W.p <<" E  "<<W.E<<std::endl;

    Etat Fx;
    if(W.rho !=0 && W.v !=0) {
        Fx.rho = W.u;
        Fx.u = G(W).v + (pow(W.u / W.rho, 2) - pow(W.v / W.rho, 2));
        Fx.v = W.u * W.v / W.rho;
        Fx.E =W.u * G(W).E / W.v;
    }
    else{Fx=GW;}

    return Fx;
}


Etat Euler::G_star(const Etat& W, const Etat& FW) {
    //FW.p=0.4 * (FW.E -0.5 *FW.rho *(pow(FW.u /FW.rho,2)+ pow(FW.v /FW.rho, 2)));
   Etat Gy;

    if(W.rho !=0 && W.u !=0) {
    Gy.rho = W.v;
    Gy.u =  W.u * W.v/W.rho;
    Gy.v =F(W).u+ (pow(W.v, 2)-pow(W.u, 2))/W.rho;
    Gy.E =  W.v*F(W).E/W.u;
    }
    else{Gy=FW;}
    return Gy;
}

//-----------------------------------------PRESSURE UPDATE--------------------------------------------------------------

void Euler::Pressure(){
#pragma omp parallel for collapse(2)
            for (unsigned int it = 0; it < m_mesh.m_nx; it++) {
                for (int j = 0; j < m_mesh.m_ny; ++j) {
                       m_mesh.Primal[it][j].next.p = std::max(0., 0.4 * (m_mesh.Primal[it][j].next.E -0.5 * m_mesh.Primal[it][j].next.rho *
                       (pow(m_mesh.Primal[it][j].next.u / m_mesh.Primal[it][j].next.rho, 2) +pow(m_mesh.Primal[it][j].next.v /m_mesh.Primal[it][j].next.rho, 2))));
                }
            }
}


Etat Euler::FHLL(const Etat &WL, const Etat &WR, double lamin, double lamax) {
    Etat WF;
 // lamin =std::min(WR.u/WR.rho-sqrt(1.4*WR.p/WR.rho),WL.u / WL.rho - sqrt(gama * WL.p / WL.rho));
 // lamax =std::max(WR.u/WR.rho+sqrt(1.4*WR.p/WR.rho),WL.u / WL.rho + sqrt(gama * WL.p / WL.rho));
   lamin =WL.u / WL.rho - sqrt(gama * WL.p / WL.rho);
    lamax =WR.u/WR.rho+sqrt(1.4*WR.p/WR.rho);

     if (lamin < 0. && 0. < lamax)
    {
       WF=(lamax*F(WL)-lamin*F(WR)+lamin*lamax*(WR-WL))/(lamax-lamin);

    }

   if (0. <= lamin) { WF = F(WL);     }
   if (lamax <= 0.) { WF = F(WR);     }
   return WF;
}


Etat Euler::GHLL(const Etat &WD, const Etat &WU, double lamin, double lamax) {
    Etat WG;
   // lamin= std::min(WU.v/ WU.rho  - sqrt(gama * WU.p / WU.rho),WD.v/ WD.rho  -sqrt(gama * WD.p / WD.rho));
   // lamax=std::max(WU.v/ WU.rho + sqrt(gama * WU.p / WU.rho),WD.v/ WD.rho  +sqrt(gama * WD.p / WD.rho));
    lamin= WD.v/ WD.rho  -sqrt(gama * WD.p / WD.rho);
    lamax=WU.v/ WU.rho + sqrt(gama * WU.p / WU.rho);

    if (lamin < 0. && 0. < lamax)
    {
        WG=(lamax*G(WD)-lamin*G(WU)+lamin*lamax*(WU-WD))/(lamax-lamin);

    }
    if (0. <= lamin) {  WG=G(WD); }
    if (lamax <= 0.) {  WG=G(WU); }
    return WG;
}

// when flow is supersonic in both x & y use upwind

std::pair<Etat, Etat> Euler::UPWIND(const Etat &WLD, const Etat &WLU, const Etat &WRD, const Etat &WRU) {
    Etat FUP,GUP;

//-------------------RIGHT PB-------------------------------
    SDR =WRD.v/ WRD.rho  - sqrt(gama * WRD.p / WRD.rho);
    SUR = WRU.v/ WRU.rho + sqrt(gama * WRU.p / WRU.rho);
//-------------------LEFT PB--------------------------------
    SDL = WLD.v/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho);
    SUL = WLU.v/ WLU.rho  + sqrt(gama * WLU.p / WLU.rho);

//--------------------BOTTOM PB-----------------------------
    SLD = WLD.u/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho);
    SRD = WRD.u/ WRD.rho + sqrt(gama * WRD.p / WRD.rho);
//-------------------TOP PB---------------------------------
    SRU =WRU.u/ WRU.rho  + sqrt(gama * WRU.p / WRU.rho);
    SLU =WLU.u/ WLU.rho - sqrt(gama * WLU.p / WLU.rho);


    SR = std::max(SRD,SRU);
    SL = std::min(SLU, SLD);
    SU = std::max(SUR, SUL);
    SD = std::min(SDR, SDL);



    if (0. <= SL && 0. <= SD)
    {
            FUP=F(WLD);
            GUP=G(WLD);
    }

    else if(0. <= SL && SU <= 0.)
    {
            FUP=F(WLU);
            GUP=G(WLU);
    }


    else if (SR <= 0. && 0. <= SD)
    {
        FUP=F(WRD);
        GUP=G(WRD);
    }
    else if(SR <= 0. && SU <= 0.)
    {
        FUP=F(WRU);
        GUP=G(WRU);
    }
    else{

      FUP=0.5*(FHLL(WLD,WRD)+FHLL(WLU,WRU));
      GUP=0.5*(GHLL(WLD,WLU)+GHLL(WRD,WRU));
      // FUP=0.25 * (F(WLU) + F(WRU))-0.25 * m_mesh.getdx() / dt * (WRU - WLU)+0.25 * (F(WLD) + F(WRD))-0.25 * m_mesh.getdx() / dt * (WRD - WLD);
       //GUP= 0.25 * (G(WLD) + G(WLU))-0.25 * m_mesh.getdy() / dt *(WLU- WLD)+ 0.25 * (G(WRD) + G(WRU))-0.25 * m_mesh.getdy() / dt *(WRU -WRD);

    }
    return {FUP, GUP};
}

//--------------------------------------------MULTIDIMENSIONNAL HLL---------------------------



std::pair<Etat, Etat> Euler::HLL_MULTID(const Etat &WLD, const Etat &WLU, const Etat &WRD, const Etat &WRU)
{

    Etat Wstar, Fstar, Gstar;
    double u_r,c_r,u_l, c_l, u_bot, c_bot, u_top, c_top;
    std::vector<double> tmp;
    
    Etat UD, UR, UU, UL, FUD, FUU, GUL, GUR;
//---------------------------------WAVES SPEEDS-------------
//-------------------RIGHT PB-------------------------------
    u_r=(WRD.v*sqrt(WRD.rho)+WRU.v*sqrt(WRU.rho))/(sqrt(WRD.rho)+sqrt(WRU.rho));
    c_r=sqrt( ((gama*WRD.p/sqrt(WRD.rho))+(gama*WRU.p/sqrt(WRU.rho)))/(sqrt(WRD.rho)+sqrt(WRU.rho))+\
    0.5*(sqrt(WRD.rho*WRU.rho)/pow(sqrt( WRU.rho)+sqrt( WRD.rho),2))*pow(WRU.v-WRD.v,2));
    
    SRD = std::min(WRD.v/ WRD.rho - sqrt(gama * WRD.p / WRD.rho), u_r-c_r);
    SRU = std::max(WRU.v/ WRU.rho + sqrt(gama * WRU.p / WRU.rho), u_r+c_r);
    
//-------------------LEFT PB--------------------------------

    u_l=(WLD.v*sqrt(WLD.rho)+WLU.v*sqrt(WLU.rho))/(sqrt(WLD.rho)+sqrt(WLU.rho));
    c_l=sqrt( ((gama*WLD.p/sqrt(WLD.rho))+(gama*WLU.p/sqrt(WLU.rho)))/(sqrt(WLD.rho)+sqrt(WLU.rho))+\
    0.5*(sqrt(WLD.rho*WLU.rho)/pow(sqrt( WLD.rho)+sqrt( WLD.rho),2))*pow(WLU.v-WLD.v,2));

    SLD = std::min(WLD.v/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho), u_l-c_l);
    SLU = std::max(WLU.v/ WLU.rho  + sqrt(gama * WLU.p / WLU.rho), u_l+c_l);
    
//--------------------BOT PB--------------------------------
       
     u_bot=(WLD.u*sqrt( WLD.rho)+WRD.u*sqrt( WRD.rho))/(sqrt( WLD.rho)+sqrt( WRD.rho));
     c_bot=sqrt((gama * WLD.p /sqrt(WLD.rho)+(gama * WRD.p / sqrt(WRD.rho)))/(sqrt( WLD.rho)+sqrt( WRD.rho))+\
    0.5*((sqrt( WLD.rho*WRD.rho))/pow(sqrt( WLD.rho)+sqrt( WRD.rho),2))*pow(WRD.u-WLD.u,2));
     
    SDL =std::min( WLD.u/ WLD.rho - sqrt(gama * WLD.p / WLD.rho), u_bot-c_bot);
    SDR =std::max( WRD.u/ WRD.rho + sqrt(gama * WRD.p / WRD.rho), u_bot+c_bot) ;

//-------------------TOP PB---------------------------------
        
    u_top=(WLU.u*sqrt( WLU.rho)+WRU.u*sqrt( WRU.rho))/(sqrt( WLU.rho)+sqrt( WRU.rho));
    c_top=sqrt((gama * WLU.p /sqrt(WLU.rho)+(gama * WRU.p / sqrt(WRU.rho)))/(sqrt( WLU.rho)+sqrt( WRU.rho))+\
    0.5*((sqrt( WLU.rho*WRU.rho))/pow(sqrt( WLU.rho)+sqrt( WRU.rho),2))*pow(WRU.u-WLU.u,2));
    
    SUL =std::min( WLU.u/ WLU.rho - sqrt(gama * WLU.p / WLU.rho), u_top-c_top);
    SUR =std::max( WRU.u/ WRU.rho + sqrt(gama * WRU.p / WRU.rho), u_top+c_top);
    
//-------------------------------------------------------------

    SR = std::max(SRU, SRD);
    SL = std::min(SLD, SLU);
    SU = std::max(SUL, SUR);
    SD = std::min(SDL, SDR);
	tmp = {abs(SD),abs(SU),abs(SL),abs(SR)};
    const auto max = std::max_element(begin(tmp), end(tmp));
    Wave.push_back(*max);
#pragma omp parallel num_threads(10) shared(Wstar,Fstar, Gstar)

    if (SL <= 0. && 0. <= SR && SD <= 0. &&  0. <= SU) // strong interaction state straddles time axis

    {
      // ---------------------------------BALSARA APPROXIMATION
 //--------------------------------------------------------------------------------------------------------------------

    Wstar =  (SR * SU * WRU + SL * SD * WLD - SR * SD * WRD - SL * SU * WLU)/( (SR - SL) * (SU - SD))
            -  (SU * (F(WRU) - F(WLU)) - SD * (F(WRD) - F(WLD)) + SR * (G(WRU) - G(WRD)) - SL * (G(WLU) - G(WLD)))/((SR - SL) * (SU - SD))
            + 0.5*(SRU * (F(WRU) - F_star(UR,GUR)) - SRD * (F(WRD) - F_star(UR,GUR)) - SLU * (F(WLU) - F_star(UL,GUL))
            + SLD * (F(WLD) - F_star(UL,GUL)))/( (SR - SL) * (SU - SD))
            +0.5* (SUR * (G(WRU) - G_star(UU,FUU)) - SUL * (G(WLU) - G_star(UU,FUU)) - SDR * (G(WRD) - G_star(UD,FUD)) +
            SDL * (G(WLD) - G_star(UD,FUD)))/( (SR - SL) * (SU - SD));

//--------------------------------------------------------------------------------------------------------------------

    Fstar= 2. * SR * Wstar - (SU * FHLL(WLU, WRU)-SD* FHLL(WLD, WRD))/ (SU - SD)
           + 2.  * (SU * F(WRU) - SD * F(WRD) + SR * (G(WRU) - G(WRD)) - SR * (SU * WRU - SD * WRD))/ (SU - SD)
           + (SRD * (F(WRD) - F_star(UR,GUR)) - SRU * (F(WRU) - F_star(UR,GUR)) - std::max(SUR, 0.) * (G(WRU) - G_star(UU,FUU))
              + std::max(SUL, 0.) * (G(WLU) - G_star(UU,FUU))+ std::max(SDR, 0.) * (G(WRD) - G_star(UD,FUD))
              - std::max(SDL, 0.) * (G(WLD) - G_star(UD,FUD)))/(SU - SD);

    Gstar=2. * SU * Wstar -(SR * GHLL(WRD, WRU) -SL* GHLL(WLD, WLU))/ (SR - SL)
          +  2. * (SR * G(WRU) - SL * G(WLU) + SU * (F(WRU) - F(WLU)) - SU * (SR * WRU - SL * WLU))/ (SR - SL)
          +   (SUL * (G(WLU) - G_star(UU,FUU)) - SUR * (G(WRU) - G_star(UU,FUU)) - std::max(SRU, 0.) * (F(WRU) - F_star(UR,GUR))
          + std::max(SRD, 0.) * (F(WRD) - F_star(UR,GUR)) + std::max(SLU, 0.)* (F(WLU) - F_star(UL,GUL))
          - std::max(SLD, 0.) * (F(WLD) -F_star(UL,GUL)))/ (SR - SL);
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
       
       
    }
   else if (((SL >= 0. && SR >= 0.) || (SL <= 0. && SR <= 0.)) &&((SU >= 0. && SD >= 0.) || (SU <= 0. && SD <= 0.)))
    {// both x,y supersonic
        Fstar = UPWIND(WLD, WLU, WRD, WRU).first;
        Gstar = UPWIND(WLD, WLU, WRD, WRU).second;
    }
    else if ((SD <= 0. && SU >= 0.) && (SL >= 0. || SR <= 0.))
    { //x supersonic, y subsonic
        Fstar = (SU * FHLL(WLD, WRD) - SD * FHLL(WLU, WRU))/ (SU - SD);
        if (SL >= 0.)
        {            Gstar = GHLL(WLD, WLU);}
        else if (SR <= 0.)
        {         Gstar = GHLL(WRD, WRU); }
    }
    else if ((SL <= 0. && SR >= 0.) && (SD >= 0. || SU <= 0.))
    {//x subsonic, y supersonic
            Gstar =(SR * GHLL(WLD, WLU) - SL * GHLL(WRD, WRU)) / (SR - SL);

        if (SD >= 0.)
        {           Fstar =FHLL(WLD, WRD);}
        else if (SU <= 0.)
        {            Fstar =FHLL(WLU, WRU);  }
    }
#pragma omp barrier
    return {Fstar, Gstar};
}


void Euler::Solve() {

    Etat FNW, FSW, FSE, FNE, GNW, GSW, GSE, GNE;


    for (unsigned it = 1; it < m_mesh.m_nx-1 ; it++) {
        for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j) {
//--------------------------------------INTERNAL CELLS---------------------------------------------------------------
//

            if(m_mesh.Primal[it][j].bord==0) {  //-------INTERNAL CELLS-----
                    FNW = HLL_MULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;

                    FSW = HLL_MULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                    m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;

                    FSE = HLL_MULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                    m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;

                    FNE = HLL_MULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                    m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;

                    GNW = HLL_MULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;

                    GSW = HLL_MULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                    m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;
                    GSE = HLL_MULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                    m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;
                    GNE = HLL_MULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                    m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;


                    m_mesh.Primal[it][j].FR =
                          (1./6.)* (FNE + FSE)+(4./6.)*FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                    m_mesh.Primal[it][j].FL =
                            (1./6.)*(FNW + FSW)+(4./6.)*FHLL(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it ][j].prev);;
                    m_mesh.Primal[it][j].GU =
                            (1./6.)* (GNW + GNE)+(4./6.)*GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);
                    m_mesh.Primal[it][j].GD =
                           (1./6.)*(GSW + GSE)+(4./6.)*GHLL(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j].prev);

                }
//-----------------------------------------BOTTOM BORDER-------------------------------------------------------------
           if(m_mesh.Primal[it][j].bord==1) {

                FNW = HLL_MULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                FNE = HLL_MULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;
                GNW = HLL_MULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                GNE = HLL_MULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;

                FSE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);

                FSW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it + 1][j].prev);

                //-------------------
                m_mesh.Primal[it][j].GU = 0.25 * (GNW + GNE);
                m_mesh.Primal[it][j].FR = 0.25 * (FNE + FSE);
                m_mesh.Primal[it][j].FL = 0.25 * (FNW + FSW);

                m_mesh.Primal[it][j].GD.rho = 0.;
                m_mesh.Primal[it][j].GD.u = 0.;
                m_mesh.Primal[it][j].GD.v = 0.;
                m_mesh.Primal[it][j].GD.E = 0.;
            }
//-------------------------------------------RIGHT BORDER------------------------------------------------------------
            if(m_mesh.Primal[it][j].bord==2) {

                    FNW = HLL_MULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                    FSW = HLL_MULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                    m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;

                    GNW = HLL_MULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                    GSW = HLL_MULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                    m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;

                    GSE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);


                    GNE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                    m_mesh.Primal[it][j].FL = 0.25 * (FNW + FSW);
                    m_mesh.Primal[it][j].GU = 0.25 * (GNW + GNE);
                    m_mesh.Primal[it][j].GD = 0.25 * (GSW + GSE);

                    m_mesh.Primal[it][j].FR.rho = 0.;
                    m_mesh.Primal[it][j].FR.u = 0.;
                    m_mesh.Primal[it][j].FR.v = 0.;
                    m_mesh.Primal[it][j].FR.E = 0.;
                }
//--------------------------------------------TOP BORDER-------------------------------------------------------------
            if(m_mesh.Primal[it][j].bord==3) {

                FSW = HLL_MULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;
                FSE = HLL_MULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;

                GSW = HLL_MULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;
                GSE = HLL_MULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;

                FNW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);

                FNE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);

                m_mesh.Primal[it][j].FR = 0.25 * (FNE + FSE);
                m_mesh.Primal[it][j].FL = 0.25 * (FNW + FSW);
                m_mesh.Primal[it][j].GD = 0.25 * (GSW + GSE);

                m_mesh.Primal[it][j].GU.rho = 0.;
                m_mesh.Primal[it][j].GU.u = 0.;
                m_mesh.Primal[it][j].GU.v = 0.;
                m_mesh.Primal[it][j].GU.E = 0.;
            }
//-------------------------------------------LEFT BORDER-------------------------------------------------------------
            if(m_mesh.Primal[it][j].bord==4) {

                FSE = HLL_MULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;
                FNE = HLL_MULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;


                GSE = HLL_MULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;
                GNE = HLL_MULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;
                GNW = GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);

                GSW = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                m_mesh.Primal[it][j].FR = 0.25 * (FNE + FSE);
                m_mesh.Primal[it][j].GU = 0.25 * (GNW + GNE);
                m_mesh.Primal[it][j].GD = 0.25 * (GSW + GSE);


                m_mesh.Primal[it][j].FL.rho = 0.;
                m_mesh.Primal[it][j].FL.u = 0.;
                m_mesh.Primal[it][j].FL.v = 0.;
                m_mesh.Primal[it][j].FL.E = 0.;
            }

        }

    }


}


void Euler::HLL_Solver(double time, int testcase) {
alpha=0.4;
    Tmax = time;
    m_mesh.save(0,"hll");

   while (T < Tmax) {
   
        Solve();
        dt = CFL();
        if (dt == 0.) {
            std::cout << "space step " << m_mesh.getdx() << std::endl;
            std::cout << "time step " << dt << ' ' << T << std::endl;
           //break;
        }
        Update();
        T += dt;
  }

    m_mesh.save(1,"hll");
}

double Euler::CFL() {
    
    std::vector<double> tmp;
    double step;
    double   xmid,xmax,ymid,ymax;

    for (unsigned it =0; it < m_mesh.m_nx-1; it++) {
        for (unsigned j =0; j < m_mesh.m_ny-1; j++) {


            xmid = abs(m_mesh.Primal[it][j].prev.u/ m_mesh.Primal[it][j].prev.rho) +sqrt(1.4* m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho) ;                      ;
            ymid = abs(m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho)+sqrt(1.4* m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho);
            xmax =abs(m_mesh.Primal[it+1][j].prev.u/ m_mesh.Primal[it+1][j].prev.rho )+sqrt(1.4* m_mesh.Primal[it+1][j].prev.p / m_mesh.Primal[it+1][j].prev.rho);
            ymax= abs(m_mesh.Primal[it][j+1].prev.v / m_mesh.Primal[it][j+1].prev.rho) +sqrt(1.4* m_mesh.Primal[it][j+1].prev.p / m_mesh.Primal[it][j+1].prev.rho);

            tmp = {xmid,xmax,ymid,ymax};
            const auto max = std::max_element(begin(tmp), end(tmp));
            Wave.push_back(*max);
        }
    }

    const auto max = std::max_element(begin(Wave), end(Wave));
    step =alpha*m_mesh.getdx() / *max;
    Wave.clear();
    return step;
}

void Euler::Bord() {
#pragma omp parallel for
    for (unsigned it = 0; it < m_mesh.m_nx; it++) // BORD HORIZONTAUX
    {
        m_mesh.Primal[it][0].next = m_mesh.Primal[it][m_mesh.m_ny -2].prev;// index -2 ??
        m_mesh.Primal[it][m_mesh.m_ny - 1].next = m_mesh.Primal[it][1].prev;
    }
#pragma omp parallel for
 for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j)// BORD VERTICAUX
    {
        m_mesh.Primal[0][j].next = m_mesh.Primal[m_mesh.m_nx - 2][j].prev;
        m_mesh.Primal[m_mesh.m_nx - 1][j].next = m_mesh.Primal[1][j].prev;
    }


}
//---------------------------------------------------------------------------------------------------------------------
//                                            HLL_split SOLVER-
//---------------------------------------------------------------------------------------------------------------------

void Euler::HLL_SPLIT(double time, int testcase) {
alpha=0.25;

    Tmax = time;
    m_mesh.save(0,"split_hll");


   while (T < Tmax) {
        dt = CFL();
        if (dt == 0.) {
            std::cout << "time step NUL "  << " Current Time " << T << "  space step " << m_mesh.getdx() <<std::endl;
           break;
        }
#pragma omp parallel num_threads(8)
#pragma omp parallel for collapse(2) schedule(dynamic)
       for (unsigned it = 1; it < m_mesh.m_nx-1 ; ++it) {
            for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j) {

            if(m_mesh.Primal[it][j].bord==0) {

                //---------------------------------------INTERNAL CELLS-----------------------------------------
                m_mesh.Primal[it][j].FR =FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                m_mesh.Primal[it][j].FL =FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);
                m_mesh.Primal[it][j].GU = GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);
                m_mesh.Primal[it][j].GD = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

            }

                        //-----------------------------------------BOTTOM BORDER-------------------------------------------------------------
            else if(m_mesh.Primal[it][j].bord==1) {
                m_mesh.Primal[it][j].FR = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                m_mesh.Primal[it][j].FL = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);

                m_mesh.Primal[it][j].GU.rho = 0.;
                m_mesh.Primal[it][j].GU.u = 0.;
                m_mesh.Primal[it][j].GU.v = 0.;
                m_mesh.Primal[it][j].GU.E = 0.;

                m_mesh.Primal[it][j].GD.rho = 0.;
                m_mesh.Primal[it][j].GD.u = 0.;
                m_mesh.Primal[it][j].GD.v = 0.;
                m_mesh.Primal[it][j].GD.E = 0.;
            }
                        //-------------------------------------------RIGHT BORDER---------------------------------------
            else if(m_mesh.Primal[it][j].bord==2) {
                m_mesh.Primal[it][j].GU = GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);
                m_mesh.Primal[it][j].GD = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                m_mesh.Primal[it][j].FR.rho = 0.;
                m_mesh.Primal[it][j].FR.u = 0.;
                m_mesh.Primal[it][j].FR.v = 0.;
                m_mesh.Primal[it][j].FR.E = 0.;

                m_mesh.Primal[it][j].FL.rho = 0.;
                m_mesh.Primal[it][j].FL.u = 0.;
                m_mesh.Primal[it][j].FL.v = 0.;
                m_mesh.Primal[it][j].FL.E = 0.;
            }

                        //--------------------------------------------TOP BORDER----------------------------------------
            else if(m_mesh.Primal[it][j].bord==3) {
                m_mesh.Primal[it][j].FR = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                m_mesh.Primal[it][j].FL = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);

                m_mesh.Primal[it][j].GU.rho = 0.;
                m_mesh.Primal[it][j].GU.u = 0.;
                m_mesh.Primal[it][j].GU.v = 0.;
                m_mesh.Primal[it][j].GU.E = 0.;

                m_mesh.Primal[it][j].GD.rho = 0.;
                m_mesh.Primal[it][j].GD.u = 0.;
                m_mesh.Primal[it][j].GD.v = 0.;
                m_mesh.Primal[it][j].GD.E = 0.;
            }
                        //-------------------------------------------LEFT BORDER----------------------------------------
            else if(m_mesh.Primal[it][j].bord==4){
                        m_mesh.Primal[it][j].GU =GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);
                        m_mesh.Primal[it][j].GD = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                        m_mesh.Primal[it][j].FR.rho =0.;
                        m_mesh.Primal[it][j].FR.u =0.;
                        m_mesh.Primal[it][j].FR.v =0.;
                        m_mesh.Primal[it][j].FR.E =0.;

                        m_mesh.Primal[it][j].FL.rho =0.;
                        m_mesh.Primal[it][j].FL.u =0.;
                        m_mesh.Primal[it][j].FL.v =0.;
                        m_mesh.Primal[it][j].FL.E =0.;
                }
            }
        }
        Update();
        T +=dt;
   }
    m_mesh.save(1,"split_hll");
}

//---------------------------------------------------------------------------------------------------------------------

void Euler::LF_Solver(double time, int testcase) {
alpha=0.7;
    Tmax = time;
    m_mesh.save(0,"LAX");
 while (T < Tmax) {
        dt = CFL();

       if (dt == 0.) {
          //  std::cout << "space step " << m_mesh.getdx() << std::endl;
            std::cout << "time step " << dt << ' ' << T << std::endl;
            //break;
        }

#pragma omp parallel for collapse(2) num_threads(8)
        for (unsigned it = 1; it < m_mesh.m_nx - 1; it++) {
            for (unsigned j = 1; j < m_mesh.m_ny - 1; ++j) {
                m_mesh.Primal[it][j].FR = 0.25 * (F(m_mesh.Primal[it][j].prev) + F(m_mesh.Primal[it + 1][j].prev))
                -0.25 * m_mesh.getdx() / dt * (m_mesh.Primal[it + 1][j].prev - m_mesh.Primal[it][j].prev);
                m_mesh.Primal[it][j].FL = 0.25 * (F(m_mesh.Primal[it][j].prev) + F(m_mesh.Primal[it - 1][j].prev))
                - 0.25 * m_mesh.getdx() / dt *(m_mesh.Primal[it][j].prev - m_mesh.Primal[it - 1][j].prev);
                m_mesh.Primal[it][j].GU =  0.25 * (G(m_mesh.Primal[it][j].prev) + G(m_mesh.Primal[it][j + 1].prev))
                -0.25 * m_mesh.getdy() / dt *(m_mesh.Primal[it][j + 1].prev - m_mesh.Primal[it][j].prev);
                m_mesh.Primal[it][j].GD = 0.25 * (G(m_mesh.Primal[it][j].prev) + G(m_mesh.Primal[it][j - 1].prev))
                -0.25 * m_mesh.getdy() / dt *(m_mesh.Primal[it][j].prev - m_mesh.Primal[it][j - 1].prev);
            }
        }
Update();
        T += dt;
 }

    m_mesh.save(1,"lax");
}
//----------------------------------------------------------------------------------------------------------------------

void Euler::UPWIND_Solver(double time, int testcase){
    alpha=0.5;
    Tmax = time;
    m_mesh.save(0,"upwind");
    while (T < Tmax) {
        dt = CFL();

        if (dt == 0.) {
            std::cout << "time step " << dt << ' ' << T << std::endl;
            break;
        }

#pragma omp parallel num_threads(8)
#pragma omp parallel for collapse(2)
        for (unsigned it = 1; it < m_mesh.m_nx - 1; ++it) {
            for (unsigned j = 1; j < m_mesh.m_ny - 1; ++j) {
                m_mesh.Primal[it][j].FR =0.25*UPWIND(m_mesh.Primal[it][j].prev, m_mesh.Primal[it ][j+1].prev,m_mesh.Primal[it+1 ][j].prev,m_mesh.Primal[it+1 ][j+1].prev).first
                        +0.25*UPWIND(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it ][j].prev,m_mesh.Primal[it+1 ][j-1].prev,m_mesh.Primal[it+1 ][j].prev).first;

                m_mesh.Primal[it][j].FL = 0.25*UPWIND(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it-1][j+1].prev,m_mesh.Primal[it ][j].prev,m_mesh.Primal[it ][j+1].prev).first
                        +0.25*UPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev).first;

                m_mesh.Primal[it][j].GU =0.25*UPWIND(m_mesh.Primal[it][j].prev, m_mesh.Primal[it ][j+1].prev,m_mesh.Primal[it+1 ][j].prev,m_mesh.Primal[it+1 ][j+1].prev).second
                        +0.25*UPWIND(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it-1 ][j+1].prev,m_mesh.Primal[it ][j].prev,m_mesh.Primal[it ][j+1].prev).second;


                m_mesh.Primal[it][j].GD =0.25*UPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev).second
                        +0.25*UPWIND(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it ][j].prev,m_mesh.Primal[it+1 ][j-1].prev,m_mesh.Primal[it+1 ][j].prev).second;
            }
        }
        Update();
        T += dt;
    }

    m_mesh.save(1,"upwind");
}

//----------------------------------------------------------------------------------------------------------------------
double Euler::minmod(double a, double b)
{
if(a*b<=0.) return 0.;
   else  return std::min(a,b);
}
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
void Euler::Update(){

#pragma omp parallel for collapse(2) num_threads(8)
    for (unsigned it =1; it < m_mesh.m_nx-1 ; ++it) {
        for (unsigned j = 1; j < m_mesh.m_ny-1; ++j) {

                    m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev
                            -dt / m_mesh.getdx() *(m_mesh.Primal[it][j].FR -m_mesh.Primal[it][j].FL)
                            -dt / m_mesh.getdy() *(m_mesh.Primal[it][j].GU - m_mesh.Primal[it][j].GD);
        }
    }
 
    Bord();
    Pressure();



#pragma omp parallel for collapse(2) num_threads(8)
    for (unsigned it = 0; it < m_mesh.m_nx; it++) {
        for (unsigned j =0; j < m_mesh.m_ny; j++) {
            m_mesh.Primal[it][j].prev = m_mesh.Primal[it][j].next;
        }
    }
}
