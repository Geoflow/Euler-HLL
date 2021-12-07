#include "euler.h"
#include "cell.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <algorithm>


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
    Gy.u =  W.u/W.rho * W.v;;
    Gy.v = pow(W.v, 2)/W.rho + W.p;
    Gy.E = ( W.v/W.rho)* (W.E + W.p );

    return Gy;
}

//---------------------------STAR FLUXES-------------------------------------------------------------

Etat Euler::F_star(const Etat& W) {

    Etat Fx;

    Fx.rho = W.u;
    Fx.u =  pow(W.v, 2)/W.rho + W.p+ (pow(W.u, 2)-pow(W.v, 2))/W.rho;
    Fx.v = W.u * W.v/W.rho;
    Fx.E = W.u*(W.E + W.p )*(W.u/W.rho)/W.v;

    return Fx;
}


Etat Euler::G_star(const  Etat& W) {

    Etat Gy;

    Gy.rho = W.v;
    Gy.u =  W.u * W.v/W.rho;;
    Gy.v = pow(W.u, 2)/W.rho + W.p+ (pow(W.u, 2)-pow(W.v, 2))/W.rho;
    Gy.E = ( W.v/W.rho)* (W.E + W.p )*(W.v/W.rho)/W.u;

    return Gy;
}

//-----------------------------------------PRESSURE UPDATE--------------------------------------------------------------

void Euler::Pressure(){
#pragma omp parallel for collapse(2)
            for (unsigned int it = 0; it < m_mesh.m_nx; it++) {
                for (int j = 0; j < m_mesh.m_ny; ++j) {
                   if(m_mesh.Primal[it][j].next.rho >=0. && m_mesh.Primal[it][j].next.E>=0.){
                        m_mesh.Primal[it][j].next.p =std::max(0., 0.4 * (m_mesh.Primal[it][j].next.E -0.5 *m_mesh.Primal[it][j].next.rho *
                         (pow(m_mesh.Primal[it][j].next.u /m_mesh.Primal[it][j].next.rho,2)+ pow(m_mesh.Primal[it][j].next.v /m_mesh.Primal[it][j].next.rho, 2))));

                       }else {m_mesh.Primal[it][j].next.p =0.;}

                }
            }
}


Etat Euler::FHLL(const Etat &WL, const Etat &WR, double lamin, double lamax) {
    Etat WF;

        lamin = std::min(WR.u / WR.rho - sqrt(gama * WR.p / WR.rho),
                         WL.u / WL.rho - sqrt(gama * WL.p / WL.rho));//WL.u/WL.rho- sqrt((1.4*WL.p)/WL.rho)
        lamax = std::max(WR.u / WR.rho + sqrt(gama * WR.p / WR.rho),
                         WL.u / WL.rho + sqrt(gama * WL.p / WL.rho));//WR.u/WR.rho+sqrt((1.4*WR.p)/WR.rho);


    if (lamin < 0. && 0. < lamax)
    {
       //  WF=1./(lamax-lamin)*(lamax*F(WL)-lamin*F(WR)+lamin*lamax*(WR-WL));
        WF=0.5*(F(WR)+F(WL))-0.5*(lamin+lamax)/(lamax-lamin)*(F(WR)+F(WL))+(lamin*lamax)/(lamax-lamin)*(WR-WL);
    }

        if (0. <= lamin) { WF = F(WL); }
        if (lamax <= 0.) { WF = F(WR); }

    return WF;
}


Etat Euler::GHLL(const Etat &WD, const Etat &WU, double lamin, double lamax) {
    Etat WG;
    lamin= std::min(WU.v/ WU.rho  - sqrt(gama * WU.p / WU.rho),WD.v/ WD.rho  -sqrt(gama * WD.p / WD.rho));// WD.v/WD.rho- sqrt((1.4*WD.p)/WD.rho);
    lamax=std::max(WU.v/ WU.rho + sqrt(gama * WU.p / WU.rho),WD.v/ WD.rho  +sqrt(gama * WD.p / WD.rho));//WU.v/WU.rho+sqrt((1.4*WU.p)/WU.rho);

    if (lamin < 0. && 0. < lamax)
    {
       // WG=1./(lamax-lamin)*(lamax*G(WD)-lamin*G(WU)+lamin*lamax*(WU-WD));
        WG=0.5*(G(WU)+G(WD))-0.5*(lamin+lamax)/(lamax-lamin)*(G(WU)+G(WD))+(lamin*lamax)/(lamax-lamin)*(WU-WD);

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



    if (0. < SL && 0. < SD)
    {
            FUP=F(WLD);
            GUP=G(WLD);
    }

    else if(0. < SL && SU < 0.)
    {
            FUP=F(WLU);
            GUP=G(WLU);
    }


    else if (SR < 0. && 0. < SD)
    {
        FUP=F(WRD);
        GUP=G(WRD);
    }
    else if(SR < 0. && SU < 0.)
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

std::pair<Etat, Etat> Euler::HLLMULTID(const Etat &WLD, const Etat &WLU, const Etat &WRD, const Etat &WRU) {

    Etat Wstar;
    Etat Fstar;
    Etat Gstar;
    Etat UD, UR, UU, UL;
//---------------------------------WAVES SPEEDS-------------
/*
//-------------------RIGHT PB-------------------------------
    SRD = WRD.v/ WRD.rho  - sqrt(gama * WRD.p / WRD.rho);
    SRU = WRU.v/ WRU.rho + sqrt(gama * WRU.p / WRU.rho);
//-------------------LEFT PB--------------------------------
    SLD = WLD.v/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho);
    SLU = WLU.v/ WLU.rho  + sqrt(gama * WLU.p / WLU.rho);
//--------------------BOT PB--------------------------------
    SDR = WRD.u/ WRD.rho + sqrt(gama * WRD.p / WRD.rho);
    SDL = WLD.u/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho);
//-------------------TOP PB---------------------------------
    SUR = WRU.u/ WRU.rho  + sqrt(gama * WRU.p / WRU.rho);
    SUL = WLU.u/ WLU.rho - sqrt(gama * WLU.p / WLU.rho);
*/
//-------------------RIGHT PB-------------------------------
    SRD =std::min(WRD.v/ WRD.rho  - sqrt(gama * WRD.p / WRD.rho), WRU.v/ WRU.rho - sqrt(gama * WRU.p / WRU.rho));
    SRU =std::max(WRD.v/ WRD.rho  + sqrt(gama * WRD.p / WRD.rho), WRU.v/ WRU.rho + sqrt(gama * WRU.p / WRU.rho));
//-------------------LEFT PB--------------------------------
    SLD =std::min( WLD.v/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho),WLU.v/ WLU.rho  -sqrt(gama * WLU.p / WLU.rho));
    SLU = std::max(WLU.v/ WLU.rho  + sqrt(gama * WLU.p / WLU.rho),WLD.v/ WLD.rho  +sqrt(gama * WLD.p / WLD.rho));
//--------------------BOTTOM PB-----------------------------

    SDL = std::min(WLD.u/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho),WRD.u/ WRD.rho - sqrt(gama * WRD.p / WRD.rho));
    SDR =std::max( WRD.u/ WRD.rho + sqrt(gama * WRD.p / WRD.rho),WLD.u/ WLD.rho +sqrt(gama * WLD.p / WLD.rho));
//-------------------TOP PB---------------------------------
    SUR = std::max(WRU.u/ WRU.rho  + sqrt(gama * WRU.p / WRU.rho),WLU.u/ WLU.rho + sqrt(gama * WLU.p / WLU.rho));
    SUL = std::min(WLU.u/ WLU.rho - sqrt(gama * WLU.p / WLU.rho),WRU.u/ WRU.rho  -sqrt(gama * WRU.p / WRU.rho));


    SR = std::max(SDR, SUR);
    SL = std::min(SDL, SUL);
    SU = std::max(SLU, SRU);
    SD = std::min(SRD, SLD);




    if (SL <= 0. && SR >= 0. && SD <= 0. && SU >= 0.)
    { // strong interaction state straddles time axis

//----------------------------------------STAR STATES-----------------------------------------------------------------
        UD =  (SDR*WRD-SDL*WLD + F(WRD) - F(WLD))/ (SDR - SDL);
        UD.p=0.4 *UD.E -0.5*UD.rho*(pow(UD.u/UD.rho, 2)+ pow(UD.v/UD.rho, 2));

        UR =   (SRU*WRU-SLD*WRD +G(WRU) - G(WRD))/(SRU - SRD);
        UR.p=0.4 *UR.E -0.5*UR.rho*(pow(UR.u/UR.rho, 2)+ pow(UR.v/UR.rho, 2)) ;

        UU = (SUR*WRU-SUL*WLU + F(WRU) - F(WLU))/ (SUR - SUL);
        UU.p=0.4 *UU.E -0.5*UU.rho*(pow(UU.u/UU.rho, 2)+ pow(UU.v/UU.rho, 2));

        UL = (SLU*WLU-SLD*WLD + G(WLU) - G(WLD)) / (SLU - SLD);
        UL.p=0.4 *UL.E -0.5*UL.rho*(pow(UL.u/UL.rho, 2)+ pow(UL.v/UL.rho, 2));

//--------------------------------------------------------------------------------------------------------------------

        Wstar = (1. / ((SR - SL) * (SU - SD))) * (SR * SU * WRU + SL * SD * WLD - SR * SD * WRD - SL * SU * WLU) \
 - (1. / (2 * (SR - SL) * (SU - SD))) *
   (SU * (F(WRU) - F(WLU)) - SD * (F(WRD) - F(WLD)) + SR * (G(WRU) - G(WRD)) - SL * (G(WLU) - G(WLD))) \
 + (SRU * (F(WRU) - F_star(UR)) - SRD * (F(WRD) - F_star(UR)) - SLU * (F(WLU) - F_star(UL)) + SLD * (F(WLD) - F_star(UL))) +
            (1. / (2 * (SR - SL) * (SU - SD)))\
 * (SUR * (G(WRU) - G_star(UU)) - SUL * (G(WLU) - G_star(UU)) - SDR * (G(WRD) - G_star(UD)) + SDL * (G(WRD) - G_star(UD)));

//--------------------------------------------------------------------------------------------------------------------

        Fstar = 2. * SR * Wstar - (SU / (SU - SD)) * FHLL(WLU, WRU) +
                (SD / (SU - SD)) * FHLL(WLD, WLU)
 + (2. / (SU - SD)) * (SU * F(WRU) - SD * F(WRD) + SR * (G(WRU) - G(WRD)) - SR * (SU * WRU - SD * WRD))
 + (1. / (SU - SD)) *(std::max(SRD, 0.) * (F(WRD) - F_star(UR)) - SRU * (F(WRU) - F_star(UR)) - std::max(SUR, 0.) * (G(WRU) - G_star(UU))
 + std::max(SLU, 0.) * (G(WLU) - G_star(UU)) + std::max(SDR, 0.) * (G(WRD) - G_star(UD)) - std::max(SLD, 0.) * (G(WLD) - G_star(UD)));

        Gstar = 2. * SU * Wstar - SR/ (SR - SL) * GHLL(WRD, WRU) +SL/ (SR - SL)* GHLL(WLD, WLU)
 +  (2. / (SU - SD))* (SR * G(WRU) - SL * G(WLU) + SU * (F(WRU) - F(WLU)) - SU * (SR * WRU - SL * WLU))
 + (1. / (SU - SD)) * (SLU * (G(WLU) - G_star(UU)) - SUR * (G(WRU) - G_star(UU)) - std::max(SRU, 0.) * (F(WRU) - F_star(UR))
 + std::max(SRD, 0.) * (F(WRD) - F_star(UR)) + std::max(SLU, 0.) * (F(WLU) - F_star(UL)) - std::max(SLD, 0.) * (F(WLD) - F_star(UL)));

    } if (((SL >= 0. && SR >= 0.) || (SL <= 0. && SR <= 0.)) &&((SU >= 0. && SD >= 0.) || (SU <= 0. && SD <= 0.)))  // both x,y supersonic
    {
        Fstar = UPWIND(WLD, WLU, WRD, WRU).first;
        Gstar = UPWIND(WLD, WLU, WRD, WRU).second;

    }  if ((SD <= 0. && SU >= 0.) && (SL >= 0. || SR <= 0.))
    { //x supersonic, y subsonic
        Fstar = (SU * FHLL(WLD, WRD) - SD * FHLL(WLU, WRU))/ (SU - SD);
        if (SL > 0.)
        {
            Gstar = GHLL(WLU, WRU);//GHLL(WLU, WRU, SUL, SUR);
        } else if (SR < 0.)
        {
            Gstar = GHLL(WRD, WLD);//GHLL(WRD, WLD, SDL, SDR);
        }
    }
    if ((SL <= 0. && SR >= 0.) && (SD >= 0. || SU <= 0.))
    {//x subsonic, y supersonic
        Gstar =(SR * GHLL(WLD, WLU) - SL * GHLL(WRD, WRU)) / (SR - SL);

        if (SD > 0.)
        {
            Fstar =FHLL(WLD, WRD);// FHLL(WLD, WRD, SDL, SDR);
        } else if (SU < 0.)
        {
            Fstar =FHLL(WLU, WRU); //FHLL(WLU, WRU, SUL, SUR);
        }
    }

    return {Fstar, Gstar};
}


void Euler::Solve() {

    Etat FNW, FSW, FSE, FNE, GNW, GSW, GSE, GNE;

#pragma omp parallel num_threads(8)
#pragma omp parallel for collapse(2)
    for (unsigned it = 1; it < m_mesh.m_nx-1 ; it++) {
        for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j) {

            switch (m_mesh.Primal[it][j].bord) {//-------INTERNAL CELLS-----

                case 0:
                    FNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;

                    FSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                    m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;

                    FSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                    m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;

                    FNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                    m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;

                    GNW = HLLMULTID(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it-1][j+1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j+1].prev).second;

                    GSW = HLLMULTID(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it - 1][j].prev,
                                    m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j].prev).second;
                    GSE = HLLMULTID(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j].prev,
                                    m_mesh.Primal[it+1][j-1].prev, m_mesh.Primal[it + 1][j].prev).second;
                    GNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j+1].prev,
                                    m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;


                    m_mesh.Primal[it][j].FR =(1./6.)* (FNE + FSE)+(4./6.)*FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                    m_mesh.Primal[it][j].FL =(1./6.)* (FNW + FSW)+(4./6.)*FHLL(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it ][j].prev);;
                    m_mesh.Primal[it][j].GU =(1./6.)* (GNW + GNE)+(4./6.)*GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);
                    m_mesh.Primal[it][j].GD =(1./6.)* (GSW + GSE)+(4./6.)*GHLL(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j].prev);


//-----------------------------------------BOTTOM BORDER-------------------------------------------------------------
                case 1:

                    FNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                    FNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                    m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;
                    GNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                    m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                    GNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                    m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;

                    FSE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);

                    FSW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it + 1][j].prev);

                    //-------------------
                    m_mesh.Primal[it][j].GU = 0.25 * (GNW + GNE);
                    m_mesh.Primal[it][j].FR = 0.25 * (FNE + FSE);
                    m_mesh.Primal[it][j].FL = 0.25 * (FNW + FSW);

                    m_mesh.Primal[it][j].GD.rho =0.;
                    m_mesh.Primal[it][j].GD.u =0.;
                    m_mesh.Primal[it][j].GD.v =0.;
                    m_mesh.Primal[it][j].GD.E =0.;

//-------------------------------------------RIGHT BORDER-------------------------------------------------------------
                case 2:

                        FNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                        m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                        FSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                        m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;

                        GNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev,
                                        m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                        GSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                        m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;

                        GSE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);


                        GNE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                        m_mesh.Primal[it][j].FL = 0.25 * (FNW + FSW);
                        m_mesh.Primal[it][j].GU = 0.25 * (GNW + GNE);
                        m_mesh.Primal[it][j].GD = 0.25 * (GSW + GSE);

                        m_mesh.Primal[it][j].FR.rho =0.;
                        m_mesh.Primal[it][j].FR.u =0.;
                        m_mesh.Primal[it][j].FR.v =0.;
                        m_mesh.Primal[it][j].FR.E =0.;

                    //--------------------------------------------TOP BORDER-------------------------------------------------------------
                case 3:

                        FSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                        m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;
                        FSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                        m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;

                        GSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev,
                                        m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;
                        GSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                        m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;

                        FNW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);

                        FNE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);

                        m_mesh.Primal[it][j].FR = 0.25 * (FNE + FSE);
                        m_mesh.Primal[it][j].FL = 0.25 * (FNW + FSW);
                        m_mesh.Primal[it][j].GD = 0.25 * (GSW + GSE);

                        m_mesh.Primal[it][j].GU.rho =0.;
                        m_mesh.Primal[it][j].GU.u =0.;
                        m_mesh.Primal[it][j].GU.v =0.;
                        m_mesh.Primal[it][j].GU.E =0.;

                    //-------------------------------------------LEFT BORDER-------------------------------------------------------------
                case 4:

                        FSE = HLLMULTID(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j].prev,
                                        m_mesh.Primal[it+1][j-1].prev, m_mesh.Primal[it + 1][j].prev).first;
                        FNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                        m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;


                        GSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev,
                                        m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;
                        GNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev,
                                        m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;
                        GNW = GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);

                        GSW = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                        m_mesh.Primal[it][j].FR = 0.25 * (FNE + FSE);
                        m_mesh.Primal[it][j].GU = 0.25 * (GNW + GNE);
                        m_mesh.Primal[it][j].GD = 0.25 * (GSW + GSE);


                        m_mesh.Primal[it][j].FL.rho =0.;
                        m_mesh.Primal[it][j].FL.u =0.;
                        m_mesh.Primal[it][j].FL.v =0.;
                        m_mesh.Primal[it][j].FL.E =0.;

            }
        }

    }


}


void Euler::HLL_Solver(double time, int testcase) {
alpha=0.5;
    Tmax = time;
    m_mesh.save(0,"hll");

 //   while (T < Tmax) {
    dt = CFL();

        if (dt == 0.) {
            std::cout << "space step " << m_mesh.getdx() << std::endl;
            std::cout << "time step " << dt << ' ' << T << std::endl;
           //break;
        }

        Solve();
        Update();


        T += dt;


  // }

    m_mesh.save(1,"hll");
}

double Euler::CFL() {
    Wave.clear();
    std::vector<double> tmp;
    double step;
    double   xmid,xmax,ymid,ymax;

      //  alpha=0.125;
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
alpha=0.15;

    Tmax = time;
    m_mesh.save(0,"split_hll");


   while (T < Tmax) {
        dt = CFL();
        if (dt == 0.) {
            std::cout << "time step NUL "  << " Current Time " << T << "  space step " << m_mesh.getdx() <<std::endl;
           break;
        }


        for (unsigned it = 1; it < m_mesh.m_nx-1 ; ++it) {
            for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j) {
                switch (m_mesh.Primal[it][j].bord) {
                    case 0:
                        //---------------------------------------INTERNAL CELLS-----------------------------------------
                        m_mesh.Primal[it][j].FR =(4./6.)*FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev)
                        +(1./6.)*FHLL(m_mesh.Primal[it][j+1].prev, m_mesh.Primal[it + 1][j+1].prev)
                        +(1./6.)*FHLL(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it + 1][j-1].prev);

                        m_mesh.Primal[it][j].FL =(4./6.)*FHLL(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it][j].prev)
                       +(1./6.)*FHLL(m_mesh.Primal[it-1][j+1].prev, m_mesh.Primal[it ][j+1].prev)
                       +(1./6.)*FHLL(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it ][j-1].prev);

                        m_mesh.Primal[it][j].GU =(4./6.)*GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev)
                      +(1./6.)*GHLL(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it-1][j + 1].prev)
                      +(1./6.)*GHLL(m_mesh.Primal[it+1][j].prev, m_mesh.Primal[it+1][j + 1].prev);


                        m_mesh.Primal[it][j].GD =(4./6.)*GHLL(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j ].prev)
                       +(1./6.)*GHLL(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1][j ].prev)
                      +(1./6.)*GHLL(m_mesh.Primal[it+1][j-1].prev, m_mesh.Primal[it+1][j].prev);

                        //-----------------------------------------BOTTOM BORDER-------------------------------------------------------------
                    case 1:
                        m_mesh.Primal[it][j].FR = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                        m_mesh.Primal[it][j].FL = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it ][j].prev);

                        m_mesh.Primal[it][j].GU.rho =0.;
                        m_mesh.Primal[it][j].GU.u =0.;
                        m_mesh.Primal[it][j].GU.v =0.;
                        m_mesh.Primal[it][j].GU.E =0.;

                        m_mesh.Primal[it][j].GD.rho =0.;
                        m_mesh.Primal[it][j].GD.u =0.;
                        m_mesh.Primal[it][j].GD.v =0.;
                        m_mesh.Primal[it][j].GD.E =0.;

                        //-------------------------------------------RIGHT BORDER---------------------------------------
                    case 2:
                        m_mesh.Primal[it][j].GU = GHLL(m_mesh.Primal[it][j ].prev, m_mesh.Primal[it][j+1].prev);
                        m_mesh.Primal[it][j].GD = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                        m_mesh.Primal[it][j].FR.rho =0.;
                        m_mesh.Primal[it][j].FR.u =0.;
                        m_mesh.Primal[it][j].FR.v =0.;
                        m_mesh.Primal[it][j].FR.E =0.;

                        m_mesh.Primal[it][j].FL.rho =0.;
                        m_mesh.Primal[it][j].FL.u =0.;
                        m_mesh.Primal[it][j].FL.v =0.;
                        m_mesh.Primal[it][j].FL.E =0.;

                        //--------------------------------------------TOP BORDER----------------------------------------
                    case 3:
                        m_mesh.Primal[it][j].FR = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                        m_mesh.Primal[it][j].FL = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);

                        m_mesh.Primal[it][j].GU.rho =0.;
                        m_mesh.Primal[it][j].GU.u =0.;
                        m_mesh.Primal[it][j].GU.v =0.;
                        m_mesh.Primal[it][j].GU.E =0.;

                        m_mesh.Primal[it][j].GD.rho =0.;
                        m_mesh.Primal[it][j].GD.u =0.;
                        m_mesh.Primal[it][j].GD.v =0.;
                        m_mesh.Primal[it][j].GD.E =0.;

                        //-------------------------------------------LEFT BORDER----------------------------------------
                    case 4:
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
//                                            LAX FRIEDRIECH SOLVER-
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

#pragma omp parallel for collapse(2)
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
    alpha=0.25;
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
void Euler::Muscl() {

    slope_x.resize(m_mesh.m_nx, std::vector<double>(m_mesh.m_ny));
    slope_y.resize(m_mesh.m_nx, std::vector<double>(m_mesh.m_ny));

#pragma omp parallel for collapse(2)
    for (unsigned it = 1; it < m_mesh.m_nx - 1; it++) {
        for (unsigned j = 1; j < m_mesh.m_ny - 1; j++) {

            slope_x[it][j] = m_mesh.getdx()* 0.5 *minmod(m_mesh.Primal[it][j].prev.u -m_mesh.Primal[it - 1][j].prev.u ,
                    m_mesh.Primal[it + 1][j].prev.u -m_mesh.Primal[it][j].prev.u );
            slope_y[it][j] =m_mesh.getdy()*0.5 *minmod(m_mesh.Primal[it][j].prev.v -m_mesh.Primal[it - 1][j].prev.v,
                                                       m_mesh.Primal[it + 1][j].prev.v -m_mesh.Primal[it][j].prev.v );
        }
    }


#pragma omp parallel for collapse(2)
    for (unsigned it = 1; it < m_mesh.m_nx-1 ; ++it) {
        for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j) {

                    m_mesh.Primal[it][j].FR =FHLL( (1.+slope_x[it][j])*m_mesh.Primal[it][j].prev, (1.-slope_x[it+1][j])*m_mesh.Primal[it+1][j].prev);

                    m_mesh.Primal[it][j].FL =FHLL((1.+slope_x[it-1][j])*m_mesh.Primal[it-1][j].prev, (1.-slope_x[it][j])*m_mesh.Primal[it][j].prev);

                    m_mesh.Primal[it][j].GU =GHLL((1.+slope_y[it][j])*m_mesh.Primal[it][j].prev, (1.-slope_x[it][j+1])*m_mesh.Primal[it][j+1].prev);

                    m_mesh.Primal[it][j].GD =GHLL((1.+slope_y[it][j-1])*m_mesh.Primal[it][j-1].prev, (1.-slope_x[it][j])*m_mesh.Primal[it][j].prev);
        }
    }



}
//----------------------------------------------------------------------------------------------------------------------
void Euler::Update(){

#pragma omp parallel for collapse(2)
    for (unsigned it =1; it < m_mesh.m_nx-1 ; it++) {
        for (unsigned j = 1; j < m_mesh.m_ny-1; ++j) {

                    m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev
                            -dt / m_mesh.getdx() *(m_mesh.Primal[it][j].FR -m_mesh.Primal[it][j].FL)
                            -dt / m_mesh.getdy() *(m_mesh.Primal[it][j].GU - m_mesh.Primal[it][j].GD);
        }
    }
    Bord();
    Pressure();



#pragma omp parallel for collapse(2)
    for (unsigned it = 0; it < m_mesh.m_nx; it++) {
        for (unsigned j =0; j < m_mesh.m_ny; j++) {
            m_mesh.Primal[it][j].prev = m_mesh.Primal[it][j].next;
        }
    }
}