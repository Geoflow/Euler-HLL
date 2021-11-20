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
    Fx.v = W.u/W.rho * W.v;
    Fx.E = (W.u /W.rho)*(W.E + W.p );

    return Fx;
}


Etat Euler::G(const  Etat& W) {

    Etat Gy;

    Gy.rho = W.v;
    Gy.u =  W.u * W.v/W.rho;;
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
                for (int j = 1; j < m_mesh.m_ny-1; ++j) {
                        m_mesh.Primal[it][j].next.p = 0.4 *(m_mesh.Primal[it][j].next.E -
                        0.5*m_mesh.Primal[it][j].next.rho*(pow(m_mesh.Primal[it][j].next.u/m_mesh.Primal[it][j].next.rho, 2)
                        + pow(m_mesh.Primal[it][j].next.v/m_mesh.Primal[it][j].next.rho, 2)) );

                }
            }
}


Etat Euler::FHLL(const Etat &WL, const Etat &WR, double lamin, double lamax) {
    lamin=WL.u/WL.rho- sqrt(1.4*WL.p/WL.rho);
    lamax=WR.u/WR.rho+sqrt(1.4*WR.p/WR.rho);

    if (lamin < 0. && 0. < lamax)
    {
   // return (lamax*F(WL)-lamin*F(WR)+lamin*lamax*(WR-WL))/(lamax-lamin);
   return 0.5*(F(WR)+F(WL))-0.5*(lamin+lamax)/(lamax-lamin)*(F(WR)+F(WL))+(lamin*lamax)/(lamax-lamin)*(WR-WL);
    }
    else
    {
        if (0. <= lamin) { return F(WL); }
        if (lamax <= 0.) { return F(WR); }
    }
}


Etat Euler::GHLL(const Etat &WD, const Etat &WU, double lamin, double lamax) {

    lamin=WD.v/WD.rho- sqrt(gama*WD.p/WD.rho);
    lamax=WU.v/WU.rho+sqrt(gama*WU.p/WU.rho);


    if (lamin < 0. && 0. < lamax)
    {
     //return (lamax*G(WD)-lamin*G(WU)+lamin*lamax*(WU-WD))/(lamax-lamin);
     return 0.5*(G(WU)+G(WD))-0.5*(lamin+lamax)/(lamax-lamin)*(G(WU)+G(WD))+(lamin*lamax)/(lamax-lamin)*(WU-WD);
    }
    else
    {
        if (0. <= lamin) { return G(WD); }
        if (lamax <= 0.) { return G(WU); }
    }
}

// when flow is supersonic in both x & y use upwind

Etat Euler::FUPWIND(const Etat &WLD, const Etat &WLU, const Etat &WRD, const Etat &WRU) {
//-------------------RIGHT PB-------------------------------
    SRD = WRD.v/ WRD.rho  - sqrt(gama * WRD.p / WRD.rho);
    SRU = WRU.v/ WRU.rho + sqrt(gama * WRU.p / WRU.rho);
//-------------------LEFT PB--------------------------------
    SLD = WLD.v/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho);
    SLU = WLU.v/ WLU.rho  + sqrt(gama * WLU.p / WLU.rho);
//--------------------BOTTOM PB-----------------------------
    SDR = WRD.u/ WRD.rho + sqrt(gama * WRD.p / WRD.rho);
    SDL = WLD.u/ WLD.rho  - sqrt(gama * WLD.p / WLD.rho);
//-------------------TOP PB---------------------------------
    SUR = WRU.u/ WRU.rho  + sqrt(gama * WRU.p / WRU.rho);
    SUL = WLU.u/ WLU.rho - sqrt(gama * WLU.p / WLU.rho);

    SR = std::max(SDR, SUR);
    SL = std::min(SDL, SUL);
    SU = std::max(SLU, SRU);
    SD = std::min(SRD, SLD);

    if (SL >= 0. && SD >= 0.) { return F(WLD); }
    else if (SL >= 0. && SU <= 0.) { return F(WLU); }
    else if (SR <= 0. && SD >= 0.) { return F(WRD); }
    else if (SR <= 0. && SU >= 0.) { return F(WRU); }

}

Etat Euler::GUPWIND(const Etat &WLD, const Etat &WLU, const Etat &WRD, const Etat &WRU) {

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

    SR = std::max(SDR, SUR);
    SL = std::min(SDL, SUL);
    SU = std::max(SLU, SRU);
    SD = std::min(SRD, SLD);

    if (SL >= 0. && SD >= 0.) { return G(WLD); }
    else if (SL >= 0. && SU <= 0.) { return G(WLU); }
    else if (SR <= 0. && SD >= 0.) { return G(WRD); }
    else if (SR <= 0. && SU >= 0.) { return G(WRU); }
}

//--------------------------------------------MULTIDIMENSIONNAL HLL---------------------------

std::pair<Etat, Etat> Euler::HLLMULTID(const Etat &WLD, const Etat &WLU, const Etat &WRD, const Etat &WRU) {

    Etat Wstar;
    Etat Fstar;
    Etat Gstar;
    Etat UD, UR, UU, UL;
//---------------------------------WAVES SPEEDS-------------
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

    SR = std::max(SDR, SUR);
    SL = std::min(SDL, SUL);
    SU = std::max(SLU, SRU);
    SD = std::min(SRD, SLD);




    if (SL <= 0. && SR >= 0. && SD <= 0. && SU >= 0.)
    { // strong interaction state straddles time axis

//----------------------------------------STAR STATES-----------------------------------------------------------------
        UD =  (SDR*WRD-SDL*WLD + F(WRD) - F(WLD))/ (SDR - SDL);
        UD.p=0.4 *UD.E -0.5*((pow(UD.u, 2)+ pow(UD.v, 2))/UD.rho );

        UR =   (SRU*WRU-SLD*WRD +G(WRU) - G(WRD))/(SRU - SRD);
        UR.p=0.4 *UR.E -0.5*((pow(UR.u, 2)+ pow(UR.v, 2))/UR.rho );

        UU = (SUR*WRU-SUL*WLU + F(WRU) - F(WLU))/ (SUR - SUL);
        UU.p=0.4 *UU.E -0.5*((pow(UU.u, 2)+ pow(UU.v, 2))/UU.rho );

        UL = (SLU*WLU-SLD*WLD + G(WLU) - G(WLD)) / (SLU - SLD);
        UL.p=0.4 *UL.E -0.5*((pow(UL.u, 2)+ pow(UL.v, 2))/UL.rho );

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

    } if (((SL >= 0. && SR >= 0.) || (SL <= 0. && SR <= 0.)) &&
               ((SU >= 0. && SD >= 0.) || (SU <= 0. && SD <= 0.)))// both x,y supersonic
    {
        Fstar = FUPWIND(WLD, WLU, WRD, WRU);
        Gstar = GUPWIND(WLD, WLU, WRD, WRU);

    }  if ((SD <= 0. && SU >= 0.) && (SL >= 0. || SR <= 0.))
    { //x supersonic, y subsonic
        Fstar = (SU * FHLL(WLD, WRD) - SD * FHLL(WLU, WRU))/ (SU - SD);
        if (SL > 0.) {
            Gstar = GHLL(WLU, WRU);//GHLL(WLU, WRU, SUL, SUR);
        } else if (SR < 0.) {
            Gstar = GHLL(WRD, WLD);//GHLL(WRD, WLD, SDL, SDR);
        }
    }  if ((SL <= 0. && SR >= 0.) && (SD >= 0. || SU <= 0.))
    {//x subsonic, y supersonic
        Gstar =(SR * GHLL(WLD, WLU) - SL * GHLL(WRD, WRU)) / (SR - SL);

        if (SD > 0.) {
            Fstar =FHLL(WLD, WRD);// FHLL(WLD, WRD, SDL, SDR);
        } else if (SU < 0.) {
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


                    m_mesh.Primal[it][j].FR = 1./6. * (FNE + FSE)+(4./6.)*FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                    m_mesh.Primal[it][j].FL =1./6.* (FNW + FSW)+(4./6.)*FHLL(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it ][j].prev);;
                    m_mesh.Primal[it][j].GU = 1./6.* (GNW + GNE)+(4./6.)*GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);
                    m_mesh.Primal[it][j].GD =1./6.* (GSW + GSE)+(4./6.)*GHLL(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it][j].prev);


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
solver_choice=2;
    Tmax = time;
    m_mesh.save(0,"hll");

    while (T < Tmax) {
    dt = CFL();

        if (dt == 0.) {
            std::cout << "space step " << m_mesh.getdx() << std::endl;
            std::cout << "time step " << dt << ' ' << T << std::endl;
           //break;
        }

        Solve();
        Update();


        T += dt;


   }

    m_mesh.save(1,"hll");
}

double Euler::CFL() {
    Wave.clear();
    std::vector<double> tmp;
    double step,alpha;
    double  horizmin,horizmax,vertimin,vertimax;

        alpha=0.35;//UNIDIMENSIONNAL  HLL CFL0
    for (unsigned it =1; it < m_mesh.m_nx-1; it++) {
        for (unsigned j =1; j < m_mesh.m_ny-1; j++) {

            horizmin = m_mesh.Primal[it-1][j].prev.u/ m_mesh.Primal[it-1][j].prev.rho -sqrt(1.4* m_mesh.Primal[it-1][j].prev.p / m_mesh.Primal[it-1][j].prev.rho) ;                      ;
            vertimin = m_mesh.Primal[it][j-1].prev.v / m_mesh.Primal[it][j-1].prev.rho-sqrt(1.4* m_mesh.Primal[it][j-1].prev.p / m_mesh.Primal[it][j-1].prev.rho);
            horizmax =m_mesh.Primal[it+1][j].prev.u/ m_mesh.Primal[it+1][j].prev.rho +sqrt(1.4* m_mesh.Primal[it+1][j].prev.p / m_mesh.Primal[it+1][j].prev.rho);
            vertimax= m_mesh.Primal[it][j+1].prev.v / m_mesh.Primal[it][j+1].prev.rho +sqrt(1.4* m_mesh.Primal[it][j+1].prev.p / m_mesh.Primal[it][j+1].prev.rho);

            tmp = {abs(horizmin), abs(horizmax), abs(vertimin), abs(vertimax)};
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
    solver_choice=1;
    Tmax = time;
    m_mesh.save(0,"split_hll");


   while (T < Tmax) {
        dt = CFL();
        if (dt == 0.) {
            std::cout << "time step NUL "  << " Current Time " << T << "  space step " << m_mesh.getdx() <<std::endl;
           break;
        }
        Etat GNW, GNE, FNE, FSE,FNW,FSW,GSW, GSE;

#pragma omp parallel for collapse(2)
        for (unsigned it = 1; it < m_mesh.m_nx-1 ; it++) {
            for (unsigned j = 1; j < m_mesh.m_ny-1 ; ++j) {

                switch (m_mesh.Primal[it][j].bord) {

                    case 0:
                        //---------------------------------------INTERNAL CELLS-----------------------------------------

                        m_mesh.Primal[it][j].FR =(4./6.)*FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev)
                        +(1./6.)*FUPWIND(m_mesh.Primal[it][j].prev, m_mesh.Primal[it ][j+1].prev,m_mesh.Primal[it+1 ][j].prev,m_mesh.Primal[it+1 ][j+1].prev)
                        +(1./6.)*FUPWIND(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it ][j].prev,m_mesh.Primal[it+1 ][j-1].prev,m_mesh.Primal[it+1 ][j].prev);

                        m_mesh.Primal[it][j].FL =(4./6.)*FHLL(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it][j].prev)
                         +(1./6.)*FUPWIND(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it-1][j+1].prev,m_mesh.Primal[it ][j].prev,m_mesh.Primal[it ][j+1].prev)
                         +(1./6.)*FUPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev);

                        m_mesh.Primal[it][j].GU =(4./6.)*GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev)
                         +(1./6.)*GUPWIND(m_mesh.Primal[it][j].prev, m_mesh.Primal[it ][j+1].prev,m_mesh.Primal[it+1 ][j].prev,m_mesh.Primal[it+1 ][j+1].prev)
                         +(1./6.)*GUPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev);


                        m_mesh.Primal[it][j].GD =(4./6.)*GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j - 1].prev)
                        +(1./6.)*GUPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev)
                        +(1./6.)*GUPWIND(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it ][j].prev,m_mesh.Primal[it+1 ][j-1].prev,m_mesh.Primal[it+1 ][j].prev);

                        //-----------------------------------------BOTTOM BORDER-------------------------------------------------------------
                    case 1:


                        FSE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);
                        FSW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it ][j].prev);


                        m_mesh.Primal[it][j].FR = FSE;
                        m_mesh.Primal[it][j].FL = FSW;

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

                            GSE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);
                            GNE = GHLL(m_mesh.Primal[it][j ].prev, m_mesh.Primal[it][j+1].prev);



                            m_mesh.Primal[it][j].GU = GNE;
                            m_mesh.Primal[it][j].GD =  GSE;

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


                            FNW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev);
                            FNE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev);

                            m_mesh.Primal[it][j].FR = FNE ;
                            m_mesh.Primal[it][j].FL = FNW ;

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



                            GNW = GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev);

                            GSW = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev);

                            m_mesh.Primal[it][j].GU = GNW ;
                            m_mesh.Primal[it][j].GD = GSW ;

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

//Muscl();
//dt = CFL();
//step+=dt;
//Update();
        T +=dt;


   }
    m_mesh.save(1,"split_hll");
}

//---------------------------------------------------------------------------------------------------------------------
//                                            LAX FRIEDRIECH SOLVER-
//---------------------------------------------------------------------------------------------------------------------

void Euler::LF_Solver(double time, int testcase) {
    solver_choice=0;
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
    solver_choice=0;
    Tmax = time;
    m_mesh.save(0,"upwind");
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
                m_mesh.Primal[it][j].FR =0.25*FUPWIND(m_mesh.Primal[it][j].prev, m_mesh.Primal[it ][j+1].prev,m_mesh.Primal[it+1 ][j].prev,m_mesh.Primal[it+1 ][j+1].prev)
                        +0.25*FUPWIND(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it ][j].prev,m_mesh.Primal[it+1 ][j-1].prev,m_mesh.Primal[it+1 ][j].prev);
                m_mesh.Primal[it][j].FL = 0.25*FUPWIND(m_mesh.Primal[it-1][j].prev, m_mesh.Primal[it-1][j+1].prev,m_mesh.Primal[it ][j].prev,m_mesh.Primal[it ][j+1].prev)
                        +0.25*FUPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev);

                m_mesh.Primal[it][j].GU =0.25*GUPWIND(m_mesh.Primal[it][j].prev, m_mesh.Primal[it ][j+1].prev,m_mesh.Primal[it+1 ][j].prev,m_mesh.Primal[it+1 ][j+1].prev)
                        +0.25*GUPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev);


                m_mesh.Primal[it][j].GD =0.25*GUPWIND(m_mesh.Primal[it-1][j-1].prev, m_mesh.Primal[it-1 ][j].prev,m_mesh.Primal[it ][j-1].prev,m_mesh.Primal[it ][j].prev)
                        +0.25*GUPWIND(m_mesh.Primal[it][j-1].prev, m_mesh.Primal[it ][j].prev,m_mesh.Primal[it+1 ][j-1].prev,m_mesh.Primal[it+1 ][j].prev);
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
    for (unsigned it = 1; it < m_mesh.m_nx-1 ; it++) {
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