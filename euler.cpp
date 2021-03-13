#include "euler.h"
#include "mesh.h"
#include "cell.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm>

Etat Euler::F(Etat W) {

    Etat Fx;

    Fx.rho = W.u;
    Fx.u = (1 / W.rho) * pow(W.u, 2) + W.p;
    Fx.v = (1 / W.rho) * (W.u * W.v);
    Fx.E = (1 / W.rho) * (W.E * W.u + W.p * W.u);

    return Fx;
}

Etat Euler::G(Etat W) {

    Etat Gy;

    Gy.rho = W.v;
    Gy.u = (1 / W.rho) * W.u * W.v;
    Gy.v = (1 / W.rho) * pow(W.v, 2) + W.p;
    Gy.E = (1 / W.rho) * (W.E * W.v + W.p * W.v);

    return Gy;
}

void Euler::Pressure(int cas) {

    switch (cas) {
        case 0:
            for (unsigned int it = 0; it < m_mesh.m_nx; it++) {
                for (int j = 0; j < m_mesh.m_ny ; ++j) {
                    m_mesh.Primal[it][j].prev.p =
                            (gama - 1) * (m_mesh.Primal[it][j].prev.E \
 - 0.5 * ((pow(m_mesh.Primal[it][j].prev.u, 2) + pow(m_mesh.Primal[it][j].prev.v, 2)) /
          pow(m_mesh.Primal[it][j].prev.rho, 2)));
                }
            }
        case 1:
            for (unsigned int it = 0; it < m_mesh.m_nx; it++) {
                for (int j = 0; j < m_mesh.m_ny; ++j) {
                    m_mesh.Primal[it][j].next.p =
                            (gama - 1) * (m_mesh.Primal[it][j].next.E - 0.5 * pow(m_mesh.Primal[it][j].next.u, 2) + pow(m_mesh.Primal[it][j].next.v, 2));
                }
            }
    }

}

Etat Euler::FHLL(Etat WL, Etat WR, double lamin, double lamax) {
    Etat Wstar;

    if (0 < lamin) {
        Wstar = WL;
    } else if (lamin < 0 && 0 < lamax) {

        Wstar = 1. / (lamax - lamin) * (lamax * WR - lamin * WL + F(WL) - F(WR));
    } else if (0 > lamax) {

        Wstar = WR;
    }
    return Wstar;
}


Etat Euler::GHLL(Etat WD, Etat WU, double lamin, double lamax) {
    Etat Wstar;

    if (0 < lamin) { Wstar = WD; }
    else if (lamin < 0 && 0 < lamax) {
        Wstar = 1. / (lamax - lamin) * (lamax * WU - lamin * WD + G(WD) - G(WU));
    } else if (0 > lamax) { Wstar = WU; }
    return Wstar;
}

// when flow is supersonic in both x & y use upwind

Etat Euler::FUPWIND(Etat WLD, Etat WLU, Etat WRD, Etat WRU) {
    Etat R;

    if (SL >= 0. && SD >= 0.) { R = F(WLD); }
    else if (SL >= 0. && SU <= 0.) { R = F(WLU); }
    else if (SR <= 0. && SD >= 0.) { R = F(WRD); }
    else if (SR <= 0. && SU >= 0.) { R = F(WRU); }
    return R;
}

Etat Euler::GUPWIND(Etat WLD, Etat WLU, Etat WRD, Etat WRU) {
    Etat R;

    if (SL >= 0. && SD >= 0.) { R = G(WLD); }
    else if (SL >= 0. && SU <= 0.) { R = G(WLU); }
    else if (SR <= 0. && SD >= 0.) { R = G(WRD); }
    else if (SR <= 0. && SU >= 0.) { R = G(WRU); }
    return R;
}

//----------------------MULTIDIMENSIONNAL HLL---------------------------

std::pair<Etat, Etat> Euler::HLLMULTID(Etat WLD, Etat WLU, Etat WRD, Etat WRU) {

    Etat Wstar;
    Etat Fstar;
    Etat Gstar;
    Etat UD, UR, UU, UL;
    //------------------------------WAVES SPEEDS--------------------------------------------

//-------------------RIGHT PB-------------------------------
    SRD = WRD.v / WRD.rho - sqrt(gama * WRD.p / WRD.rho);
    SRU = WRU.v / WRU.rho + sqrt(gama * WRU.p / WRU.rho);
//-------------------LEFT PB--------------------------------
    SLD = WLD.v / WLD.rho - sqrt(gama * WLD.p / WLD.rho);
    SLU = WLU.v / WLU.rho + sqrt(gama * WLU.p / WLU.rho);
//--------------------BOT PB--------------------------------
    SDR = WRD.u / WRD.rho + sqrt(gama * WRD.p / WRD.rho);
    SDL = WLD.u / WRD.rho - sqrt(gama * WLD.p / WLD.rho);
//-------------------TOP PB---------------------------------
    SUR = WRU.u / WRU.rho + sqrt(gama * WRU.p / WRU.rho);
    SUL = WLU.u / WLU.rho - sqrt(gama * WLU.p / WLU.rho);

    SR = std::max(SDR, SUR);
    SL = std::min(SDL, SUL);
    SU = std::max(SLU, SRU);
    SD = std::min(SRD, SLD);


//----------------------------------------STAR STATES---------------------------------
    UD = (1 / (SDR - SDL)) * (SDR * WRD - SDL * WLD + F(WLD) - F(WRD));
    UR = (1 / (SRU - SRD)) * (SRU * WRU - SRD * WRD + F(WRD) - F(WRU));
    UU = (1 / (SLU - SLD)) * (SLU * WLU - SLD * WLD + F(WLD) - F(WLU));
    UL = (1 / (SRU - SRD)) * (SRU * WRU - SRD * WRD + F(WRD) - F(WRU));
//--------------------------------------------------------------------------------------------------------------------

    Wstar = (1. / ((SR - SL) * (SU - SD))) * (SR * SU * WRU + SL * SD * WLD - SR * SD * WRD - SL * SU * WLU) \
 - (1. / (2 * (SR - SL) * (SU - SD))) *
   (SU * (F(WRU) - F(WLU)) - SD * (F(WRD) - F(WLD)) + SR * (G(WRU) - G(WRD)) - SL * (G(WLU) - G(WLD))) \
 + (SRU * (F(WRU) - F(UR)) - SRD * (F(WRD) - F(UR)) - SLU * (F(WLU) - F(UL)) + SLD * (F(WLD) - F(UL))) \
 + (1. / (2 * (SR - SL) * (SU - SD))) *
   (SUR * (G(WRU) - G(UU)) - SUL * (G(WLU) - G(UU)) - SDR * (G(WRD) - G(UD)) + SDL * (G(WRD) - G(UD)));

//--------------------------------------------------------------------------------------------------------------------
    if (SL < 0. && SR > 0. && SD < 0. && SU > 0.) { // straddles time axis
        Fstar = 2 * SR * Wstar - (SU / (SU - SD)) * FHLL(WLU, WRU, SUL, SUR) +
                (SD / (SU - SD)) * FHLL(WLD, WLU, SDL, SDR) \
 + (2 / (SU - SD)) * (SU * F(WRU) - SD * F(WRD) + SR * (G(WRU) - G(WRD)) - SR * (SU * WRU - SD * WRD)) \
 + (1 / (SU - SD)) *
   (std::max(SRD, 0.) * (F(WRD) - F(UR)) - SRU * (F(WRU) - F(UR)) - std::max(SUR, 0.) * (G(WRU) - G(UU)) +
    std::max(SLU, 0.) * (G(WLU) - G(UU)) +
    std::max(SDR, 0.) * (G(WRD) - G(UD)) - std::max(SLD, 0.) * (G(WLD) - G(UD)));

        Gstar = 2 * SU * Wstar - (SR / (SR - SL)) * GHLL(WRD, WRU, SRD, SRU) +
                (SL / (SR - SL)) * GHLL(WLD, WLU, SLD, SLU) \
 + (2 / (SR - SL)) * (SR * G(WRU) - SL * G(WLU) + SU * (F(WRU) - F(WLU)) - SU * (SR * WRU - SL * WLU)) \
 + (1 / (SR - SL)) *
   (SLU * (G(WLU) - G(UU)) - SUR * (G(WRU) - G(UU)) - std::max(SRU, 0.) * (F(WRU) - F(UR)) +
    std::max(SRD, 0.) * (F(WRD) - F(UR)) +
    std::max(SUL, 0.) * (F(WLD) - F(UL)) + std::max(SLD, 0.) * (F(WLD) - F(UL)));
    } else if (((SL > 0. && SR > 0.) || (SL < 0. && SR < 0.)) &&
               ((SU > 0. && SD > 0.) || (SU < 0. && SD < 0.)))// both x,y supersonic
    {
        Fstar = FUPWIND(WLD, WLU, WRD, WRU);
        Gstar = GUPWIND(WLD, WLU, WRD, WRU);

    } else if ((SD < 0. && SU > 0.) && (SL > 0. || SR < 0.)) { //x supersonic, y subsonic
        Fstar = (1. / (SU - SD)) * (SU * FHLL(WLD, WRD, SDL, SDR) - SD * FHLL(WLU, WRU, SUL, SUR));
        if (SL > 0.) {
            Gstar = GHLL(WLU, WRU, SUL, SUR);
        } else if (SR < 0.) {
            Gstar = GHLL(WRD, WLD, SDL, SDR);
        }
    } else if ((SL < 0. && SR > 0.) && (SD > 0. || SU < 0.)) {//x subsonic, y supersonic
        Gstar = (1. / (SR - SL)) * (SR * GHLL(WLD, WLU, SLD, SLU) + SL * GHLL(WRD, WRU, SRD, SRU));

        if (SD > 0.) {
            Fstar = FHLL(WLD, WRD, SDL, SDR);
        } else if (SU < 0.) {
            Fstar = FHLL(WLU, WRU, SUL, SUR);
        }
    }

    return {Fstar, Gstar};
}


void Euler::Solve() {


    int border;
    Etat FNW, FSW, FSE, FNE;
    Etat GNW, GSW, GSE, GNE;
    int nx(m_mesh.m_nx);
    int ny(m_mesh.m_ny);
    double lspeed, rspeed, dspeed, uspeed;

    for (unsigned it = 0; it < m_mesh.m_nx-1 ; it++) {
        for (unsigned j = 0; j < m_mesh.m_ny-1 ; ++j) {

            border = m_mesh.Primal[it][j].bord;

            switch (border) {//-------INTERNAL CELLS-----

                case 0:
                    FNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev, \
                  m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                    FSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev, \
                  m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;
                    FSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, \
                  m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;
                    FNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, \
                  m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;

                    GNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev, \
                  m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                    GSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev, \
                  m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;
                    GSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, \
                  m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;
                    GNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, \
                  m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;


                    m_mesh.Primal[it][j].FR = 0.5 * (FNE + FSE);
                    m_mesh.Primal[it][j].FL = 0.5 * (FNW + FSW);
                    m_mesh.Primal[it][j].GU = 0.5 * (GNW + GNE);
                    m_mesh.Primal[it][j].GD = 0.5 * (GSW + GSE);
//-----------------------------------------BOTTOM BORDER-------------------------------------------------------------
                case 1:

                    FNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev, \
                  m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                    FNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, \
                  m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;
                    GNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev, \
                  m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                    GNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, \
                  m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;

                    //-----------
                    lspeed = std::min(0., m_mesh.Primal[it][j].prev.u / m_mesh.Primal[it][j].prev.rho -
                                          sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho));
                    rspeed = std::max(0., m_mesh.Primal[it + 1][j].prev.u / m_mesh.Primal[it + 1][j].prev.rho +
                                          sqrt(gama * m_mesh.Primal[it + 1][j].prev.p /
                                               m_mesh.Primal[it + 1][j].prev.rho));
                    FSE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev, lspeed, rspeed);


                    //-----------
                    lspeed = std::min(0., m_mesh.Primal[it - 1][j].prev.u / m_mesh.Primal[it - 1][j].prev.rho -
                                          sqrt(gama * m_mesh.Primal[it - 1][j].prev.p /
                                               m_mesh.Primal[it - 1][j].prev.rho));
                    rspeed = std::max(0., m_mesh.Primal[it][j].prev.u / m_mesh.Primal[it - 1][j].prev.rho +
                                          sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho));

                    FSW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it + 1][j].prev, lspeed, rspeed);

                    //-------------------
                    m_mesh.Primal[it][j].GU = 0.5 * (GNW + GNE);
                    m_mesh.Primal[it][j].FR = 0.5 * (FNE + FSE);
                    m_mesh.Primal[it][j].FL = 0.5 * (FNW + FSW);

//-------------------------------------------RIGHT BORDER-------------------------------------------------------------
                case 2:
                    if (it != 0 && j != 0) {
                        FNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev, \
                          m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).first;
                        FSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev, \
                     m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;

                        GNW = HLLMULTID(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it - 1][j + 1].prev, \
                       m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev).second;
                        GSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev, \
                       m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;

                        dspeed = std::min(0., m_mesh.Primal[it][j - 1].prev.v / m_mesh.Primal[it][j - 1].prev.rho -
                                              sqrt(gama * m_mesh.Primal[it][j - 1].prev.p /
                                                   m_mesh.Primal[it][j - 1].prev.rho));
                        uspeed = std::max(0., m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho +
                                              sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho));
                        GSE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, dspeed, uspeed);

                        dspeed = std::min(0., m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho -
                                              sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho));
                        uspeed = std::max(0., m_mesh.Primal[it][j + 1].prev.v / m_mesh.Primal[it][j + 1].prev.rho +
                                              sqrt(gama * m_mesh.Primal[it][j + 1].prev.p /
                                                   m_mesh.Primal[it][j + 1].prev.rho));

                        GNE = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, dspeed, uspeed);


                        //m_mesh.Primal[it][j].FR = 0.5 * (FNE + FSE);
                        m_mesh.Primal[it][j].FL = 0.5 * (FNW + FSW);
                        m_mesh.Primal[it][j].GU = 0.5 * (GNW + GNE);
                        m_mesh.Primal[it][j].GD = 0.5 * (GSW + GSE);
                    }
                    //--------------------------------------------TOP BORDER-------------------------------------------------------------
                case 3:
                    if (it != 0 && j != 0) {
                        FSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev, \
                       m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).first;
                        FSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, \
                       m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;

                        GSW = HLLMULTID(m_mesh.Primal[it - 1][j - 1].prev, m_mesh.Primal[it - 1][j].prev, \
                          m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev).second;
                        GSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, \
                          m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;

                        //-----------
                        lspeed = std::min(0., m_mesh.Primal[it - 1][j].prev.u / m_mesh.Primal[it - 1][j].prev.rho -
                                              sqrt(gama * m_mesh.Primal[it - 1][j].prev.p /
                                                   m_mesh.Primal[it - 1][j].prev.rho));
                        rspeed = std::max(0., m_mesh.Primal[it][j].prev.u / m_mesh.Primal[it][j].prev.rho +
                                              sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho));
                        FNW = FHLL(m_mesh.Primal[it - 1][j].prev, m_mesh.Primal[it][j].prev, lspeed, rspeed);


                        //--------------
                        lspeed = std::min(0., m_mesh.Primal[it][j].prev.u / m_mesh.Primal[it][j].prev.rho -
                                              sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho));
                        rspeed = std::max(0., m_mesh.Primal[it + 1][j].prev.u / m_mesh.Primal[it + 1][j].prev.rho +
                                              sqrt(gama * m_mesh.Primal[it + 1][j].prev.p /
                                                   m_mesh.Primal[it + 1][j].prev.rho));

                        FNE = FHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it + 1][j].prev, lspeed, rspeed);

                        m_mesh.Primal[it][j].FR = 0.5 * (FNE + FSE);
                        m_mesh.Primal[it][j].FL = 0.5 * (FNW + FSW);
                        // m_mesh.Primal[it][j].GU = 0.5 * (GNW + GNE);
                        m_mesh.Primal[it][j].GD = 0.5 * (GSW + GSE);
                    }
                    //-------------------------------------------LEFT BORDER-------------------------------------------------------------
                case 4:
                    if (it != 0 && j != 0) {
                        FSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, \
                           m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).first;
                        FNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, \
                           m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).first;


                        GSE = HLLMULTID(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, \
                       m_mesh.Primal[it + 1][j - 1].prev, m_mesh.Primal[it + 1][j].prev).second;
                        GNE = HLLMULTID(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, \
                       m_mesh.Primal[it + 1][j].prev, m_mesh.Primal[it + 1][j + 1].prev).second;
                        dspeed = std::min(0., m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho -
                                              sqrt(gama * m_mesh.Primal[it][j].prev.p /
                                                   m_mesh.Primal[it][j].prev.rho));
                        uspeed = std::max(0., m_mesh.Primal[it][j + 1].prev.v / m_mesh.Primal[it][j + 1].prev.rho +
                                              sqrt(gama * m_mesh.Primal[it][j + 1].prev.p /
                                                   m_mesh.Primal[it][j + 1].prev.rho));
                        GNW = GHLL(m_mesh.Primal[it][j].prev, m_mesh.Primal[it][j + 1].prev, dspeed, uspeed);

                        dspeed = std::min(0., m_mesh.Primal[it][j - 1].prev.v / m_mesh.Primal[it][j - 1].prev.rho -
                                              sqrt(gama * m_mesh.Primal[it][j - 1].prev.p /
                                                   m_mesh.Primal[it][j - 1].prev.rho));
                        uspeed = std::max(0., m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho +
                                              sqrt(gama * m_mesh.Primal[it][j].prev.p /
                                                   m_mesh.Primal[it][j].prev.rho));
                        GSW = GHLL(m_mesh.Primal[it][j - 1].prev, m_mesh.Primal[it][j].prev, dspeed, uspeed);

                        m_mesh.Primal[it][j].FR = 0.5 * (FNE + FSE);
                        //  m_mesh.Primal[it][j].FL = 0.5 * (FNW + FSW);
                        m_mesh.Primal[it][j].GU = 0.5 * (GNW + GNE);
                        m_mesh.Primal[it][j].GD = 0.5 * (GSW + GSE);
                    }
            }
        }

    }


}


void Euler::Update(double time) {
    int border;
    Tmax = time;
    m_mesh.uinit(4);

    while (T < Tmax) {
        dt=CFL();

        if (dt == 0.) {
            std::cout << "time step " << dt << ' ' << T << std::endl;
            break;
        }
        Wave.clear();
        Solve();


        for (unsigned it = 1; it < m_mesh.m_nx-1; it++) {
            for (unsigned j = 1; j < m_mesh.m_ny-1; ++j) {

                border = m_mesh.Primal[it][j].bord;
                switch (border) {
                    case 0: // INTERNAL CELLS
                        m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev - \
                        dt / m_mesh.getdx() * (m_mesh.Primal[it][j].FR - m_mesh.Primal[it][j].FL) + \
                        dt / m_mesh.getdy() * (m_mesh.Primal[it][j].GU - m_mesh.Primal[it][j].GD);
                    case 1: // BOTTOM BORDER
                        m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev - dt / m_mesh.getdx() \
 * (m_mesh.Primal[it][j].FR - m_mesh.Primal[it][j].FL) /*- dt / m_mesh.getdy() * m_mesh.Primal[it][j].GU*/;
                    case 2: // RIGHT BORDER
                        m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev + dt / m_mesh.getdx() \
 * m_mesh.Primal[it][j].FL - dt / m_mesh.getdx() * (m_mesh.Primal[it][j].GU - m_mesh.Primal[it][j].GD);
                    case 3: // TOP BORDER
                        m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev - dt / m_mesh.getdx() \
 * (m_mesh.Primal[it][j].FR - m_mesh.Primal[it][j].FL) + dt / m_mesh.getdx() * m_mesh.Primal[it][j].GD;
                    case 4: // LEFT BORDER
                        m_mesh.Primal[it][j].next = m_mesh.Primal[it][j].prev - dt / m_mesh.getdx() \
 * m_mesh.Primal[it][j].FR - dt / m_mesh.getdx() * (m_mesh.Primal[it][j].GU - m_mesh.Primal[it][j].GD);

                }

            }
        }

        Bord();
        Pressure(1);

        T += dt;
    }

    m_mesh.save(1);
}

double Euler::CFL() {
    std::vector<double> tmp;
    double step;
    double horizmin, horizmax, vertimin, vertimax;
    for (unsigned it = 0; it < m_mesh.m_nx; it++) {
        for (unsigned j = 0; j < m_mesh.m_ny; ++j) {
            horizmin = m_mesh.Primal[it][j].prev.u / m_mesh.Primal[it][j].prev.rho -
                       sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho);
            vertimin = m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho -
                       sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho);
            horizmax = m_mesh.Primal[it][j].prev.u / m_mesh.Primal[it][j].prev.rho +
                       sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho);
            vertimax = m_mesh.Primal[it][j].prev.v / m_mesh.Primal[it][j].prev.rho +
                       sqrt(gama * m_mesh.Primal[it][j].prev.p / m_mesh.Primal[it][j].prev.rho);
            //
            tmp = {std::abs(horizmin), std::abs(horizmax), std::abs(vertimin), std::abs(vertimax)};
            const auto max = std::max_element(begin(tmp), end(tmp));
            Wave.push_back(*max);
        }
    }

    const auto max = std::max_element(begin(Wave), end(Wave));
    step =  0.5*m_mesh.getdx()/ *max;
    return step;
}

void Euler::Bord()
{

    for (unsigned it = 0; it < m_mesh.m_nx; it++) // BORD HORIZONTAUX
    {
        m_mesh.Primal[it][0].next = m_mesh.Primal[it][m_mesh.m_ny - 2].prev;
        m_mesh.Primal[it][m_mesh.m_ny - 1].next = m_mesh.Primal[it][1].prev;
    }
    for (unsigned j = 1; j < m_mesh.m_ny-1; ++j)// BORD VERTICAUX
    {
        m_mesh.Primal[0][j].next = m_mesh.Primal[m_mesh.m_nx-2][j].prev;
        m_mesh.Primal[m_mesh.m_nx-1][j].next = m_mesh.Primal[1][j].prev;
    }
}