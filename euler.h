#ifndef EULER_H
#define EULER_H

#include "mesh.h"
#include "cell.h"
#include <math.h>
#include <utility>
#include <algorithm>

const double gama = 1.4;

class Euler {
public:
    Euler(Mesh &M) : m_mesh(M){};


    Mesh &m_mesh;
    void HLL_Solver(double time,int testcase);
    void LF_Solver(double time, int testcase);
    void HLL_SPLIT(double time, int testcase);
    void UPWIND_Solver(double time, int testcase);
private:

    Etat F(const Etat& W);

    Etat F_star(const Etat& W);
    Etat G_star(const Etat& W);

    Etat G(const Etat& W);

    Etat FHLL(const Etat& WL, const Etat& WR, double lamin=0., double lamax=0.);

    Etat GHLL(const Etat& WD, const Etat& WU,double lamin=0., double lamax=0.);

    Etat FUPWIND(const Etat& WLD, const Etat& WLU, const Etat& WRD, const Etat& WRU) ;

    Etat GUPWIND(const Etat& WLD, const Etat& WLU, const Etat& WRD, const Etat& WRU) ;

    std::pair<Etat,Etat>  HLLMULTID(const Etat& WLD, const Etat& WLU, const Etat& WRD, const Etat& WRU) ;

    void Solve();
    void Update();
    void Pressure();
    void Bord();
    double CFL();
    double minmod(double a, double b);
    void  Muscl();
    double SRD,SLD, SRU, SLU, SDR, SDL, SUR, SUL;  // vitesses d'ondes
    double SU, SD, SL, SR;
    double T{0.};
    double Tmax{0.};
    double dt{0.};
    int solver_choice;
    std::vector<double> Wave;
    std::vector < std::vector<double> > slope_x, slope_y;
};


#endif
