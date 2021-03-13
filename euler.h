#ifndef EULER_H
#define EULER_H

#include "mesh.h"
#include "cell.h"
#include <math.h>
#include <utility>



const double gama = 1.4;

class Euler {
public:
    Euler(Mesh &M) : m_mesh(M){};
    Mesh &m_mesh;
    void Update(double time);

private:
    //void solvePas(double dt);
    //void cfl() const;
    Etat F(Etat W);

    Etat G(Etat W);

    Etat FHLL(Etat WL, Etat WR, double lamin, double lamax);

    Etat GHLL(Etat WD, Etat WU, double lamin, double lamax);

    Etat FUPWIND(Etat WLD, Etat WLU, Etat WRD, Etat WRU);

    Etat GUPWIND(Etat WLD, Etat WLU, Etat WRD, Etat WRU);

    std::pair<Etat,Etat>  HLLMULTID(Etat WLD, Etat WLU, Etat WRD, Etat WRU);

    void Solve();
    void Pressure(int cas);
    void Bord();
    double CFL();
    double SRD,SLD, SRU, SLU, SDR, SDL, SUR, SUL;  // vitesses d'ondes
    double SU, SD, SL, SR;
    double T{0.000000001};
    double Tmax{0.};
    double dt{0.};
    std::vector<double> Wave;
};


#endif
