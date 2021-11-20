#ifndef MESH_H
#define MESH_H

#include  <iostream>
#include <vector>
#include "cell.h"


class Mesh {

public:
    Mesh();//constructeur
    Mesh(int i, int j, double L, double H);

    double getdx();
    double getdy();


    void uinit(int choice);

    void save(int a, std::string filename);
    int m_nx; // nombre de points en x et y
    int m_ny;
    std::vector < std::vector<cell_p> > Primal;

private:

    double m_Long;
    double m_Haut;

    double m_dx;
    double m_dy;


};

#endif 
