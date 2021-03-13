#include <cstdlib>
#include <iostream>
#include <fstream>
#include "mesh.h"
#include "cell.h"
#include <vector>
#include <cmath>

Mesh::Mesh() {
    m_nx = 10;
    m_ny = 10;
    m_Long = 1.0;
    m_Haut = 1.0;
}


//Mesh :: ~Mesh() {}

Mesh::Mesh(int i, int j, double L, double H) {

    m_nx = i;
    m_ny = j;
    m_Long = L;
    m_Haut = H;
    m_dx = 2 * m_Long / m_nx ;
    m_dy = 2 * m_Haut / m_ny ;

    m_dxd = m_Long / m_nx;
    m_dyd = m_Haut / m_ny;



    Primal.resize(m_nx, std::vector<cell_p>(m_ny));

    //--------------PRIMAL MESH-----------------------------------

    for (int i = 0; i < m_nx ; ++i) {
        for (int j = 0; j < m_ny ; ++j) {

            //INTERNAL CELLS--
            Primal[i][j].bord = 0;
            Primal[i][j].x = -m_Long +(i+0.5) * m_dx;
            Primal[i][j].y = -m_Haut +j * m_dy;


        }
    }
   for (int i = 0; i < m_nx; ++i) {
        // bottom border
        Primal[i][0].bord = 1;
        Primal[i][0].x = -m_Long + (i+0.5) * m_dx;
        Primal[i][0].y = -m_Haut+0.5 *m_dy;

    }
    for (int i = 0; i < m_nx; ++i) {
        //top border
        Primal[i][m_nx - 1].bord = 3;
        Primal[i][m_nx - 1].x = -m_Long +(i+0.5) * m_dx;
        Primal[i][m_nx - 1].y = m_Haut -  0.5 *m_dy;

    }


    for (int j = 0; j < m_ny; ++j) {
        //left border
        Primal[0][j].bord = 4;
        Primal[0][j].x = -m_Long + 0.5 * m_dx;
        Primal[0][j].y = -m_Haut +j * m_dy;

    }
    //right border
    for (int j = 0; j < m_ny; ++j) {
        Primal[m_nx - 1][j].bord = 2;
        Primal[m_nx - 1][j].x = m_Long - 0.5*m_dx;
        Primal[m_nx - 1][j].y = -m_Haut+ j*m_dy;

    }


    Primal[0][0].bord = 5;
    Primal[m_nx - 1][0].bord = 5;
    Primal[0][m_ny-1].bord = 5;
    Primal[m_nx - 1][m_ny-1].bord = 5;


    std::cout << "Primal Mesh Size " << m_nx*m_ny << std::endl;
    for (unsigned int it = 0; it < Primal[1].size(); it++) {
        for (unsigned int j = 0; j < Primal.size(); j++) {
if(Primal[it][j].bord==5) std::cout << Primal[it][j].x << " " << Primal[it][j].y << "  " << Primal[it][j].bord << std::endl;
        }
    }
}





double Mesh::getdx() { return m_dx; }

double Mesh::getdy() { return m_dy; }

void Mesh::uinit(int choice) {

    Etat K;
switch(choice) {
// TEST CASE NUMERO 3
    case 3:
    for (unsigned int it = 0; it < Primal[1].size(); it++) {
        for (unsigned int j = 0; j < Primal.size(); j++) {

            if (Primal[it][j].x > 0.000000001 && Primal[it][j].y > 0.000000001) {// UPPER RIGHT
                K.rho = 1.5;
                K.u = 0.;
                K.v = 0.;
                K.p = 1.5;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;
            }
            if (Primal[it][j].x < 0.000000001 && Primal[it][j].y > 0.000000001) { //UPPER LEFT
                K.rho = 0.5323;
                K.u = 1.206;
                K.v = 0.;
                K.p = 0.3;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;

            }
            if (Primal[it][j].x < 0.000000001 && Primal[it][j].y < 0.000000001) { //LOWER LEFT
                K.rho = 0.138;
                K.u = 1.206;
                K.v = 1.206;
                K.p = 0.029;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;

            }

        if (Primal[it][j].x > 0.000000001 && Primal[it][j].y < 0.000000001) { //LOWER RIGHT
            K.rho = 0.5323;
            K.u = 0.;
            K.v = 1.206;
            K.p = 0.3;
            K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
            Primal[it][j].prev = K;

        }

    }
    }
    case 4:
    //----------------------- TEST CASE NUMERO 4------------------------------------------------------------------------
    for (unsigned int it = 0; it < Primal[1].size(); it++) {
        for (unsigned int j = 0; j < Primal.size(); j++) {

            if (Primal[it][j].x > 0.000000001 && Primal[it][j].y > 0.000000001) {// UPPER RIGHT
                K.rho = 1.1;
                K.u = 0.;
                K.v = 0.;
                K.p = 1.1;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;
            }
            if (Primal[it][j].x < 0.000000001 && Primal[it][j].y > 0.000000001) { //UPPER LEFT
                K.rho = 0.5065;
                K.u = 0.8939;
                K.v = 0.;
                K.p = 0.35;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;
            }
            if (Primal[it][j].x < 0.000000001 && Primal[it][j].y < 0.000000001) { //LOWER LEFT
                K.rho = 1.1;
                K.u = 0.8939;
                K.v = 0.8939;
                K.p = 1.1;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;
            }

            if (Primal[it][j].x > 0.000000001 && Primal[it][j].y < 0.000000001) { //LOWER RIGHT
                K.rho = 0.5065;
                K.u = 0.8939;
                K.v = 0.8939;
                K.p = 0.35;
                K.E = K.p / 0.4 + 0.5 * (pow(K.u, 2) + pow(K.v, 2));
                Primal[it][j].prev = K;
            }

        }
    }
    }
std::cout<<"initilisation terminee!"<<std::endl;

}

void Mesh::save(int a) {
      std::ofstream outdata;
        // opens the file
    //outvalue.open("value.txt");
    if (a==0){
        outdata.open("uinit.txt");

    }
    else if(a==1) {
        outdata.open("maj.txt");
    }
            for (unsigned int it = 0; it < Primal[1].size(); it++) {
                for (unsigned int j =0; j < Primal.size(); j++) {

                        outdata << Primal[it][j].x << "  " << Primal[it][j].y << "  " <<Primal[it][j].prev.rho<< \
                        "  "<<Primal[it][j].prev.u/Primal[it][j].prev.rho<<"  " <<Primal[it][j].prev.p<<"  " <<Primal[it][j].prev.E<< std::endl;

                }

            }

       outdata.close();
}

//---------------------------------------------------------------------
