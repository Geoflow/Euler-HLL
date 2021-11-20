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




Mesh::Mesh(int i, int j, double L, double H) {

    m_nx = i;
    m_ny = j;
    m_Long = L;
    m_Haut = H;
    m_dx =  m_Long / m_nx ;
    m_dy =  m_Haut / m_ny ;




    Primal.resize(m_nx, std::vector<cell_p>(m_ny));

    //--------------PRIMAL MESH-----------------------------------

    for (int i = 1; i < m_nx-1 ; ++i) {
        for (int j = 1; j < m_ny-1 ; ++j) {

            //INTERNAL CELLS--
            Primal[i][j].bord = 0;
            Primal[i][j].x = (i+0.5) * m_dx;//-m_Long +
            Primal[i][j].y = (j+0.5) * m_dy;//-m_Haut +


        }
    }
   for (int i = 0; i < m_nx; ++i) {
        //BOTTOM BORDER
        Primal[i][0].bord = 1;
        Primal[i][0].x =  (i+0.5) * m_dx;
        Primal[i][0].y = 0.5 *m_dy;
        //TOP BORDER
        Primal[i][m_nx - 1].bord = 3;
        Primal[i][m_nx - 1].x = (i+0.5) * m_dx;
        Primal[i][m_nx - 1].y = m_Haut- 0.5 *m_dy;

    }

    for (int j = 0; j < m_ny; ++j) {
        //LEFT BORDER
        Primal[0][j].bord = 4;
        Primal[0][j].x =  0.5 * m_dx;
        Primal[0][j].y = (j+0.5) * m_dy;
        //RIGHT BORDER
        Primal[m_nx - 1][j].bord = 2;
        Primal[m_nx - 1][j].x = m_Long - 0.5*m_dx;
        Primal[m_nx - 1][j].y =  (j+0.5)*m_dy;

    }

    Primal[0][0].bord = 5;
    Primal[m_nx - 1][0].bord = 5;
    Primal[0][m_ny-1].bord = 5;
    Primal[m_nx - 1][m_ny-1].bord = 5;


    std::cout << "Primal Mesh Size " << m_nx*m_ny << std::endl;
    /*for (unsigned int it = 0; it < Primal[1].size(); it++) {
        for (unsigned int j = 0; j < Primal.size(); j++) {
if(Primal[it][j].bord==5) std::cout << Primal[it][j].x << " " << Primal[it][j].y << "  " << Primal[it][j].bord << std::endl;
        }
    }*/
}


double Mesh::getdx() { return m_dx; }
double Mesh::getdy() { return m_dy; }

void Mesh::uinit(int choice) {


//-------------------------------------------SOD SHOCK -----------------------------------------------------------------

    if (choice == 0) {
        for (unsigned int it = 0; it < Primal.size(); it++) {
            for (unsigned int j = 0; j < Primal.size(); j++) {

                if (Primal[it][j].x <= 0.5) {// UPPER RIGHT
                    Primal[it][j].prev.rho = 1.;
                    Primal[it][j].prev.u = 0.;
                    Primal[it][j].prev.v = 0.;
                    Primal[it][j].prev.p = 1.;

                } else { //UPPER LEFT
                    Primal[it][j].prev.rho = 0.125;
                    Primal[it][j].prev.u = 0.;
                    Primal[it][j].prev.v = 0.;
                    Primal[it][j].prev.p = 0.1;


                }
            }
        }
    }
//-------------------------------------TEST CASE NUMERO 3-----------------------------------------------------------
        if (choice == 3) {
            for (unsigned int it = 0; it < Primal.size(); it++) {
                for (unsigned int j = 0; j < Primal.size(); j++) {

                    if (Primal[it][j].x >= 0.4 && Primal[it][j].y >= 0.4 ) {// UPPER RIGHT
                        Primal[it][j].prev.rho = 1.5;
                        Primal[it][j].prev.u = 0.;
                        Primal[it][j].prev.v = 0.;
                        Primal[it][j].prev.p = 1.5;

                    }
                    if (Primal[it][j].x < 0.4 && Primal[it][j].y > 0.4 ) { //UPPER LEFT
                        Primal[it][j].prev.rho = 0.532258064516129;
                        Primal[it][j].prev.u = 1.206045378311055;
                        Primal[it][j].prev.v = 0.;
                        Primal[it][j].prev.p = 0.3;


                    }
                    if (Primal[it][j].x <= 0.4 && Primal[it][j].y <= 0.4 ) { //LOWER LEFT
                        Primal[it][j].prev.rho = 0.137992831541219 ;
                        Primal[it][j].prev.u = 1.206045378311055;
                        Primal[it][j].prev.v = 1.206045378311055;
                        Primal[it][j].prev.p = 0.029032258064516;


                    }

                    if (Primal[it][j].x > 0.4 && Primal[it][j].y < 0.4 ) { //LOWER RIGHT
                        Primal[it][j].prev.rho = 0.532258064516129;
                        Primal[it][j].prev.u = 0.;
                        Primal[it][j].prev.v = 1.206045378311055;
                        Primal[it][j].prev.p = 0.;


                    }

                }
            }
        }

        if (choice == 4) {
            //----------------------- TEST CASE NUMERO 4------------------------------------------------------------------------
            for (unsigned int it = 0; it < Primal[1].size(); it++) {
                for (unsigned int j = 0; j < Primal.size(); j++) {

                    if (Primal[it][j].x > 0.4 && Primal[it][j].y > 0.4 ) {// UPPER RIGHT
                        Primal[it][j].prev.rho = 1.1;
                        Primal[it][j].prev.u = 0.;
                        Primal[it][j].prev.v = 0.;
                        Primal[it][j].prev.p = 1.1;

                    }
                    if (Primal[it][j].x < 0.4 && Primal[it][j].y > 0.4 ) { //UPPER LEFT
                        Primal[it][j].prev.rho = 0.5065;
                        Primal[it][j].prev.u = 0.8939;
                        Primal[it][j].prev.v = 0.;
                        Primal[it][j].prev.p = 0.35;

                    }
                    if (Primal[it][j].x < 0.4 && Primal[it][j].y < 0.4 ) { //LOWER LEFT
                        Primal[it][j].prev.rho = 1.1;
                        Primal[it][j].prev.u = 0.8939;
                        Primal[it][j].prev.v = 0.8939;
                        Primal[it][j].prev.p = 1.1;

                    }

                    if (Primal[it][j].x > 0.4 && Primal[it][j].y < 0.4 ){ //LOWER RIGHT
                        Primal[it][j].prev.rho = 0.5065;
                        Primal[it][j].prev.u = 0.;
                        Primal[it][j].prev.v = 0.8939;
                        Primal[it][j].prev.p = 0.35;


                    }

                }
            }
        }

    if (choice == 6) {
        //----------------------- TEST CASE NUMERO 6------------------------------------------------------------------------
        for (unsigned int it = 0; it < Primal[1].size(); it++) {
            for (unsigned int j = 0; j < Primal.size(); j++) {

                if (Primal[it][j].x > 0.5 && Primal[it][j].y > 0.5 ) {// UPPER RIGHT
                    Primal[it][j].prev.rho = 1.;
                    Primal[it][j].prev.u = .75;
                    Primal[it][j].prev.v = -.5;
                    Primal[it][j].prev.p = 1.;

                }
                if (Primal[it][j].x < 0.5 && Primal[it][j].y > 0.5 ) { //UPPER LEFT
                    Primal[it][j].prev.rho = 2.;
                    Primal[it][j].prev.u = 0.75;
                    Primal[it][j].prev.v = 0.5;
                    Primal[it][j].prev.p = 1.;

                }
                if (Primal[it][j].x < 0.5 && Primal[it][j].y < 0.5 ) { //LOWER LEFT
                    Primal[it][j].prev.rho = 1.;
                    Primal[it][j].prev.u = -.75;
                    Primal[it][j].prev.v = .5;
                    Primal[it][j].prev.p = 1.;

                }

                if (Primal[it][j].x > 0.5 && Primal[it][j].y < 0.5 ){ //LOWER RIGHT
                    Primal[it][j].prev.rho = 3.;
                    Primal[it][j].prev.u = 1.;
                    Primal[it][j].prev.v = -.75;
                    Primal[it][j].prev.p = -.5;


                }

            }
        }
    }

 //---------------------------------------------------------------------------------------------------------------------
 //---------------------------VARIABLES CONSERVATIVES-------------------------------------------------------------------
    for (unsigned int it = 0; it < Primal[1].size(); it++) {
        for (unsigned int j = 0; j < Primal.size(); j++) {
            Primal[it][j].prev.E = Primal[it][j].prev.p / 0.4 +0.5*Primal[it][j].prev.rho * (pow(Primal[it][j].prev.u, 2) + pow(Primal[it][j].prev.v, 2));//OK DONT TOUCH
            Primal[it][j].prev.u *= Primal[it][j].prev.rho ;
            Primal[it][j].prev.v *= Primal[it][j].prev.rho ;
          Primal[it][j].prev.E*= Primal[it][j].prev.rho ;


        }
        }
        std::cout << "initilisation terminee!" << std::endl;
}

void Mesh::save(int a, std::string filename) {

      std::ofstream outdata;
      std::string nom_init;
      nom_init="init_"+filename+".txt";
      std::string nom_maj;
      nom_maj="maj_"+filename+".txt";

    switch (a)
    {
    case 0:
        outdata.open(nom_init);
            for (unsigned int it = 0; it < Primal.size(); it++) {
                for (unsigned int j =0; j < Primal.size(); j++) {

                    outdata << Primal[it][j].x << " " << Primal[it][j].y << " " <<Primal[it][j].prev.rho<< \
                        "  "<<Primal[it][j].prev.u/Primal[it][j].prev.rho<<" " <<Primal[it][j].prev.p<<" " <<Primal[it][j].prev.E<< std::endl;

                }

            }
            outdata.close();
    case 1:
        outdata.open(nom_maj);
            for (unsigned int it = 0; it < Primal.size(); it++) {
                for (unsigned int j =0; j < Primal.size(); j++) {

                    outdata << Primal[it][j].x << " " << Primal[it][j].y << " " <<Primal[it][j].prev.rho<< \
                        " "<<Primal[it][j].prev.u/Primal[it][j].prev.rho<<" " <<Primal[it][j].prev.p<<" " <<Primal[it][j].prev.E/Primal[it][j].prev.rho<< std::endl;

                }

            }
            outdata.close();
    }
}

