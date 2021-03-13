#ifndef CELL_H
#define CELL_H


struct Etat  // ATTENTION ON RESOUT EN VAR CONSERVATIVES
{

    friend Etat operator*(double scal, Etat A);

    friend Etat operator+( Etat A, Etat B);

    friend Etat operator-( Etat A, Etat B);

    double rho;
    double u;
    double v;
    double p;
    double E;

};


struct cell_p {
    int bord;
    double x;
    double y;
    Etat prev;
    Etat next;
    Etat FL, FR, GU, GD;

};

struct cell_bord {
    int ind;
    double x;
    double y;
    Etat prev;
    Etat next;
};

struct cell_d {
    int ind;
    double x;
    double y;
    Etat prev;
    Etat next;

};

struct cell_d_bord {
    int ind;
    double x;
    double y;
    Etat prev;
    Etat next;
    double mes;

};




#endif
