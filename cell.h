#ifndef CELL_H
#define CELL_H


struct Etat  // ATTENTION ON RESOUT EN VAR CONSERVATIVES
{

    friend Etat operator*(double scal, Etat A);
    friend Etat operator/(Etat A,double scal );
    friend Etat operator+(Etat A,  const Etat& B);
    friend Etat operator-( Etat A,  const Etat& B);

    double rho;
    double u;
    double v;
    double p;
    double E;

};

Etat operator*(double scal, Etat A);
Etat operator+(Etat A,  const Etat& B);
Etat operator-(Etat A,  const Etat& B);
Etat operator/(Etat A,double scal );


struct cell_p {
    int bord;
    double x;
    double y;
    Etat prev;
    Etat next;
    Etat FL, FR, GU, GD;

};





#endif
