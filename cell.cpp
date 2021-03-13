#include "cell.h"

Etat operator*(double scal, Etat A) {
    Etat C;
    C.rho = scal*A.rho;
    C.u = scal*A.u;
    C.v = scal*A.v;
    C.E = scal*A.E;
    return C;
}

Etat operator+( Etat A,  Etat B) {
   Etat C;
    C.rho= A.rho + B.rho;
    C.u=A.u + B.u;
    C.v=A.v + B.v;
   // C.p=A.p + B.p;
    C.E=A.E + B.E;
    return C;
}

Etat operator-(Etat A, Etat B) {
    Etat C;
    C.rho= A.rho - B.rho;
    C.u=A.u - B.u;
    C.v=A.v - B.v;
    //C.p=A.p - B.p;
    C.E=A.E - B.E;
    return C;
}