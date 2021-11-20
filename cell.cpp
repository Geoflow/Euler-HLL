#include "cell.h"

Etat operator*(double scal,  Etat A) {

    A.rho *= scal;
    A.u *= scal;
    A.v *= scal;
    A.E *= scal;
    return A;
}

Etat operator+( Etat A,  const Etat& B) {
 A.rho += B.rho;
 A.u += B.u;
 A.v += B.v;
 //A.p += B.p;
   A.E +=B.E;
    return A;
}

Etat operator-(Etat A,  const Etat& B){

     A.rho -= B.rho;
    A.u -= B.u;
    A.v -= B.v;
    //A.p-= B.p;
    A.E -= B.E;
    return A;
}

Etat operator/( Etat A, double scal) {

    A.rho /= scal;
    A.u /= scal;
    A.v /= scal;
    A.E /= scal;
    return A;
}
