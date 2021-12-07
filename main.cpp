
#include "mesh.h"
#include "euler.h"

int main() {

    Mesh Maillage(100,100, 1., 1.);
    Maillage.uinit(3);
    Euler PB(Maillage);
//PB.LF_Solver(0.16, 3); //OK
//PB.UPWIND_Solver(0.16, 3);//OK
PB.HLL_SPLIT(0.08, 3); //not OK

//PB.HLL_Solver(0.002, 3);

    return 0;
}

