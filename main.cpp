
#include "mesh.h"
#include "euler.h"

int main() {

    Mesh Maillage(128,128, 1., 1.);
    Maillage.uinit(4);
    Euler PB(Maillage);
//PB.LF_Solver(0.3, 3);
//PB.HLL_SPLIT(0.32, 4);
PB.HLL_Solver(0.001, 4);
//PB.UPWIND_Solver(0.05, 0);
    return 0;
}

