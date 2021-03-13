#include <cstdlib>
#include <iostream>
#include "mesh.h"
#include "euler.h"
int main() {


    Mesh Maillage(200, 200, .5, .5);
    Maillage.uinit(4);
    Maillage.save(0);
    Euler PB(Maillage);
    PB.Update(0.25);

    return 0;
}
