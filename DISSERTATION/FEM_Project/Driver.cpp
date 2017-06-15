#include <iostream>
#include <cmath>

#include "FE_Solution.hpp"
#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Quadrature.hpp"

int main(int argc, char* argv[])
{
    Quadrature* Quad = new Quadrature(0,2,4);

    Vector* functionPoints = new Vector(Quad->GetGQPoints());

    for (int i=0; i<functionPoints->GetSize(); i++)
    {
        (*functionPoints)[i] = pow((*functionPoints)[i],2.0);
    }

    std::cout << Quad->GaussQuadrature(*functionPoints) << std::endl;

    delete Quad;
    return 0;
}
