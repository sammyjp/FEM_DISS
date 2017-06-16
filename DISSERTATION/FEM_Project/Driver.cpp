#include <iostream>
#include <cmath>

#include "Element.hpp"
#include "FE_Solution.hpp"
#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Quadrature.hpp"

int main(int argc, char* argv[])
{
    Element* E = new Element(1);

    Vector* nodes = new Vector(2);
    (*nodes)[0] = 0;
    (*nodes)[1] = 1;

    Vector* localCoords = new Vector(3);
    (*localCoords)[0] = -1;
    (*localCoords)[1] = 0;
    (*localCoords)[2] = 1;

    Vector* globalCoords = new Vector(localCoords->GetSize());

    E->MapLocalToGlobal(*nodes, *localCoords, *globalCoords);

    std::cout << *globalCoords << std::endl;

    *localCoords = 0;
    E->MapGlobalToLocal(*nodes, *globalCoords, *localCoords);
    std::cout << *localCoords << std::endl;

    Element* EMat = new Element(2);

    Matrix* nodeMatrix = new Matrix(2,3);
    (*nodeMatrix)(1,1) = 0.5;
    (*nodeMatrix)(1,2) = 1;
    (*nodeMatrix)(1,3) = 0;
    (*nodeMatrix)(2,1) = 0;
    (*nodeMatrix)(2,2) = 0;
    (*nodeMatrix)(2,3) = 1;

    Matrix* localMatrix = new Matrix(2,3);
    (*localMatrix)(1,1) = 0.5;
    (*localMatrix)(1,2) = 1;
    (*localMatrix)(1,3) = 0;
    (*localMatrix)(2,1) = 0;
    (*localMatrix)(2,2) = 0;
    (*localMatrix)(2,3) = 1;

    Matrix* globalMatrix = new Matrix(localMatrix->GetNumberOfRows(), localMatrix->GetNumberOfColumns());

    EMat->MapLocalToGlobal(*nodeMatrix, *localMatrix, *globalMatrix);

    std::cout << *globalMatrix << std::endl;

    for (int i=1; i<=localMatrix->GetNumberOfRows(); i++)
    {
        for (int j=1; j<=localMatrix->GetNumberOfColumns(); j++)
        {
            (*localMatrix)(i,j) = 0;
        }
    }

    EMat->MapGlobalToLocal(*nodeMatrix, *globalMatrix, *localMatrix);
    std::cout << *localMatrix << std::endl;


    delete E;
    delete EMat;
    delete localCoords;
    delete globalCoords;
    delete localMatrix;
    delete globalMatrix;


    return 0;
}
