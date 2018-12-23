#include "slae.h"

int main()
{
	InitSolution();
	X = solveSLAE(Values, VectorF, X, DiagonalElements, PointerNumberOfValuesInRow, PointerColumnOfValueElement);
	writeToConsole(X);
	system("pause");
}