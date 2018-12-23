#pragma region head
#pragma once
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
typedef double real;



extern vector <real> DiagonalElements; 
extern vector <real> VectorF; 
extern vector <real> X;
extern vector <real> Values; // contains elements of matrix A line by line (bottom triangle only without main diagonal)
extern vector <int> PointerNumberOfValuesInRow; // contains pointers number of values in a row (see row sparse format)
extern vector <int> PointerColumnOfValueElement;// [i] element shows column of Values[i] in matrix A (see sparse format)
extern size_t maxitter, N, epsilon;



#pragma endregion
#pragma region methods

// it works in compact format
vector <real> multiplyVectorByMatrix(vector<vector<real>>, vector<real>);

void writeToConsole(vector<real>);


void writeToConsole(vector<vector<real>>);



real CalcScalarProduct(vector<real>, vector<real>);


// AlphaK = (p(k-1),r(k-1))/(p(k-1),p(k-1))
real CalcAlphaK(vector<real>, vector<real> );


template <typename T>
void ReadVectorFromFile(ifstream&, vector<T>& );


vector<real> VectorMinusVector(vector<real>, vector<real>);

//works in row sparse format
vector<real> multiplyVectorByMatrix(vector<real>, vector<int> , vector<int>, vector<real>, vector<real> );

//  reads data from files, allocates memory for vectors,matrix, sets first approximation
void InitSolution();

// solves SLAE Ax=f 

vector<real> solveSLAE(vector<real>, vector<real>, vector<real>, vector<real>, vector<int>, vector<int> );

//calcs BettaK  by p, A and r(argument p is p(k-1) and  r is r(k-1), A is values,ig,jg,diagonal - row sparse format )
real CalcBettaK(vector<real>, vector<real>, vector<real>, vector<int> , vector<int> , vector<real>);


vector<real> VectorPlusVector(vector<real> , vector<real> );


vector<real> VectorMultiplyByConst(vector<real>, real);


bool areVectorsEqual(vector<real> , vector<real>);
#pragma endregion





