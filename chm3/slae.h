#pragma region head
#pragma once
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
typedef double real;



extern vector <real> DiagonalElements; // diagonal elements of the matrix A
extern vector <real> VectorF; // right side vector
extern vector <real> X;// first approximation and solution
extern vector <real> Values; // contains elements of matrix A line by line (bottom triangle only without main diagonal)
extern vector <int> PointerNumberOfValuesInRow; // contains pointers number of values in a row (see row sparse format)
extern vector <int> PointerColumnOfValueElement;// [i] element shows column of Values[i] in matrix A (see sparse format)
extern size_t maxitter, N, epsilon; //maxitter is number of maximum iterations  N is matrix size, epsilon is residual



#pragma endregion
#pragma region methods

// it works in compact format, and returns matrix*vectorForMultiply
// count of columns in matrix should equals count of rows in vectorForMultiply
// matrix should be square
vector <real> multiplyVectorByMatrix(vector<vector<real>>, vector<real>);

//writes your vector to Console
void writeToConsole(vector<real>);

//writes your matrix to Console
void writeToConsole(vector<vector<real>>);


// returns the scalar product of (a,b)
// vectors A and B should have same sizes
real CalcScalarProduct(vector<real>, vector<real>);

//calcs AlphaK  by p and r(argument p is p(k-1) and  r is r(k-1))
// AlphaK = (p(k-1),r(k-1))/(p(k-1),p(k-1))
// vectors p and r should have same sizes
real CalcAlphaK(vector<real>, vector<real> );

//reads vector from file
template <typename T>
void ReadVectorFromFile(ifstream&, vector<T>& );

//returns result of  vector A - vector B
// vectors A and B should have same sizes
vector<real> VectorMinusVector(vector<real>, vector<real>);

// it works in row sparse format, and returns matrix*vectorForMultiply
// values -- should contain values, [i] element of columns should show column of value[i], diagonal should contain diagonal elements
// PointerOfNumberOfValuesInRow is ig array
// for more info see row sparse format
vector<real> multiplyVectorByMatrix(vector<real>, vector<int> , vector<int>, vector<real>, vector<real> );

// this method reads data from files, allocates memory for vectors,matrix, sets first approximation
void InitSolution();

// solves SLAE Ax=f 
//returns X
vector<real> solveSLAE(vector<real>, vector<real>, vector<real>, vector<real>, vector<int>, vector<int> );

//calcs BettaK  by p, A and r(argument p is p(k-1) and  r is r(k-1), A is values,ig,jg,diagonal - row sparse format )
// vectors p and r should have same sizes
real CalcBettaK(vector<real>, vector<real>, vector<real>, vector<int> , vector<int> , vector<real>);

//returns result of  vector A + vector B
// vectors A and B should have same sizes
vector<real> VectorPlusVector(vector<real> , vector<real> );

//multiplies a vector by a constant
vector<real> VectorMultiplyByConst(vector<real>, real);

// returns true if A == B, false if A!=B
//this function needs only for debug
bool areVectorsEqual(vector<real> , vector<real>);
#pragma endregion





