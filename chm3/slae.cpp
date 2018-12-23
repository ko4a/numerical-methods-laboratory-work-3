#pragma once
#include "slae.h"

vector <real> DiagonalElements; // diagonal elements of the matrix A
vector <real> VectorF; // right side vector
vector <real> X;// first approximation and solution
vector <real> Values; // contains elements of matrix A line by line (bottom triangle only without main diagonal)
vector <int> PointerNumberOfValuesInRow; // contains pointers number of values in a row (see row sparse format)
vector <int> PointerColumnOfValueElement;// [i] element shows column of Values[i] in matrix A (see sparse format)
size_t maxitter, N, epsilon; //maxitter is number of maximum iterations  N is matrix size, epsilon is residual

#pragma region methods

// it works in compact format, and returns matrix*vectorForMultiply
// count of columns in matrix should equals count of rows in vectorForMultiply
// matrix should be square
vector <real> multiplyVectorByMatrix(vector<vector<real>>matrix, vector<real> vectorForMultiply)
{
	vector<real>result(vectorForMultiply.size(), 0.0);

	for (size_t i = 0; i < vectorForMultiply.size(); i++)
		for (size_t j = 0; j < vectorForMultiply.size(); j++)
			result[i] += matrix[i][j] * vectorForMultiply[j];
	return result;
}


//writes your vector to Console
void writeToConsole(vector<real> vectorToWrite)
{
	cout.precision(17);
	for (auto element : vectorToWrite)
		cout << element << ' ';
	cout << endl;
}

//writes your matrix to Console
void writeToConsole(vector<vector<real>> matrixToWrite)
{
	for (auto v : matrixToWrite)
	{
		for (auto x : v) cout << x << ' ';
		cout << endl;
	}
}




// returns the scalar product of (a,b)
// vectors A and B should have same sizes
real CalcScalarProduct(vector<real>A, vector<real> B)
{
	if (A.size() != B.size()) throw invalid_argument("A and B should have same sizes");

	real result = 0;
	for (size_t i = 0; i < A.size(); i++)
		result += A[i] * B[i];

	return result;

}

//calcs AlphaK  by p and r(argument p is p(k-1) and  r is r(k-1))
// AlphaK = (p(k-1),r(k-1))/(p(k-1),p(k-1))
// vectors p and r should have same sizes
real CalcAlphaK(vector<real>p, vector<real> r)
{
	return (CalcScalarProduct(p, r) / CalcScalarProduct(p, p));
}


//reads vector from file
template <typename T>
void ReadVectorFromFile(ifstream& file, vector<T>& myVector)
{
	for (auto &element : myVector)
		file >> element;
}

//returns result of  vector A - vector B
// vectors A and B should have same sizes
vector<real> VectorMinusVector(vector<real>A, vector<real>B)
{
	if (A.size() != B.size()) throw invalid_argument("A  and B should have same sizes");

	vector<real> result(A.size(), 0.0);
	for (size_t i = 0; i < A.size(); i++)
		result[i] = A[i] - B[i];
	return result;

}


// it works in row sparse format, and returns matrix*vectorForMultiply
// values -- should contain values, [i] element of columns should show column of value[i], diagonal should contain diagonal elements
// PointerOfNumberOfValuesInRow is ig array
// for more info see row sparse format
vector<real> multiplyVectorByMatrix(vector<real>Values, vector<int> ColumnsPointer, vector<int>PointerOfNumberOfValuesInRow, vector<real>diagonal, vector<real> multiplicator)
{
	vector <real> result(diagonal.size(), 0.0);
	int k1, k2;
	for (size_t i = 0; i < diagonal.size(); i++)
	{
		result[i] = diagonal[i] * multiplicator[i];
		k1 = PointerOfNumberOfValuesInRow[i];
		k2 = PointerOfNumberOfValuesInRow[i + 1];
		for (int k = k1; k < k2; k++)
		{
			result[i] += Values[k] * multiplicator[ColumnsPointer[k]];
			result[ColumnsPointer[k]] += Values[k] * multiplicator[i];
		}
	}
	return result;

}

// this method reads data from files, allocates memory for vectors,matrix, sets first approximation
void InitSolution()
{
	ifstream kuslauFile("kuslau.txt");
	ifstream igFile("ig.txt");
	ifstream jgFile("jg.txt");
	ifstream ggFile("gg.txt");
	ifstream diFile("di.txt");
	ifstream prFile("pr.txt");

	kuslauFile >> N;
	kuslauFile >> maxitter;
	kuslauFile >> epsilon;

	X.resize(N, 0.0);
	VectorF.resize(N, 0.0);
	DiagonalElements.resize(N, 0.0);
	PointerNumberOfValuesInRow.resize(N + 1, 0.0);


	ReadVectorFromFile(prFile, VectorF);
	ReadVectorFromFile(diFile, DiagonalElements);
	ReadVectorFromFile(igFile, PointerNumberOfValuesInRow);
	Values.resize(PointerNumberOfValuesInRow[N], 0.0);
	PointerColumnOfValueElement.resize((size_t)PointerNumberOfValuesInRow[N], 0.0);
	ReadVectorFromFile(ggFile, Values);
	ReadVectorFromFile(jgFile, PointerColumnOfValueElement);



}


real VectorNorm(vector<real> myVector)
{
	real result=0;
	for (auto v : myVector)
		result += v * v;
	return (sqrt(result));
	
}

// calcs norm(A)/norm(B) 
real CalcRealtiveResidual(vector<real> A, vector<real> B)
{
	return (VectorNorm(A) / VectorNorm(B));
}


// solves SLAE Ax=f 
//returns X
vector<real> solveSLAE(vector<real>gg, vector<real>f, vector<real>x, vector<real> diagonal, vector<int> ig, vector<int> jg)
{
	vector <real> r = VectorMinusVector(f, multiplyVectorByMatrix(gg, jg, ig, diagonal, x));
	vector<real> z = r;
	vector<real> p = multiplyVectorByMatrix(gg, jg, ig, diagonal, x);
	real alphaK,bettaK;
	
	vector <real> result(diagonal.size(),0.0);
	int itter = 0;
	while (itter<maxitter && CalcRealtiveResidual(r,f)>epsilon)
	{
		alphaK = CalcAlphaK(p, r);
		result = VectorPlusVector(result, VectorMultiplyByConst(z, alphaK));
		r = VectorMinusVector(r, VectorMultiplyByConst(p, alphaK));
		bettaK = CalcBettaK(p, r, gg, ig, jg, diagonal);
		z = VectorPlusVector(r, VectorMultiplyByConst(z, bettaK));
		p = VectorPlusVector(multiplyVectorByMatrix(gg, jg, ig, diagonal, r), VectorMultiplyByConst(p, bettaK));
	}

	return result;
}

//calcs BettaK  by p, A and r(argument p is p(k-1) and  r is r(k-1), A is values,ig,jg,diagonal - row sparse format )
// vectors p and r should have same sizes
real CalcBettaK(vector<real>p, vector<real>r, vector<real>values, vector<int> ig, vector<int> jg, vector<real>diagonal)
{
	vector <real> tmp = multiplyVectorByMatrix(values, jg, ig, diagonal, r);
	return -(CalcScalarProduct(p, tmp) / (CalcScalarProduct(p, p)));
}

//returns result of  vector A + vector B
// vectors A and B should have same sizes
vector<real> VectorPlusVector(vector<real> A, vector<real> B)
{
	if (A.size() != B.size()) throw invalid_argument("A and B should have same sizes");

	vector <real> result(A.size(), 0.0);
	for (size_t i = 0; i < A.size(); i++)
		result[i] = A[i] + B[i];
	return result;
}

//multiplies a vector by a constant
vector<real> VectorMultiplyByConst(vector<real> multiplyVector, real constant)
{
	for (auto& v : multiplyVector)
		v *= constant;
	return multiplyVector;
}

// returns true if A == B, false if A!=B
//this function needs only for debug
bool areVectorsEqual(vector<real> A, vector<real>B)
{
	for (size_t i = 0; i < A.size(); i++)
		if (A[i] != B[i])return false;
	return true;
}



// creates Hilbert Matrix n size;
void CreateHilbertMatrix(int n)
{
	PointerNumberOfValuesInRow.resize(n + 1);
	PointerNumberOfValuesInRow[0] = 0;
	for (int i = 1; i < n + 1; i++)
		PointerNumberOfValuesInRow[i] = PointerNumberOfValuesInRow[i - 1] + i - 1;
	DiagonalElements.resize(n);
	Values.resize(PointerNumberOfValuesInRow[n]);
	PointerColumnOfValueElement.resize(PointerNumberOfValuesInRow[n]);
	int k1, k2;
	for (int i = 0; i < n; i++)
	{
		DiagonalElements[i] = double(1.0 / (2.0*i + 1.0));
		k1 = PointerNumberOfValuesInRow[i];
		k2 = PointerNumberOfValuesInRow[i + 1];
		int l = 0;
		for (int k = k1; k < k2; k++, l++)
		{
			DiagonalElements[k] = double(1 / (i + l + 1.0));
			PointerColumnOfValueElement[k] = l;
		}
	}
}

#pragma endregion
