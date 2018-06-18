#ifndef SVR_H
#define SVR_H

#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib> 
#include <ctime> 
#include<fstream>
#include<string>

using namespace std;

class Svr {

//Atributos
public:
vector< vector<float> > X;
vector<float> Y;
vector< vector<float> > m;
vector< vector<float> > n;
vector<float> alphas;
vector<float> errors;
vector<float> w;

float b;
float eps; //epsilon
float C;
float tol;

int kernelType; // 1 = linear
int useLinearOptim;
	

//Metodos
public:
	//For the moment X and Y wont be parametrized for testing, but later they will
    Svr(float C, float tol, int kernelType, int useLinearOptim, vector<vector<float>> XData, vector<float> YData);
    virtual ~Svr();
    float Output(int i);
    int TakeStep(int i1, int i2, float a2, float y2, float E2, vector<float> x2);
    float GetError(int i1);

    float Kernel(vector<float> v1, vector<float> v2);
    float ComputeB (float E1, float a1, float a1New, float a2New, float k11, float k12, float k22, float y1, float y2, float a2, float E2);
    void ComputeW(vector<float> multipliers, vector< vector<float> > X, vector<float> y);
    float FirstHeuristic();
    float SecondHeuristic(vector<int> nonBoundIndices, float E2);
    int ExamineExample (int i2);
    float Error(int i2);
    vector<int> GetNonBoundIndexes();
    void MainRoutine();
    float Predict(vector<float> newX);


    float DotProduct(vector<float> v1, vector<float> v2);
    vector<float> VectorByScalar(vector<float> v1, float scalar);
    vector<float> VectorSum(vector<float> v1, vector<float> v2);
    int FillWithCeros(int size, vector<float> &vector);
    float GetMax(float n1, float n2);
    float GetMin(float n1, float n2);  
    int RandNumGenerator(int n1, int n2);

    void PrintVector(vector<float> vector);
    void PrintMatrix(vector< vector<float> > array);
    void GiveSizeMatrix(int i, int j, vector< vector<float> > &array);
};

#endif // SVR_H