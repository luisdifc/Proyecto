#ifndef SVR_H
#define SVR_H

#include <cmath>
#include <iostream>
#include <vector>

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
    Svr(float C, float tol, int kernelType, int useLinearOptim, int size);
    virtual ~Svr();
    float Output(int i);
    int TakeStep(int i1, int i2);
    float GetError(int i1);

    float Kernel(vector<float> v1, vector<float> v2);
    float ComputeB (float E1, float a1, float a1New, float a2New, float k11, float k12, float k22, float y1, float y2, float a2, float E2);
    float SecondHeuristic(vector<int> nonBoundIndices, float E2);
    float Error(int i2);


    float DotProduct(vector<float> v1, vector<float> v2);
    int FillWithCeros(int size, vector<float> &vector);
    float GetMax(float n1, float n2);
    float GetMin(float n1, float n2);
};

#endif // SVR_H
