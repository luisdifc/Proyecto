#ifndef SVR_H
#define SVR_H

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



    float DotProduct(vector<float> v1, vector<float> v2);
    int FillWithCeros(int size, vector<float> &vector);
};

#endif // SVR_H
