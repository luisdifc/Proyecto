#include "Svr.h"


//Constructor
Svr::Svr(float C, float tol, int kernelType, int useLinearOptim, int size) 
{
    cout << "\nConstruyendo un SVR..!\n" << endl;

	this->FillWithCeros(size, this->alphas);
	this->FillWithCeros(size, this->errors);
	this->FillWithCeros(size, this->w);

    this->b = 0;
    this->eps = 0.001;
    this->C = C;

    this->kernelType = kernelType;
    this->useLinearOptim = useLinearOptim;
}


//Destructor
Svr::~Svr() 
{

}	


float Svr::Output(int i) 
{
	int error = -1;
	float sum = 0;

	if (useLinearOptim == 1) 
	{
		return ( this->DotProduct(this->w, this->X[i]) - this->b ) ;
	} 
	else 
	{
		for (int index = 0; index < this->m.size(); index++) 
		{
			sum += (this->alphas[index] * this->Y[index] * this->Kernel(this->X[index], this->X[i]));
		}
		return sum - this->b;
	}

	return error;
}


int Svr::TakeStep(int i1, int i2) 
{
	//return = 0 false, return 1 = true, return errror = error
	int error = -1;

	float a1 = 0.0;
	float y1 = 0.0; 
	vector<float> x1;
	float E1 = 0.0;
	float s = 0.0;

	if (i1 == i2) 
	{
		return 0;
	}

	a1 = this->alphas[i1];
	y1 = this->Y[i1];
	x1 = this->X[i1];
	E1 = this->GetError(i1);

	s = y1 * this->Y[i2];

}


float Svr::GetError(int i1) 
{
	int error = -1;

	if (0 < this->alphas[i1] && this->alphas[i1] < this->C) 
	{
		return this->errors[i1];
	}
	else 
	{
		return this->Output(i1) - this->Y[i1];
	}

	return (float)error;
}


float Svr::Kernel(vector<float> v1, vector<float> v2) 
{
	int error = -1;

	if (this->kernelType == 1) 
	{
		return this->DotProduct(v1, v2); //linear kernel
	}

	return error;
}



//Dot product between 2 vectors
float Svr::DotProduct (std::vector<float> v1, std::vector<float> v2) 
{
	float result = 0;

	for (int index = 0; index < v1.size(); index++) 
	{
		result += v1[index] * v2[index];
	}

	return result;
}

int Svr::FillWithCeros(int size, vector<float> &vector) 
{
	for (int index = 0; index < size; index++) 
	{
		vector.push_back(0);
	}

	return 0;
}

