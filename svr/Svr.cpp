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
	//return = 0 false, return 1 = true, return error = error
	int error = -1;

	float a1 = 0.0;
	float y1 = 0.0; 
	vector<float> x1;
	float E1 = 0.0;
	float s = 0.0;
	float L = 0.0;
	float H = 0.0;
	float k11 = 0.0;
	float k12 = 0.0;
	float k22 = 0.0;
	float eta = 0.0;
	float a1New = 0.0;
	float a2New = 0.0;
	float f1 = 0.0;
	float f2 = 0.0;
	float L1 = 0.0;
	float H1 = 0.0;
	float Lobj = 0.0;
	float Hobj = 0.0;
	float newB = 0.0;
	float deltaB = 0.0;
	float delta1 = 0.0;
	float delta2 = 0.0;

	if (i1 == i2) 
	{
		return 0;
	}

	a1 = this->alphas[i1];
	y1 = this->Y[i1];
	x1 = this->X[i1];
	E1 = this->GetError(i1);

	s = y1 * this->Y[i2];

	if (y1 != this->Y[i2]) 
	{
		L = GetMax(0.0, this->alphas[i2] - a1);
		H = GetMin(this->C, this->C + this->alphas[i2] - a1);
	}
	else 
	{
		L = GetMax(0.0, this->alphas[i2] + a1 - this->C);
		H = GetMin(this->C, this->alphas[i2] + a1);
	}

	if (L == H) 
	{
		return 0;
	}

	k11 = this->Kernel(x1, x1);
	k12 = this->Kernel(x1, this->X[i2]);
	k22 = this->Kernel(this->X[i2], this->X[i2]);

	eta = k11 + k22 - (2 * k12);

	if (eta > 0 ) 
	{
		a2New = this->alphas[i2] + this->Y[i2] * (E1 - this->GetError(i2)) / eta;
		if (a2New < L) 
		{
			a2New = L;
		}
		else 
		{
			if (a2New > H) 
			{
				a2New = H;
			}
		}
	}
	else 
	{
		f1 = y1 * (E1 + this->b) - a1 * k11 - s * this->alphas[i2] * k12;
		f2 = this->Y[i2] * (this->GetError(i2) + this->b) - s * a1 * k12 - this->alphas[i2] * k22;
		L1 = a1 + s * (this->alphas[i2] - L);
		H1 = a1 + s * (this->alphas[i2] - H);
		Lobj = L1 * f1 + L * f2 + 0.5 * pow(L1, 2) * k11 + 0.5 * pow(L, 2) * k22 + s * L * L1 * k12;
		Hobj = H1 * f1 + H * f2 + 0.5 * pow(H1, 2) * k11 + 0.5 * pow(H, 2) * k22 + s * L * H1 * k12;
		if (Lobj < Hobj - this->eps) 
		{
			a2New = L;
		}
		else 
		{
			if (Lobj > Hobj + this->eps) 
			{
				a2New = H;
			}
			else 
			{
				a2New = this->alphas[i2];
			}
		}
	}

	if (abs(a2New - this->alphas[i2]) < this->eps * (a2New + this->alphas[i2] + this->eps)) 
	{
		return 0;
	}

	a1New = a1 + s * (this->alphas[i2] - a2New);
	newB = this->ComputeB(E1, a1, a1New, a2New, k11, k12, k22, y1, this->Y[i2], this->alphas[i2], this->GetError(i2));
	deltaB = newB = this->b;
	this->b = newB;

	if (this->useLinearOptim) 
	{
		// TO DO TO DO TO DO TO DO
	}

	delta1 = y1 * (a1New - a1);
	delta2 = this->Y[i2] * (a2New - this->alphas[i2]);

	for (int index = 0; index < this->m.size(); index++) 
	{
		if (0 < this->alphas[index] && this->alphas[index] < this->C) 
		{
			this->errors[index] += delta1 * this->Kernel(x1, this->X[index]) + delta2 * this->Kernel(this->X[i2], this->X[index]) - deltaB;
		}
	}

	this->errors[i1] = 0;
	this->errors[i2] = 0;

	this->alphas[i1] = a1New;
	this->alphas[i2] = a2New;

	return 1;
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


float Svr::ComputeB (float E1, float a1, float a1New, float a2New, float k11, float k12, float k22, float y1, float y2, float a2, float E2) 
{
	float b1 = 0.0;
	float b2 = 0.0;
	float newB = 0.0;

	b1 = E1 + y1 * (a1New - a1) * k11 + y2 * (a2New - a2) * k12 + this->b;
	b2 = E2 + y1 * (a1New - a1) * k12 + y2 * (a2New - a2) * k22 + this->b;

	if (0 < a1New && this->C > a1New) 
	{
		newB = b1;
	}
	else 
	{
		if (0 < a2New && this->C > a2New) 
		{
			newB = b2;
		} 
		else 
		{
			newB = (b1 + b2) / 2.0;
		}
	}

	return newB;
}


float Svr::SecondHeuristic(vector<int> nonBoundIndices, float E2) 
{
	float i1 = -1;
	float E1 = 0.0;
	float step = 0.0;
	float max = 0.0;

	if (nonBoundIndices.size() > 1) 
	{
		max = 0.0;
	}

	for (int index = 0; index < nonBoundIndices.size(); index++) 
	{
		E1 = this->errors[index] - this->Y[index];
		step = abs(E1 - E2);
		if (step > max) 
		{
			max = step;
			i1 = index;
		}
	}

	return i1;
}

float Svr::Error(int i2) 
{
	return this->Output(i2);
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


float Svr::GetMax(float n1, float n2) 
{
	int error = -1;

	if (n1 > n2 ) 
	{
		return n1;
	}
	else 
	{
		return n2;
	}

	return error;
}


float Svr::GetMin(float n1, float n2) 
{
	int error = -1;

	if (n1 < n2 ) 
	{
		return n1;
	}
	else 
	{
		return n2;
	}

	return error;
}
