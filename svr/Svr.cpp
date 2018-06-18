#include "Svr.h"


//Constructor
Svr::Svr(float C, float tol, int kernelType, int useLinearOptim, int size) 
{
    cout << "\nConstruyendo un SVR..!\n" << endl;

    GiveSizeMatrix(14, 2, this->X);
    this->X = {{8,7},{4,10},{9,7},{7,10},{9,6},{4,8},{10,10},{2,7},{8,3},{7,5},{4,4},{4,6},{1,3},{2,5}};
    this->PrintMatrix(this->X);

    this->Y.resize(14);
    this->Y = {1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1};
    this->PrintVector(this->Y);

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

	if (this->useLinearOptim == 1) 
	{
		//Equation 1
		return (this->DotProduct(this->w, this->X[i]) - this->b);
	} 
	else 
	{
		for (int index = 0; index < this->m.size(); index++) 
		{
			//Equation 10
			sum += (this->alphas[index] * this->Y[index] * this->Kernel(this->X[index], this->X[i]));
		}
		return sum - this->b;
	}

	return error;
}


//Try to solve the problem analitically 
int Svr::TakeStep(int i1, int i2, float a2, float y2, float E2, vector<float> x2) 
{
	//return = 0 false, return 1 = true

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

	s = y1 * y2;

	//Compute the bounds of the new alpha2
	if (y1 != y2) 
	{
		//Equation 13
		L = GetMax(0.0, a2 - a1);
		H = GetMin(this->C, this->C + a2 - a1);
	}
	else 
	{
		//Equation 14
		L = GetMax(0.0, a2 + a1 - this->C);
		H = GetMin(this->C, a2 + a1);
	}

	if (L == H) 
	{
		return 0;
	}

	k11 = this->Kernel(x1, x1);
	k12 = this->Kernel(x1, x2);
	k22 = this->Kernel(x2, x2);

	//Compute the second derivative of the objective function along the diagonal
	//Equation 15
	eta = k11 + k22 - 2 * k12;

	if (eta > 0 ) 
	{
		//Equation 16
		a2New = a2 + y2 * (E1 - E2) / eta;

		//Clip the new alpha so that is stays at the end of the line
		//Equation 17
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
		//Under unusual circumstances eta will not be positive
		//Equation 19
		f1 = y1 * (E1 + this->b) - a1 * k11 - s * a2 * k12;
		f2 = y2 * (E2 + this->b) - s * a1 * k12 - a2 * k22;
		L1 = a1 + s * (a2 - L);
		H1 = a1 + s * (a2 - H);
		Lobj = L1 * f1 + L * f2 + 0.5 * pow(L1, 2) * k11 + 0.5 * pow(L, 2) * k22 + s * L * L1 * k12;
		Hobj = H1 * f1 + H * f2 + 0.5 * pow(H1, 2) * k11 + 0.5 * pow(H, 2) * k22 + s * H * H1 * k12;
		
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
				a2New = a2;
			}
		}
	}

	//If alpha2 did not change enough the algorithm returns without updating the multipliers
	if (abs(a2New - a2) < (this->eps * (a2New + a2 + this->eps))) 
	{
		return 0;
	}

	//Equation 18
	a1New = a1 + s * (a2 - a2New);
	newB = this->ComputeB(E1, a1, a1New, a2New, k11, k12, k22, y1, y2, a2, E2);
	deltaB = newB - this->b;
	this->b = newB;

	//Equation 22
	if (this->useLinearOptim) 
	{
		this->w = VectorSum(this->w, VectorByScalar(x1, (y1*(a1New - a1))));
		this->w = VectorSum(this->w, VectorByScalar(x2, (y2*(a2New - a2))));
	}

	//Update the error cache using the new Lagrange multipliers 
	delta1 = y1 * (a1New - a1);
	delta2 = y2 * (a2New - a2);

	//Update the error cache
	for (int index = 0; index < this->m.size(); index++) 
	{
		if (0 < this->alphas[index] && this->alphas[index] < this->C) 
		{
			this->errors[index] += delta1 * this->Kernel(x1, this->X[index]) + delta2 * this->Kernel(x2, this->X[index]) - deltaB;
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

	//Equation 20
	b1 = E1 + y1 * (a1New - a1) * k11 + y2 * (a2New - a2) * k12 + this->b;

	//Equation 21
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


vector<float> Svr::ComputeW(vector<float> multipliers, vector< vector<float> > X, vector<float> y) 
{
	vector<float> w;
	float sum = 0.0;

	for (int index = 0; index < y.size(); index++) 
	{
		w = VectorSum(w, VectorByScalar(X[index], y[index]));
	}

	return w;
}


float Svr::FirstHeuristic() 
{
	int numChanged = 0;
	vector<int> nonBoundIndexes = GetNonBoundIndexes();
	
	for (int index = 0; index < nonBoundIndexes.size(); index++) 
	{
		numChanged += this->ExamineExample(nonBoundIndexes[index]);
	}

	return numChanged;
}


float Svr::SecondHeuristic(vector<int> nonBoundIndices, float E2) 
{
	float i1 = -1;
	float E1 = 0.0;
	float step = 0.0;
	float max = 0.0;

	if (nonBoundIndices.size() > 1) 
	{
		max = 0;
		for (int index = 0; index < nonBoundIndices.size(); index++) 
		{
			E1 = this->errors[nonBoundIndices[index]] - this->Y[nonBoundIndices[index]];
			step = abs(E1 - E2);
			if (step > max) 
			{
				max = step;
				i1 = nonBoundIndices[index];
			}
		}
	}

	return i1;
}


vector<int> Svr::GetNonBoundIndexes() 
{
	vector<int> result;
	for (int index = 0; index < this->alphas.size(); index++) 
	{
		if (this->alphas[index] > 0 && this->alphas[index] < this->C) 
		{
			result.push_back(index);
		}
		else 
		{
			//continue
		}
	}

	return result;
}

int Svr::ExamineExample (int i2) 
{
	float y2 = this->Y[i2];
	float a2 = this->alphas[i2];
	vector<float> x2 = this->X[i2];
	float E2 = this->GetError(i2);

	float r2 = 0.0;
	vector<int> nonBoundIndexes;
	int i1 = 0;
	int rand = 0;

	r2 = E2 * y2;

	if (! ( (r2 < (-1*this->tol) && a2 < this->C) || (r2 > this->tol && a2 > 0) ) ) 
	{
		//The KKT conditions are met, SMO looks at another example
		return 0;
	}

	//Second Heuristic A: choose the Lagrange multiplier which maximezes the absolute error
	nonBoundIndexes = GetNonBoundIndexes();
	i1 = SecondHeuristic(nonBoundIndexes, E2);

	if (i1 >= 0 && this->TakeStep(i1, i2, a2, y2, E2, x2) == 1) 
	{
		return 1;
	}

	
	//Second Heuristic B: look for examples making positive progress by looping over all non zero and non C alpha
	//starting at a random point 
	if (nonBoundIndexes.size() > 0) 
	{
		rand = RandNumGenerator(0, nonBoundIndexes.size());
		for (int index = rand; index < nonBoundIndexes.size(); index++) 
		{
			if (this->TakeStep(nonBoundIndexes[index], i2, a2, y2, E2, x2)) 
			{
				return 1;
			}
		}

		for (int index = 0; index < rand; index++) 
		{
			if (this->TakeStep(nonBoundIndexes[index], i2, a2, y2, E2, x2)) 
			{
				return 1;
			}
		}
	}

	rand = RandNumGenerator(0, this->m.size());
	for (int index = rand; index < m.size(); index++) 
	{
		if (this->TakeStep(index, i2, a2, y2, E2, x2)) 
		{
			return 1;
		}
	}

	for (int index = 0; index < rand; index++) 
	{
		if (this->TakeStep(index, i2, a2, y2, E2, x2)) 
		{
			return 1;
		}
	}
	//In extremely degenerate circumstances the SMO skips the first example
	return 0;
}

void Svr::MainRoutine() 
{
	float numChanged = 0.0;
	int examineAll = 1;

	while (numChanged > 0 || examineAll) 
	{
		numChanged = 0.0;

		if (examineAll == 1) 
		{
			for (int index = 0; this->m.size(); index++) 
			{
				numChanged += this->ExamineExample(index);
			}
		}
		else 
		{
			numChanged += this->FirstHeuristic();
		}

		if (examineAll == 1) 
		{
			examineAll = 0;
		}
		else 
		{
			if (numChanged == 0) 
			{
				examineAll = 1;
			}
		}
	}
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


vector<float> Svr::VectorByScalar (vector<float> v1, float scalar) 
{
	vector<float> result;

	for (int index = 0; index < v1.size(); index++) 
	{
		result[index] = v1[index] * scalar;
	}	

	return result;	
}


vector<float> Svr::VectorSum(vector<float> v1, vector<float> v2) 
{
	vector<float> result;

	for (int index = 0; index < v1.size(); index++) 
	{
		result[index] = v1[index] + v2[index];
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


int Svr::RandNumGenerator(int n1, int n2) 
{
	int randNum = 0;

	srand((unsigned)time(0)); 
    randNum = (rand()%n2)+n1; 

    return randNum;	
}


void Svr::PrintVector(vector<float> vector) 
{
	cout << "{";
	for (int index = 0; index < vector.size(); index++) 
	{
		cout << vector[index] << ",";
	}
	cout << "}\n";
}


void Svr::PrintMatrix(vector< vector<float> > array) 
{
	for (int index1 = 0; index1 < array.size(); index1++) 
	{
		for (int index2 = 0; index2 < array[index1].size(); index2+=2) 
		{
				cout << "{" << array[index1][index2] << ",";
				cout << array[index1][index2+1] << "}" << ",";	
		}
	}
	cout << "\n";
}

void Svr::GiveSizeMatrix(int i, int j, vector< vector<float> > &array) 
{
	array.resize(i);

	for (int index1 = 0; index1 < i; index1++) 
	{
		array[index1].resize(j);
	}
}