#include "Controlador.h"


//Constructor
Controlador::Controlador() {

}

//Destructor
Controlador::~Controlador() {

}


int Controlador::run() {
	Svr* maquina = new Svr(10, 0.001, 1, 0);
	maquina->MainRoutine();

	maquina->PrintVector(maquina->alphas);
	maquina->PrintVector(maquina->Y);
	maquina->PrintMatrix(maquina->X);

	maquina->PrintVector(maquina->w);

	vector<float> weights = maquina->ComputeW(maquina->alphas, maquina->X, maquina->Y);
	maquina->PrintVector(weights);
	cout << "b: " << maquina->b << endl;

	delete maquina;
	

    return 0;
}


