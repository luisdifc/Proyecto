#include "Controlador.h"


//Constructor
Controlador::Controlador() {

}

//Destructor
Controlador::~Controlador() {

}


int Controlador::run() {

	CSVReader reader;
	vector<vector<float>> XData = reader.parse2DCsvFile("XData.csv");
	vector<float> YData = reader.parse1DCsvFile("YData.csv");

	Svr* maquina = new Svr(10, 0.001, 1, 0, XData, YData);

	maquina->MainRoutine();

	vector<float> weights = maquina->ComputeW(maquina->alphas, maquina->X, maquina->Y);
	maquina->PrintVector(weights);
	cout << "b: " << maquina->b << endl;

	delete maquina;


    return 0;
}


