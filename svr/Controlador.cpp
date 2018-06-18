#include "Controlador.h"


//Constructor
Controlador::Controlador() {

}

//Destructor
Controlador::~Controlador() {

}


int Controlador::run() {

	CSVReader reader;
	vector<vector<float>> XData = reader.parse2DCsvFile("Data/XDataPdf.csv");
	vector<float> YData = reader.parse1DCsvFile("Data/YDataPdf.csv");

	//                     C    tol   K L.O
	Svr* maquina = new Svr(10, 0.001, 1, 1, XData, YData);

	maquina->MainRoutine();

	maquina->ComputeW(maquina->alphas, maquina->X, maquina->Y);

	cout << "w: ";
	maquina->PrintVector(maquina->w);
	cout << "b: " << maquina->b << endl;

	// cout << "New y1: " << maquina->Predict({4.8,3.4,1.9,0.2}) << endl;
	// cout << "New y2: " << maquina->Predict({6.8,2.8,4.8,1.4}) << endl;

	delete maquina;


    return 0;
}


