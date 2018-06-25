#include "Controlador.h"


//Constructor
Controlador::Controlador() {

}

//Destructor
Controlador::~Controlador() {

}


int Controlador::run() {

	CSVReader reader;

	// vector<vector<float>> XData = reader.parse2DCsvFile("Data/flowers/XData.csv");
	// vector<float> YData = reader.parse1DCsvFile("Data/flowers/YData.csv");

	// vector<vector<float>> XData = reader.parse2DCsvFile("Data/pdfDataset/XDataPdf.csv");
	// vector<float> YData = reader.parse1DCsvFile("Data/pdfDataset/YDataPdf.csv");

	vector<vector<float>> XData = reader.parse2DCsvFile("Data/salary/Data_X.csv");
	vector<float> YData = reader.parse1DCsvFile("Data/salary/Data_Y.csv");

	//                     C  tol  K L.O
	Svr* maquina = new Svr(5, 2, 1, 1, XData, YData); //100 

	maquina->MainRoutine();
	maquina->ComputeW(maquina->alphas, maquina->X, maquina->Y);
	cout << "w: ";
	maquina->PrintVector(maquina->w);
	cout << "b: " << maquina->b << endl;
	cout << "Alpas: ";
	maquina->PrintVector(maquina->alphas);
	cout << endl;


	// cout << "New y1(+): " << maquina->Predict({6.2,2.9,4.3,1.3}) << endl;
	// cout << "New y2(+): " << maquina->Predict({6.1,3,4.6,1.4}) << endl;
	// cout << "New y3(+): " << maquina->Predict({6.8,2.8,4.8,1.4}) << endl;
	// cout << "New y4(-): " << maquina->Predict({4.9,3,1.4,0.2}) << endl;
	// cout << "New y5(-): " << maquina->Predict({4.7,3.2,1.3,0.2}) << endl;
	// cout << "New y6(-): " << maquina->Predict({4.3,3,1.1,0.1}) << endl;

	// cout << "New y1(+): " << maquina->Predict({1,11}) << endl;
	// cout << "New y2(+): " << maquina->Predict({2,9}) << endl;
	// cout << "New y3(+): " << maquina->Predict({8,8}) << endl;
	// cout << "New y4(+): " << maquina->Predict({6,10}) << endl;
	// cout << "New y5(-): " << maquina->Predict({11,2}) << endl;
	// cout << "New y5(-): " << maquina->Predict({5,5}) << endl;
	// cout << "New y5(-): " << maquina->Predict({9,3}) << endl;
	// cout << "New y5(-): " << maquina->Predict({4,3}) << endl;

	// cout << "Predicted y: " << maquina->PredictRegression({29,19}) << endl;


	delete maquina;


    return 0;
}


