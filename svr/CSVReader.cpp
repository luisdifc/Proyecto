#include "CSVReader.h"

//Constructor
CSVReader::CSVReader() {

}

//Destructor
CSVReader::~CSVReader() {

}
 
/**
 * Reads csv file into table, exported as a vector of vector of floats.
 * @param inputFileName input file name (full path).
 * @return data as vector of vector of floats.
 */
vector<vector<float>> CSVReader::parse2DCsvFile(string inputFileName) {
 
    vector<vector<float> > data;
    ifstream inputFile(inputFileName);
    int l = 0;
 
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<float> record;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
 
            data.push_back(record);
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}
 

vector<float> CSVReader::parse1DCsvFile(string inputFileName) {
 
    vector<float> data;
    ifstream inputFile(inputFileName);
    int l = 0;
 
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<float> record;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    data.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}