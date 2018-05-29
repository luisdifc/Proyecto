#ifndef CONTROLADOR_H
#define CONTROLADOR_H

#include <iostream>
#include "Svr.h"

using namespace std;

class Controlador {

//Aributos
public:
	

//Metodos
public:
    Controlador();
    virtual ~Controlador();
    int run(); //corre el juego
};

#endif // CONTROLADOR_H
