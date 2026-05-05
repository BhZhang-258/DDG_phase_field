#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include input file and option
#include "setInput.h"

// include elastic Plate class
#include "elasticPlate.h"

// include time stepper
#include "timeStepper.h"

// include force
#include "inertialForce.h"
#include "externalGravityForce.h"
#include "elasticForce.h"
#include "dampingForce.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	void setPlateStepper();

	void updateTimeStep();

	int simulationRunning();

	int numPair();
	MatrixXd getScaledCoordinate(int i);

	double getPhi(int i);
		
private:

	// physical parameters
	bool render;
	bool saveData;
	double deltaTime;
	double totalTime;
	double YoungM;
	double density;
	double thickness;
	double Possion;
	double stol;
	double forceTol;
	double scaleRendering;
	int maxIter;
	Vector3d gVector;
	double viscosity;
	
	double Gc;
	double ell;
	double eta;

	int outputFreq;

	int Nstep;
	int timeStep;

	double characteristicForce;

	double currentTime;

	void plateBoundaryCondition();

	// Plate
	elasticPlate *plate;

	// stepper
	timeStepper *stepper_e;
	timeStepper *stepper_phi;

	// force
	inertialForce *m_inertialForce;
	externalGravityForce *m_gravityForce;
	elasticForce *m_elasticForce;
	dampingForce *m_dampingForce;

	void computeRactionForce();
	std::vector<int> loadid;
	VectorXd reForce;
	double maxVal = 0.0;

	double stretchingDistance;
};

#endif
