#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");			
	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");
	Possion = m_inputData.GetScalarOpt("Possion");
	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");

	alpha = m_inputData.GetScalarOpt("alpha");
	beta = m_inputData.GetScalarOpt("beta");
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDiscretePlate";
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		return;
	}

	if (timeStep == Nstep)
	{
		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);

			outfile << xCurrent(0) << " " << xCurrent(1) << endl;
		}
	}
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, Possion, deltaTime);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_elasticForce = new elasticForce(*plate, *stepper, alpha, beta);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	
	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_elasticForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();

	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;

}

void world::plateBoundaryCondition()
{

	for (int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		if (xCurrent(0) < -0.09)
		{
			plate->setVertexBoundaryCondition(xCurrent, i);
		}

		if (xCurrent(0) > 0.89)
		{
			//plate->setVertexBoundaryCondition(xCurrent, i);
		}
	}

	//plate->setVertexBoundaryCondition(plate->getVertex(0), 0);
	//plate->setVertexBoundaryCondition(plate->getVertex(3), 3);
}

void world::updateTimeStep()
{
	for (int i = 0; i < plate->nv; i++)
	{
		Vector2d xS = plate->getVertexStart(i);

		if (xS(0) < -0.09)
		{
			Vector2d xCurrent = plate->getVertex(i);
			//xCurrent(1) = xCurrent(1) - 0.1 * deltaTime;

			//plate->setVertexBoundaryCondition(xCurrent, i);
		}

	}

	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	// Start with a trial solution for our solution x
	plate->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_gravityForce->computeFg();
		m_elasticForce->computeFe();
		m_dampingForce->computeFd();

		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			m_gravityForce->computeJg();
			m_elasticForce->computeJe();
			m_dampingForce->computeJd();

			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter=" << iter << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		//timeStep = Nstep; // we are exiting
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

MatrixXd world::getScaledCoordinate(int i)
{
	MatrixXd xCurrent;

	xCurrent = MatrixXd::Zero(2,3);
	
	xCurrent.col(0) = plate->getVertex(plate->v_triangularElement[i].nv_1);
	xCurrent.col(1) = plate->getVertex(plate->v_triangularElement[i].nv_2);
	xCurrent.col(2) = plate->getVertex(plate->v_triangularElement[i].nv_3);
	
	return xCurrent * scaleRendering;
}

int world::numPair()
{
	return plate->ne;
}