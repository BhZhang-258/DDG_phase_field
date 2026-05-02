#include "world.h"
#include <chrono>

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
	Gc = m_inputData.GetScalarOpt("Gc");
	ell = m_inputData.GetScalarOpt("ell");
	eta = m_inputData.GetScalarOpt("eta");

	stretchingDistance = 0.0;
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

	if ( timeStep % 20 != 0)
	{
		return;
	}

	
	// computeRactionForce();

	// double totalForce;

	// totalForce = 0.0;

	// for (int i = 0; i < plate->nv; i++)
	// {
	// 	Vector2d xS = plate->getVertexStart(i);

	// 	if (abs(xS(1)-1.0) < 1e-8)
	// 	{
	// 		totalForce = totalForce + reForce(2 * i + 1);
	// 	}
	// }

	// outfile << stretchingDistance << " " << totalForce << endl;


	/*
	//if (timeStep == Nstep)
	{
		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);
			double phi_local = plate->getPhi(i);

			outfile << stretchingDistance << " " << xCurrent(0) << " " << xCurrent(1) << " " << phi_local << endl;
		}
	}
	*/

}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, Possion, deltaTime);
	plateBoundaryCondition();

	plate->setup();

	// set up stepper
	SolverConfig config_e;
	config_e.total_dof = plate->ndof_u;
	config_e.unconstrained_dof =  plate->uncons_u;
	config_e.fullToUncons = plate->fullToUnconsMap_u;
	config_e.isConstrained = plate->isConstrained_u;


	stepper_e = new timeStepper(config_e);

	SolverConfig config_phi;

	// stepper_phi = new timeStepper(config_phi);

	// set up force
	m_inertialForce = new inertialForce(*plate);
	m_gravityForce = new externalGravityForce(*plate, gVector);
	m_elasticForce = new elasticForce(*plate, Gc, ell, eta);
	m_dampingForce = new dampingForce(*plate, viscosity);
	
	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian(*stepper_e);
	m_elasticForce->setFirstJacobian(*stepper_e);
	m_dampingForce->setFirstJacobian(*stepper_e);
	stepper_e->finalizePattern();
	// set up pardiso
	stepper_e->first_time_PARDISO_setup();


	for (int i=0; i < plate->nv ; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);
		Vector2d xtemp;
		xtemp(0) = xCurrent(0) /2; 
		xtemp(1) = xCurrent(1) *2 + xCurrent(0); 
		plate->x.segment(2*i, 2) = xtemp;
	}

	m_elasticForce->computeFe(*stepper_e);
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

		if ( abs(xCurrent(1) - 0.0) < 1e-8 )
		{
			plate->setVertexBoundaryCondition(xCurrent, i);

			plate->setPhiBoundaryCondition(0.0, i);
			
		}

		if ( abs(xCurrent(1) - 1.0) < 1e-8 )
		{
			plate->setVertexBoundaryCondition(xCurrent, i);

			plate->setPhiBoundaryCondition(0.0, i);
			loadid.push_back(i);
		}
		// plate->setPhiBoundaryCondition(0.0, i);

		
	}

}

void world::updateTimeStep()
{
	
	stretchingDistance = stretchingDistance + 1e-2 * deltaTime;
		
	

	for (int i = 0; i < loadid.size(); i++)
	{
		
			Vector2d xCurrent = plate->getVertex(loadid[i]);
			Vector2d xStart = plate->getVertexStart(loadid[i]);
			xCurrent(1) = xStart(1) + stretchingDistance;

			plate->setVertexBoundaryCondition(xCurrent, loadid[i]);
		
	}
	// cerr<< "stretchingDistance: " << stretchingDistance << endl;
	// cerr<< "current position: " << loadid.size() << endl;
	// Vector2d xCurrent =  plate->getVertex(loadid[0]);
	// cerr<< "current position: " << xCurrent.transpose() << endl;
	
	

	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	// Start with a trial solution for our solution x
	plate->updateGuess(); // x = x0 + u * dt
		
	while (solved == false)
	{
		plate->prepareForIteration();

		stepper_e->setZero();

		m_inertialForce->computeFi(*stepper_e);
		m_elasticForce->computeFe(*stepper_e);
		m_dampingForce->computeFd(*stepper_e);

		normf = stepper_e->GlobalForceVec.norm();

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
			m_inertialForce->computeJi(*stepper_e);
			//m_gravityForce->computeJg();
			m_elasticForce->computeJe(*stepper_e);
			m_dampingForce->computeJd(*stepper_e);

			stepper_e->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper_e->GlobalMotionVec); // new q = old q + Delta q
			int start = plate->nv * 2;
			int len   = plate->nv;   // 因为 3nv - 1 - 2nv + 1 = nv

			maxVal = plate->x.segment(start, len).maxCoeff();

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
		cout << "maxVal: " << maxVal << endl;
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

	xCurrent(0,0) = xCurrent(0,0) -0.5;
	xCurrent(0,1) = xCurrent(0,1) -0.5;
	xCurrent(0,2) = xCurrent(0,2) -0.5;
	xCurrent(1,0) = xCurrent(1,0) -0.5;
	xCurrent(1,1) = xCurrent(1,1) -0.5;
	xCurrent(1,2) = xCurrent(1,2) -0.5;
	
	return xCurrent * scaleRendering;
}

int world::numPair()
{
	return plate->ne;
}

// void world::computeRactionForce()
// {
// 	m_elasticForce->computeFe(*stepper_e);

// 	reForce = VectorXd::Zero(plate->ndof_u);

// 	reForce = m_elasticForce->reForce;
// }

double world::getPhi(int i)
{
	double phi_1 = plate->getPhi(plate->v_triangularElement[i].nv_1);
	double phi_2 = plate->getPhi(plate->v_triangularElement[i].nv_2);
	double phi_3 = plate->getPhi(plate->v_triangularElement[i].nv_3);

	return (phi_1 + phi_2 + phi_3) / 3;
}