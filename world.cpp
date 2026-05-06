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
	viscosityPhi = m_inputData.GetScalarOpt("viscosityPhi");
	Gc = m_inputData.GetScalarOpt("Gc");
	ell = m_inputData.GetScalarOpt("ell");
	eta = m_inputData.GetScalarOpt("eta");
	outputFreq = m_inputData.GetIntOpt("outputFreq");

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

	int systemRet = system("mkdir -p datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}
	// Per-step output files are handled in CoutData()
}

void world::CloseFile(ofstream &outfile)
{
	// Per-step files are closed in CoutData()
}

void world::CoutData(ofstream &outfile)
{
	if (saveData == false) 
	{
		return;
	}

	if (outputFreq <= 0) return;

	if (timeStep % outputFreq != 0)
	{
		return;
	}

	// Build filename with step number
	ostringstream name;
	name << "datafiles/step_" << timeStep << ".txt";

	ofstream stepFile(name.str().c_str());
	if (!stepFile.is_open())
	{
		cerr << "Failed to open output file: " << name.str() << endl;
		return;
	}
	stepFile.precision(10);

	// Header: time step, current time, number of nodes
	stepFile << "# timeStep " << timeStep << " currentTime " << currentTime << endl;
	stepFile << "# nv " << plate->nv << endl;
	stepFile << "# x y phi" << endl;

	// Write node data: current position (x, y) and phase field (phi)
	for (int i = 0; i < plate->nv; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);
		double phi_local = plate->getPhi(i);

		stepFile << xCurrent(0) << " " << xCurrent(1) << " " << phi_local << endl;
	}

	stepFile.close();
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, Possion, deltaTime);
	plateBoundaryCondition();

	plate->setup();
	plate->updateTimeStep();

	// set up stepper
	SolverConfig config_e;
	config_e.total_dof = plate->ndof_u;
	config_e.unconstrained_dof =  plate->uncons_u;
	config_e.fullToUncons = plate->fullToUnconsMap_u;
	config_e.isConstrained = plate->isConstrained_u;
	stepper_e = new timeStepper(config_e);

	SolverConfig config_phi;
	config_phi.total_dof = plate->ndof_phi;
	config_phi.unconstrained_dof =  plate->uncons_phi;
	config_phi.fullToUncons = plate->fullToUnconsMap_phi;
	config_phi.isConstrained = plate->isConstrained_phi;
	stepper_phi = new timeStepper(config_phi);

	// set up force
	// m_inertialForce = new inertialForce(*plate);
	// m_gravityForce = new externalGravityForce(*plate, gVector);
	m_elasticForce = new elasticForce(*plate, Gc, ell, eta);
	m_dampingForce = new dampingForce(*plate, viscosity, viscosityPhi);
	

	// set up first jacobian
	// m_inertialForce->setFirstJacobian(*stepper_e);
	m_elasticForce->setFirstJacobian(*stepper_e);
	m_dampingForce->setFirstJacobian(*stepper_e);
	stepper_e->finalizePattern();
	// set up pardiso
	stepper_e->first_time_PARDISO_setup();

	if (plate->uncons_phi > 0)
	{
		m_elasticForce->setFirstJacobian_phi(*stepper_phi);
		m_dampingForce->setFirstJacobian_phi(*stepper_phi);
		stepper_phi->finalizePattern();
		stepper_phi->first_time_PARDISO_setup();
	}

	// a trial solution for our solution x
	// for (int i=0; i < plate->nv ; i++)
	// {
	// 	Vector2d xCurrent = plate->getVertex(i);
	// 	Vector2d xtemp;
	// 	xtemp(0) = xCurrent(0) /2; 
	// 	xtemp(1) = xCurrent(1) *2 + xCurrent(0); 
	// 	plate->x.segment(2*i, 2) = xtemp;
	// }
	// // test computeFe
	// m_elasticForce->computeFe(*stepper_e);
	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;

	// Write mesh connectivity once for post-processing
	if (saveData)
	{
		ofstream meshFile("datafiles/mesh.txt");
		meshFile.precision(10);
		meshFile << "# nv ne" << endl;
		meshFile << plate->nv << " " << plate->ne << endl;
		meshFile << "# element connectivity (0-based node indices)" << endl;
		for (int i = 0; i < plate->ne; i++)
		{
			meshFile << plate->v_triangularElement[i].nv_1 << " "
					 << plate->v_triangularElement[i].nv_2 << " "
					 << plate->v_triangularElement[i].nv_3 << endl;
		}
		meshFile.close();
	}


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
	if (maxVal <= 0.7)
	{
		stretchingDistance += 1e-1 * deltaTime;
	}
	else
	{
		stretchingDistance += 0.5e-2 * deltaTime;
	}
	
	for (int i = 0; i < loadid.size(); i++)
	{
		Vector2d xCurrent = plate->getVertex(loadid[i]);
		Vector2d xStart = plate->getVertexStart(loadid[i]);
		xCurrent(1) = xStart(1) + stretchingDistance;

		plate->setVertexBoundaryCondition(xCurrent, loadid[i]);
	}
	
	
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;
	int iter_phi = 0;

	// Start with a trial solution for our solution x
	// We use phi instead u to advance the displacement update.
	// plate->updateGuess(); // x = x0 + u * dt

	const double armijo_c = 1e-4;
	const int    max_ls   = 4;
		
	while (solved == false)
	{
		plate->prepareForIteration();

		stepper_e->setZero();

		// m_inertialForce->computeFi(*stepper_e);
		m_elasticForce->computeFe(*stepper_e);
		// m_dampingForce->computeFd(*stepper_e); 
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

		if (solved == false)
		{
			double R0 = normf;  // residual before step

			// m_inertialForce->computeJi(*stepper_e);
			m_elasticForce->computeJe(*stepper_e);
			// m_dampingForce->computeJd(*stepper_e);

			stepper_e->integrator(); // Solve equations of motion
			VectorXd direction = stepper_e->GlobalMotionVec;

			// --- Backtracking line search ---
			plate->backupState_u();
			double alpha = 1.0;
			bool ls_ok = false;
			for (int ls = 0; ls < max_ls; ls++)
			{
				plate->restoreState_u();
				plate->updateNewtonMethod(direction, alpha);

				stepper_e->setZero();
				// m_inertialForce->computeFi(*stepper_e);
				m_elasticForce->computeFe(*stepper_e);
				// m_dampingForce->computeFd(*stepper_e);
				double Rnew = stepper_e->GlobalForceVec.norm();

				if (Rnew <= (1.0 - armijo_c * alpha) * R0)
				{
					ls_ok = true;
					break;
				}
				alpha *= 0.5;
			}
			// If no sufficient decrease, still keep last trial (best effort)
			// State is at the accepted α (or last trial)
			
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	if (plate->uncons_phi > 0)
	{
		m_elasticForce->updateHistoryField();

		double normphi = forceTol * 10.0;
		double normphi0 = 0.0;
		bool solved_phi = false;
		iter_phi = 0;

		while (solved_phi == false)
		{
			plate->prepareForIteration();

			stepper_phi->setZero();
			m_elasticForce->computeFphi(*stepper_phi);
			m_dampingForce->computeFd_phi(*stepper_phi);

			normphi = stepper_phi->GlobalForceVec.norm();

			if (iter_phi == 0)
			{
				normphi0 = normphi;
			}

			if (normphi <= forceTol)
			{
				solved_phi = true;
			}
			else if (iter_phi > 0 && normphi <= normphi0 * stol)
			{
				solved_phi = true;
			}

			if (solved_phi == false)
			{
				double R0_phi = normphi;

				m_dampingForce->computeJd_phi(*stepper_phi);
				stepper_phi->integrator();
				VectorXd direction_phi = stepper_phi->GlobalMotionVec;

				// --- Backtracking line search for phase-field ---
				plate->backupState_phi();
				double alpha_phi = 1.0;
				for (int ls = 0; ls < max_ls; ls++)
				{
					plate->restoreState_phi();
					plate->updateNewtonMethod_phi(direction_phi, alpha_phi);

					stepper_phi->setZero();
					m_elasticForce->computeFphi(*stepper_phi);
					m_dampingForce->computeFd_phi(*stepper_phi);
					double Rnew_phi = stepper_phi->GlobalForceVec.norm();
					if (Rnew_phi <= (1.0 - armijo_c * alpha_phi) * R0_phi)
					{
						break;
					}
					alpha_phi *= 0.5;
				}

				iter_phi++;
			}
			
			if (iter_phi > maxIter)
			{
				cout << "Error. Phase-field could not converge. Exiting.\n";
				break;
			}
		}
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << " iter_u=" << iter << " iter_phi=" << iter_phi << endl;
		maxVal = plate->phi.maxCoeff();
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



double world::getPhi(int i)
{
	double phi_1 = plate->getPhi(plate->v_triangularElement[i].nv_1);
	double phi_2 = plate->getPhi(plate->v_triangularElement[i].nv_2);
	double phi_3 = plate->getPhi(plate->v_triangularElement[i].nv_3);

	return (phi_1 + phi_2 + phi_3) / 3;
}