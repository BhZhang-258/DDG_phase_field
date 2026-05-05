#include "externalGravityForce.h"

externalGravityForce::externalGravityForce(elasticPlate &m_plate,Vector3d m_gVector)
{
	plate = &m_plate;
	gVector = m_gVector;
	setGravity();
}

externalGravityForce::~externalGravityForce()
{
	;
}

void externalGravityForce::computeFg(timeStepper &m_stepper)
{
	for (int i=0; i < plate->ndof_u; i++)
	{
		m_stepper.addForce(i, -massGravity[i]); // subtracting gravity force
	}	
}

void externalGravityForce::computeJg(timeStepper &m_stepper)
{
	;
}

void externalGravityForce::setGravity()
{
	massGravity = VectorXd::Zero(plate->ndof_u);

	for (int i = 0; i < plate->nv; i++)
	{
		for (int k = 0; k < 2; k++)
		{
			int ind = 2 * i + k;
			massGravity[ind] = gVector[k] * plate->massArray[ind];
		}
	}
}
