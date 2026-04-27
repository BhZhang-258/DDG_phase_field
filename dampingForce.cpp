#include "dampingForce.h"

dampingForce::dampingForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_viscosity)
{
	plate = &m_plate;
    stepper = &m_stepper;

	viscosity = m_viscosity;
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	for (int i = 0; i < plate->ndof; i++)
	{
		double localForce = - viscosity * plate->massArray(i) * (plate->x(i) - plate->x0(i)) / plate->dt;

		stepper->addForce(i, - localForce);
	}
}

void dampingForce::computeJd()
{
	for (int i = 0; i < plate->ndof; i++)
	{
		double localJacob = - viscosity * plate->massArray(i) / plate->dt;

		stepper->addJacobian(i, i, - localJacob);
	}
	
}

void dampingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->ndof; i++)
	{
		stepper->addJacobian(ind, ind, 1);
	}
}
