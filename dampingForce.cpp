#include "dampingForce.h"

dampingForce::dampingForce(elasticPlate &m_plate,  double m_viscosity)
{
	plate = &m_plate;

	viscosity = m_viscosity;
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd(timeStepper &m_stepper)
{
	for (int i = 0; i < plate->ndof_u; i++)
	{
		double localForce = - viscosity * plate->massArray(i) * (plate->x(i) - plate->x0(i)) / plate->dt;

		m_stepper.addForce(i, - localForce);
	}
}

void dampingForce::computeJd(timeStepper &m_stepper)
{
	for (int i = 0; i < plate->ndof_u; i++)
	{
		double localJacob = - viscosity * plate->massArray(i) / plate->dt;

		m_stepper.addJacobian(i, i, - localJacob);
	}
	
}

void dampingForce::setFirstJacobian(timeStepper &m_stepper)
{
	for (int i = 0; i < plate->ndof_u; i++)
	{
		m_stepper.addFirstJacobian(i, i);
	}
}
