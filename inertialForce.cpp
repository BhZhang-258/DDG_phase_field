#include "inertialForce.h"

inertialForce::inertialForce(elasticPlate &m_plate)
{
	plate = &m_plate;
}

inertialForce::~inertialForce()
{
	;
}

void inertialForce::computeFi( timeStepper &m_stepper)
{
	for (int i = 0; i < plate->ndof_u; i++)
	{
		f = plate->massArray[i] * (plate->x[i] - plate->x0[i]) / ((plate->dt) *(plate->dt))
				- (plate->massArray[i] * plate->u[i])/(plate->dt);

		m_stepper.addForce(i, f);
	}
}

void inertialForce::computeJi(timeStepper &m_stepper)
{
	for (int i=0; i<plate->ndof_u; i++)
    {
		jac = plate->massArray(i)/ ((plate->dt) *(plate->dt));
		m_stepper.addJacobian(i, i, jac);
	}
}

void inertialForce::setFirstJacobian(timeStepper &m_stepper)
{
	for (int i=0; i<plate->ndof_u; i++)
	{
		m_stepper.addFirstJacobian(i, i);
	}
}
