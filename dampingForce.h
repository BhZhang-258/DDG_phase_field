#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticPlate &m_plate, double m_viscosity);
	~dampingForce();
	void computeFd(timeStepper &m_stepper);
	void computeJd(timeStepper &m_stepper);

	void setFirstJacobian(timeStepper &m_stepper);

private:
	elasticPlate *plate;
    

    double viscosity;

    Vector2d u, f;

    Vector2d jac;

    double ind;

};

#endif
