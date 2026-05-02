#ifndef EXTERNALGRAVITYFORCE_H
#define EXTERNALGRAVITYFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalGravityForce
{
public:
	externalGravityForce(elasticPlate &m_plate,  Vector3d m_gVector);
	~externalGravityForce();
	void computeFg(timeStepper &m_stepper);
	void computeJg(timeStepper &m_stepper);

	void setGravity();

	Vector3d gVector;
	
private:
	elasticPlate *plate;
	
    VectorXd massGravity;
};

#endif
