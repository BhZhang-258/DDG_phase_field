#ifndef INERTIALFORCE_H
#define INERTIALFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class inertialForce
{
public:
	inertialForce(elasticPlate &m_plate);
	~inertialForce();
	void computeFi( timeStepper &m_stepper);
	void computeJi( timeStepper &m_stepper);

	void setFirstJacobian( timeStepper &m_stepper);
	
private:
	elasticPlate *plate;
    			
    int ind1, ind2, mappedInd1, mappedInd2;	
    double f, jac;
};

#endif
