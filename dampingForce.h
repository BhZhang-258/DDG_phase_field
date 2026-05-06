#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticPlate &m_plate, double m_viscosity, double m_viscosityPhi);
	~dampingForce();
	void computeFd(timeStepper &m_stepper);
	void computeJd(timeStepper &m_stepper);

	void setFirstJacobian(timeStepper &m_stepper);

	// --- Phase-field damping ---
	void computeFd_phi(timeStepper &m_stepper);
	void computeJd_phi(timeStepper &m_stepper);
	void setFirstJacobian_phi(timeStepper &m_stepper);

private:
	elasticPlate *plate;
    
    double viscosity;
    double viscosityPhi;

    VectorXd phiNodalArea;   // lumped nodal area for phase-field damping

    void computePhiNodalArea();

    Vector2d u, f;
    Vector2d jac;
    double ind;
};

#endif
