#include "dampingForce.h"

dampingForce::dampingForce(elasticPlate &m_plate,  double m_viscosity, double m_viscosityPhi)
{
	plate = &m_plate;

	viscosity    = m_viscosity;
	viscosityPhi = m_viscosityPhi;

	computePhiNodalArea();
}

void dampingForce::computePhiNodalArea()
{
	phiNodalArea = VectorXd::Zero(plate->nv);

	for (int i = 0; i < plate->ne; i++)
	{
		double A3 = plate->v_triangularElement[i].area / 3.0;
		phiNodalArea(plate->v_triangularElement[i].nv_1) += A3;
		phiNodalArea(plate->v_triangularElement[i].nv_2) += A3;
		phiNodalArea(plate->v_triangularElement[i].nv_3) += A3;
	}
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

// ===================== Phase-field viscous damping =====================
// F_damp_phi = η_phi * A_i * (φ_i - φ_i_old) / dt
// J_damp_phi = η_phi * A_i / dt   (diagonal)

void dampingForce::computeFd_phi(timeStepper &m_stepper)
{
	const double inv_dt = 1.0 / plate->dt;

	for (int i = 0; i < plate->nv; i++)
	{
		double dphi = plate->phi(i) - plate->phi_old(i);
		double Fd   = viscosityPhi * phiNodalArea(i) * dphi * inv_dt;

		m_stepper.addForce(i, Fd);
	}
}

void dampingForce::computeJd_phi(timeStepper &m_stepper)
{
	const double inv_dt = 1.0 / plate->dt;

	for (int i = 0; i < plate->nv; i++)
	{
		double Jd = viscosityPhi * phiNodalArea(i) * inv_dt;

		m_stepper.addJacobian(i, i, Jd);
	}
}

void dampingForce::setFirstJacobian_phi(timeStepper &m_stepper)
{
	for (int i = 0; i < plate->nv; i++)
	{
		m_stepper.addFirstJacobian(i, i);
	}
}
