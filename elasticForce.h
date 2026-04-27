#ifndef ELASTICFORCE_H
#define ELASTICFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticForce
{
public:
	elasticForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_Gc, double m_ell, double m_eta);
	~elasticForce();
	void computeFe();
	void computeJe();
    void setFirstJacobian();

    VectorXd reForce;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double alpha_input;
    double beta_input;

    double Gc_input;
    double ell_input;
    double eta_input;

    double cotangent2D(const Vector2d& a, const Vector2d& b);
    void closestRotationAndY2D(
        const Matrix2d& F,
        Matrix2d& R,
        Matrix2d& Y);
    void setBlock(MatrixXd& H, int bi, int bj, const Matrix2d& A);
    Matrix2d secondPartMat(
        const Vector2d& si,
        const Vector2d& sj,
        double A0,
        double trY);
    Matrix2d abSplitMat(
        const Vector2d& ri,
        const Vector2d& rj,
        double A0);

    void triElementElasticOnly(
        const Eigen::Matrix<double, 2, 3>& x,
        const Eigen::Matrix<double, 2, 3>& X,
        double alpha,
        double beta,
        double& Eel,
        VectorXd& fel,
        MatrixXd& Hel);
    void phaseFieldSurfaceEnergyDerivatives(
        const MatrixXd& referenceShape,   // 2x3
        const Vector3d& phi,              // [phi1, phi2, phi3]
        double Gc,
        double ell,
        double& Efrac,
        Vector3d& gradPhi,                // dE/dphi
        Matrix3d& hessPhi);                // d2E/dphi2
};

#endif
