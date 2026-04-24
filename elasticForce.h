#ifndef ELASTICFORCE_H
#define ELASTICFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticForce
{
public:
	elasticForce(elasticPlate &m_plate, timeStepper &m_stepper, 
    double m_alpha, double m_beta);
	~elasticForce();
	void computeFe();
	void computeJe();
    void setFirstJacobian();

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double alpha;
    double beta;

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

    void triElement(
    const Eigen::Matrix<double, 2, 3>& x,
    const Eigen::Matrix<double, 2, 3>& X,
    double materialParameter1,
    double materialParameter2,
    double& E,
    VectorXd& f,
    MatrixXd& H);
};

#endif
