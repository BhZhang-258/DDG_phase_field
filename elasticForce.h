#ifndef ELASTICFORCE_H
#define ELASTICFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "GeometryUtils.h"
#include "timeStepper.h"

class elasticForce
{
public:
	elasticForce(elasticPlate &m_plate,  double m_Gc, double m_ell, double m_eta);
	~elasticForce();
	void computeFe(timeStepper &m_stepper);
    void computeFphi(timeStepper &m_stepper);
    void updateHistoryField();
    void computeJe(timeStepper &m_stepper);
    void setFirstJacobian(timeStepper &m_stepper);
    
    void setFirstJacobian_phi(timeStepper &m_stepper);
    
    void testFe();

    VectorXd reForce;

private:
	elasticPlate *plate;
    VectorXd lastEelPlus;


    double Gc_input;
    double ell_input;
    double eta_input;

    double Eel_plus;
    VectorXd fel_plus;
    MatrixXd Hel_plus;
    double Eel_minus;
    VectorXd fel_minus;
    MatrixXd Hel_minus; 

   void triElementElasticOnly(int idx, double &Eel_plus, VectorXd &fel_plus, MatrixXd &Hel_plus, double &Eel_minus, VectorXd &fel_minus, MatrixXd &Hel_minus);

    void triElementElasticTotal(
        int idx,
        double& Eel_total,
        VectorXd& fel_total,
        MatrixXd& Hel_total);

    void triElementElasticOnly(
        int idx,
        double& Eel,
        VectorXd& fel,
        MatrixXd& Hel);
    void phaseFieldSurfaceEnergyDerivatives(
        int idx,  
        const Vector3d& phi,              // [phi1, phi2, phi3]
        double Gc,
        double ell,
        double& Efrac,
        Vector3d& gradPhi,                // dE/dphi
        Matrix3d& hessPhi);                // d2E/dphi2
    double stvk(const Eigen::Matrix4d &C, 
        const Eigen::Matrix4d &T, 
        const Eigen::Vector4d &A_, const Eigen::Vector4d &B_, 
        const Eigen::Matrix<double, 4, 6> &gradA, 
        const Eigen::Matrix<double, 4, 6> &gradB, 
        const std::vector<Eigen::Matrix<double, 6, 6>> &hessA, 
        const std::vector<Eigen::Matrix<double, 6, 6>> &hessB, 
        Eigen::Matrix<double, 6, 1> *derivative, 
        Eigen::Matrix<double, 6, 6> *hessian);
};

#endif
