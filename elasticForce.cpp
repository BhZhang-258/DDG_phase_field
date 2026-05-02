#include "elasticForce.h"
#include <iostream>

elasticForce::elasticForce(elasticPlate &m_plate,  double m_Gc, double m_ell, double m_eta)
{
	plate = &m_plate;


    Gc_input = m_Gc;
    ell_input = m_ell;
    eta_input = m_eta;
}

elasticForce::~elasticForce()
{
    ;
}

void elasticForce::computeFe(timeStepper &m_stepper)
{
    
    for (int i = 0; i < plate->ne; i++)
    {

        Vector3d Phi;
        Phi(0) = plate->getPhi(plate->v_triangularElement[i].nv_1);
        Phi(1) = plate->getPhi(plate->v_triangularElement[i].nv_2);
        Phi(2) = plate->getPhi(plate->v_triangularElement[i].nv_3);        

        Eel_plus = 0.0;
        fel_plus.resize(6);
        fel_plus.setZero();
        Hel_plus.resize(6, 6);
        Hel_plus.setZero();

        Eel_minus = 0.0;
        fel_minus.resize(6);
        fel_minus.setZero();
        Hel_minus.resize(6, 6);
        Hel_minus.setZero();

        triElementElasticOnly(i, Eel_plus, fel_plus, Hel_plus, Eel_minus, fel_minus, Hel_minus);

        double phiBar = (Phi(0) + Phi(1) + Phi(2)) / 3.0;
        double oneMinusPhiBar = 1.0 - phiBar;

        double g = oneMinusPhiBar * oneMinusPhiBar + eta_input;

        double totalEnergy = Eel_plus * g + Eel_minus;

        VectorXd grad;
        grad.resize(6);
        grad.setZero();

        // dE/dX = g * grad1
        grad.segment(0, 6) = g * fel_plus + fel_minus;
        
        MatrixXd hessian;
        hessian.resize(6, 6);
        hessian.setZero();

        // Hxx = g * hessian1
        hessian.block(0, 0, 6, 6) = g * Hel_plus + Hel_minus;

        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

        for (int j = 0; j < 6; j++)
        {
            m_stepper.addForce( arrayIndex(j), grad(j) );

            reForce(arrayIndex(j)) = reForce(arrayIndex(j)) + grad(j);
        }

        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                m_stepper.addJacobian(arrayIndex(j), arrayIndex(k), hessian(j,k) );
            }
        }
    

    }
}
void elasticForce::computeFphi( timeStepper &m_stepper)
{
    

    for (int i = 0; i < plate->ne; i++)
    {

        // H历史变量在哪里
        Vector3d Phi;
        Phi(0) = plate->getPhi(plate->v_triangularElement[i].nv_1);
        Phi(1) = plate->getPhi(plate->v_triangularElement[i].nv_2);
        Phi(2) = plate->getPhi(plate->v_triangularElement[i].nv_3);        


        double Ephi;
        Vector3d fphi;
        Matrix3d Hphi;
        phaseFieldSurfaceEnergyDerivatives(i, Phi, Gc_input, ell_input, Ephi, fphi, Hphi);


        double phiBar = (Phi(0) + Phi(1) + Phi(2)) / 3.0;
        double oneMinusPhiBar = 1.0 - phiBar;

        double g = oneMinusPhiBar * oneMinusPhiBar + eta_input;
        
        // dg/dPhi
        Vector3d dg;
        dg.setConstant(-2.0 * oneMinusPhiBar / 3.0);

        // d2g/dPhi2
        Matrix3d d2g;
        d2g.setConstant(2.0 / 9.0); 

        double totalEnergy = Eel_plus * g + Eel_minus + Ephi;
        double scale = ell_input / Gc_input;


        VectorXd grad;
        grad.resize(3);
        grad.setZero();

        // dE/dPhi = E1 * dg + grad2
        grad.segment(0, 3) =  Eel_plus * dg + fphi;


        MatrixXd hessian;
        hessian.resize(3, 3);
        hessian.setZero();

        // Hxx = g * hessian1
        hessian.block(0, 0, 3, 3) = Eel_plus* d2g + Hphi;


        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex_phi;

        for (int j = 0; j < 3; j++)
        {
            m_stepper.addForce( arrayIndex(j), grad(j) );

            reForce(arrayIndex(j)) = reForce(arrayIndex(j)) + grad(j);
        }

        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                m_stepper.addJacobian(arrayIndex(j), arrayIndex(k), hessian(j,k) );
            }
        }
        

    }
}
void elasticForce::computeJe(timeStepper &m_stepper)
{
	;
}

void elasticForce::setFirstJacobian(timeStepper &m_stepper)
{
	for (int i = 0; i < plate->ne; i++)
	{
		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

		for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                m_stepper.addFirstJacobian(arrayIndex(j), arrayIndex(k));
            }
        }
	}
}

void elasticForce::triElementElasticOnly(
        int idx,
        double& Eel_plus,
        VectorXd& fel_plus,
        MatrixXd& Hel_plus,
        double& Eel_minus,
        VectorXd& fel_minus,
        MatrixXd& Hel_minus)
    {
        Matrix2d A;
        Matrix<double, 4, 6> fFFderiv;
        std::vector<Matrix<double, 6, 6> > fFFhess;
        Vector2d x1 = plate->getVertex(plate->v_triangularElement[idx].nv_1);
        Vector2d x2 = plate->getVertex(plate->v_triangularElement[idx].nv_2);
        Vector2d x3 = plate->getVertex(plate->v_triangularElement[idx].nv_3);

        
		A = Geometry::firstFundamentalForm( x1 , x2 , x3 , &fFFderiv,  &fFFhess );
        
        
		Matrix2d M = plate->v_triangularElement[idx].M;
		Matrix2d Minv = M.inverse();
		
		Matrix2d abar = plate->v_triangularElement[idx].abar;
		Matrix2d abarinv = plate->v_triangularElement[idx].abarinv;

		Matrix4d C = plate->v_triangularElement[idx].C;

		
		
    	double dA = plate->v_triangularElement[idx].area;
	
		Matrix<double, 6, 1> derivative;
		Matrix<double, 6, 6> hessian;
		
		Matrix2d s =  A - abar;
		Vector4d eps = Eigen::Map<Vector4d>(s.data());
        // Matrix2d S_edge = abarinv *(A - abar);
		Matrix2d S_mat = M * abarinv *(A - abar) * Minv;
        Matrix2d S_sym = 0.5 * (S_mat + S_mat.transpose());

        SelfAdjointEigenSolver<Matrix2d> solver(S_sym);
        // cout << "eig of S_mat:\n" << solver.eigenvalues() << endl;

        if (solver.info() != Success) {
            cerr << "特征值分解失败！" << endl;
            exit(1);
        }
        Vector2d lambdas = solver.eigenvalues(); 
        Matrix2d Q = solver.eigenvectors();      

        // cout << "特征值 L (向量形式):\n" << lambdas << "\n\n";
        // cout << "特征向量 Q:\n" << Q << "\n\n";

        // 
        Vector2d lambdas_plus,lambdas_minus;
        lambdas_plus(0) = max(lambdas(0), 0.0);
        lambdas_plus(1) = max(lambdas(1), 0.0);

        lambdas_minus(0) = min(lambdas(0), 0.0);
        lambdas_minus(1) = min(lambdas(1), 0.0);

        // 
        // Matrix2d S_mat_plus = Q * lambdas_plus.asDiagonal() * Q.transpose();
        // Matrix2d S_mat_minus = Q * lambdas_minus.asDiagonal() * Q.transpose();

        // cout << "拉伸部分 S_mat_plus:\n" << S_mat_plus << endl;
        // cout << "压缩部分 S_mat_minus:\n" << S_mat_minus << endl;

        
        Vector2d I_plus,I_minus;
        I_plus(0) = (lambdas(0) > 0.0) ? 1.0 : 0.0;
        I_plus(1) = (lambdas(1) > 0.0) ? 1.0 : 0.0;
        I_minus(0) = (lambdas(0) < 0.0) ? 1.0 : 0.0;
        I_minus(1) = (lambdas(1) < 0.0) ? 1.0 : 0.0;

        Matrix2d P1 = Q * I_plus.asDiagonal() * Q.transpose();
        Matrix2d P2 = Q * I_minus.asDiagonal() * Q.transpose();

        Matrix4d kronc_pl = Geometry::kron( Minv.transpose(), P1 * M * abarinv); 
        Matrix4d kronc_mi = Geometry::kron( Minv.transpose(), P2 * M * abarinv);
        Matrix4d kronc = Geometry::kron( Minv.transpose(), M * abarinv);
		

		// // stvk energy/grad/hessian  for epsilon inner product <eps,eps> with C as the stiffness matrix
        // lambda energy
        Matrix<double, 6, 1> derivative1;
        Matrix<double, 6, 6> hessian1;
        double energy1 = stvk(plate->v_triangularElement[idx].C1, kronc, eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess, &derivative1, &hessian1);
        // cout << "\nplus part energy:"  << std::endl;
        // cout << "energy:\n" << energy1 << std::endl;
        // cout << "derivative:\n" << derivative1.transpose() << std::endl;
        // cout << "hessian:\n" << hessian1 << std::endl;

        // mu energy for plus part
        Matrix<double, 6, 1> derivative2;
        Matrix<double, 6, 6> hessian2;
        double energy2 = stvk(plate->v_triangularElement[idx].C2, kronc_pl, eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess, &derivative2, &hessian2);
        // cout << "\nminus part energy:"  << std::endl;
        // cout << "energy:\n" << energy2 << std::endl;
        // cout << "derivative:\n" << derivative2.transpose() << std::endl;
        // cout << "hessian:\n" << hessian2 << std::endl;
        // mu energy for minus part
        Matrix<double, 6, 1> derivative3;
        Matrix<double, 6, 6> hessian3;        
        double energy3 = stvk(plate->v_triangularElement[idx].C2, kronc_mi, eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess, &derivative3, &hessian3);
        // cout << "\nminus part energy:"  << std::endl;
        // cout << "energy:\n" << energy3 << std::endl;
        // cout << "derivative:\n" << derivative3.transpose() << std::endl;
        // cout << "hessian:\n" << hessian3 << std::endl;

        double coeff = 1.0/4.0 * 1;  // plane stress/strain coefficient

        Eel_plus = coeff * dA * energy2;
        fel_plus = coeff * dA * derivative2;
        Hel_plus = coeff * dA * hessian2;

        Eel_minus = coeff * dA * energy3;
        fel_minus = coeff * dA * derivative3;
        Hel_minus = coeff * dA * hessian3;

        if (S_mat.trace() > 0 )
        {
            Eel_plus += coeff * dA * energy1;
            fel_plus += coeff * dA * derivative1;
            Hel_plus += coeff * dA * hessian1;

        }
        else 
        {

            Eel_minus += coeff * dA * energy1;
            fel_minus += coeff * dA * derivative1;
            Hel_minus += coeff * dA * hessian1;

        }
        
        
    }

void elasticForce::phaseFieldSurfaceEnergyDerivatives(
    int idx,
    const Vector3d& phi,              // [phi1, phi2, phi3]
    double Gc,
    double ell,
    double& Efrac,
    Vector3d& gradPhi,                // dE/dphi
    Matrix3d& hessPhi)                // d2E/dphi2
{
   
    const double A0 = plate->v_triangularElement[idx].area;
    

    hessPhi = (Gc / ell) * plate->v_triangularElement[idx].Mphi + (Gc * ell) * plate->v_triangularElement[idx].Kphi;

    // Energy
    Efrac = 0.5 * phi.dot(hessPhi * phi);

    // Gradient
    gradPhi = hessPhi * phi;
}

double elasticForce::stvk(
    const Eigen::Matrix4d &C,
    const Eigen::Matrix4d &T,
    const Eigen::Vector4d &A_,
    const Eigen::Vector4d &B_,
    const Eigen::Matrix<double,4,6> &gradA,
    const Eigen::Matrix<double,4,6> &gradB,
    const std::vector <Eigen::Matrix<double, 6, 6> > &hessA,
    const std::vector <Eigen::Matrix<double, 6, 6> > &hessB,
    Eigen::Matrix<double,6,1>* derivative,
    Eigen::Matrix<double,6,6>* hessian)
{

    Eigen::Vector4d eps1 = T * A_;
    Eigen::Vector4d eps2 = T * B_;

    double energy = 0.5 * eps2.transpose() * C * eps1;

    Eigen::MatrixXd B1 = T * gradA;
    Eigen::MatrixXd B2 = T * gradB;

    if (derivative)
    {
        *derivative =
            0.5 * B1.transpose() * C.transpose() * eps2 +
            0.5 * B2.transpose() * C.transpose() * eps1;
    }

    if (hessian)
    {
        Eigen::Matrix<double,6,6> Hmat =
            0.5 * B1.transpose() * C * B2 +
            0.5 * B2.transpose() * C * B1;

        Eigen::Matrix<double,6,6> Hgeo =
            Eigen::Matrix<double,6,6>::Zero();

        Eigen::RowVectorXd coeff1 = 0.5 * eps2.transpose() * C * T;
        Eigen::RowVectorXd coeff2 = 0.5 * eps1.transpose() * C * T;

        for(int i=0;i<4;i++)
        {
            Hgeo += coeff1(i) * hessA[i];
            Hgeo += coeff2(i) * hessB[i];
        }

        *hessian = Hmat + Hgeo;
    }

    return energy;
}