#include "elasticForce.h"
#include <iostream>

elasticForce::elasticForce(elasticPlate &m_plate,  double m_Gc, double m_ell, double m_eta)
{
	plate = &m_plate;
    lastEelPlus = VectorXd::Zero(plate->ne);


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

        triElementElasticAmor(i, Eel_plus, fel_plus, &Hel_plus, Eel_minus, fel_minus, &Hel_minus);
		lastEelPlus(i) = Eel_plus;

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

        // Hxx = g * Hel_plus + Hel_minus  (Amor split)
        hessian.block(0, 0, 6, 6) = g * Hel_plus + Hel_minus;


        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

        for (int j = 0; j < 6; j++)
        {
            m_stepper.addForce( arrayIndex(j), grad(j) );
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
void elasticForce::computeFeOnly(timeStepper &m_stepper)
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

        Eel_minus = 0.0;
        fel_minus.resize(6);
        fel_minus.setZero();

        triElementElasticAmor(i, Eel_plus, fel_plus, nullptr, Eel_minus, fel_minus, nullptr);
		lastEelPlus(i) = Eel_plus;

        double phiBar = (Phi(0) + Phi(1) + Phi(2)) / 3.0;
        double oneMinusPhiBar = 1.0 - phiBar;

        double g = oneMinusPhiBar * oneMinusPhiBar + eta_input;

        VectorXd grad;
        grad.resize(6);
        grad.setZero();

        // dE/dX = g * grad1
        grad.segment(0, 6) = g * fel_plus + fel_minus;

        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

        for (int j = 0; j < 6; j++)
        {
            m_stepper.addForce( arrayIndex(j), grad(j) );
        }

    }
}
void elasticForce::updateHistoryField()
{
    for (int i = 0; i < plate->ne; i++)
    {
        double H_max = 500.0 * (Gc_input / ell_input);
        double current_H = lastEelPlus(i) / plate->v_triangularElement[i].area;
        current_H = std::min(current_H, H_max);
        Vector3d Phi;
        Phi(0) = plate->getPhi(plate->v_triangularElement[i].nv_1);
        Phi(1) = plate->getPhi(plate->v_triangularElement[i].nv_2);
        Phi(2) = plate->getPhi(plate->v_triangularElement[i].nv_3);
        double phiBar = (Phi(0) + Phi(1) + Phi(2)) / 3.0;
        if (phiBar< 0.99)
        {
            plate->v_triangularElement[i].historyH = max(
			plate->v_triangularElement[i].historyH,current_H);
        }
        
    }
}

void elasticForce::computeFphi(timeStepper &m_stepper)
{
    for (int i = 0; i < plate->ne; i++)
    {
        Vector3d Phi;
        Phi(0) = plate->getPhi(plate->v_triangularElement[i].nv_1);
        Phi(1) = plate->getPhi(plate->v_triangularElement[i].nv_2);
        Phi(2) = plate->getPhi(plate->v_triangularElement[i].nv_3);

        double Ephi = 0.0;
        Vector3d fphi;
        Matrix3d Hphi;

        fphi.setZero();
        Hphi.setZero();

        phaseFieldSurfaceEnergyDerivatives(
            i,
            Phi,
            Gc_input,
            ell_input,
            Ephi,
            fphi,
            Hphi
        );

        const double area = plate->v_triangularElement[i].area;
        const double H = plate->v_triangularElement[i].historyH;

        const double phiBar = Phi.mean();
        const double oneMinusPhiBar = 1.0 - phiBar;

        Vector3d dg;
        dg.setConstant(-2.0 * oneMinusPhiBar / 3.0);

        Matrix3d d2g;
        d2g.setConstant(2.0 / 9.0);

        Vector3d grad = H * area * dg + fphi;
        Matrix3d hessian = H * area * d2g + Hphi;

        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex_phi;

        for (int j = 0; j < 3; j++)
        {
            m_stepper.addForce(arrayIndex(j), grad(j));
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
void elasticForce::setFirstJacobian_phi(timeStepper &m_stepper)
{
	for (int i = 0; i < plate->ne; i++)
	{
		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex_phi;

		for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                m_stepper.addFirstJacobian(arrayIndex(j), arrayIndex(k));
            }
        }
	}
}

void elasticForce::triElementElasticMiehe(
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
		Matrix2d S_mat = M * abarinv *(A - abar) * Minv;
        Matrix2d S_sym = 0.5 * (S_mat + S_mat.transpose());

        SelfAdjointEigenSolver<Matrix2d> solver(S_sym);

        if (solver.info() != Success) {
            cerr << "特征值分解失败！" << endl;
            exit(1);
        }
        Vector2d lambdas = solver.eigenvalues(); 
        Matrix2d Q = solver.eigenvectors();

        

        // 
        Vector2d lambdas_plus,lambdas_minus;
        lambdas_plus(0) = max(lambdas(0), 0.0);
        lambdas_plus(1) = max(lambdas(1), 0.0);

        lambdas_minus(0) = min(lambdas(0), 0.0);
        lambdas_minus(1) = min(lambdas(1), 0.0);

        // 
        

        
        Vector2d I_plus,I_minus;
        I_plus(0) = (lambdas(0) > 0.0) ? 1.0 : 0.0;
        I_plus(1) = (lambdas(1) > 0.0) ? 1.0 : 0.0;
        I_minus(0) = (lambdas(0) < 0.0) ? 1.0 : 0.0;
        I_minus(1) = (lambdas(1) < 0.0) ? 1.0 : 0.0;

        Matrix2d P1 = Q * I_plus.asDiagonal() * Q.transpose();
        Matrix2d P2 = Q * I_minus.asDiagonal() * Q.transpose();

        Matrix4d kronc_pl = Geometry::kron( Minv.transpose(), P1 * M * abarinv); 
        Matrix4d kronc_mi = Geometry::kron( Minv.transpose(), P2 * M * abarinv);
        Matrix4d kronc = Geometry::kron( Minv.transpose(),  M * abarinv);
		

		// // stvk energy/grad/hessian  for epsilon inner product <eps,eps> with C as the stiffness matrix
        // lambda energy
        Matrix<double, 6, 1> derivative1;
        Matrix<double, 6, 6> hessian1;
        double energy1 = stvk(plate->v_triangularElement[idx].C1, kronc, eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess, &derivative1, &hessian1);
       

        // mu energy for plus part
        Matrix<double, 6, 1> derivative2;
        Matrix<double, 6, 6> hessian2;
        double energy2 = stvk(plate->v_triangularElement[idx].C2, kronc_pl, eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess, &derivative2, &hessian2);
        
        // mu energy for minus part
        Matrix<double, 6, 1> derivative3;
        Matrix<double, 6, 6> hessian3;        
        double energy3 = stvk(plate->v_triangularElement[idx].C2, kronc_mi, eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess, &derivative3, &hessian3);
        

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

void elasticForce::triElementElasticAmor(
        int idx,
        double& Eel_plus, VectorXd& fel_plus, MatrixXd* Hel_plus,
        double& Eel_minus, VectorXd& fel_minus, MatrixXd* Hel_minus)
    {
        const bool needHess = (Hel_plus != nullptr) || (Hel_minus != nullptr);

        Matrix2d A;
        Matrix<double, 4, 6> fFFderiv;
        std::vector<Matrix<double, 6, 6> > fFFhess;
        Vector2d x1 = plate->getVertex(plate->v_triangularElement[idx].nv_1);
        Vector2d x2 = plate->getVertex(plate->v_triangularElement[idx].nv_2);
        Vector2d x3 = plate->getVertex(plate->v_triangularElement[idx].nv_3);

        A = Geometry::firstFundamentalForm( x1 , x2 , x3 , &fFFderiv,
                                            needHess ? &fFFhess : nullptr );

        Matrix2d M     = plate->v_triangularElement[idx].M;
        Matrix2d Minv  = M.inverse();
        Matrix2d abar    = plate->v_triangularElement[idx].abar;
        Matrix2d abarinv = plate->v_triangularElement[idx].abarinv;

        double dA = plate->v_triangularElement[idx].area;

        Matrix2d s = A - abar;
        Vector4d eps = Eigen::Map<Vector4d>(s.data());
        Matrix2d S_mat = M * abarinv * s * Minv;
        double trS = S_mat.trace();

        // Full (unsplit) projection
        Matrix4d kronc = Geometry::kron( Minv.transpose(), M * abarinv);

        // --- Volumetric energy  ---
        Matrix<double, 6, 1> f_vol;
        Matrix<double, 6, 6> H_vol;
        double E_vol = stvk(plate->v_triangularElement[idx].C_vol, kronc,
                           eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess,
                           &f_vol, needHess ? &H_vol : nullptr);

        // --- Deviatoric energy ---
        Matrix<double, 6, 1> f_dev;
        Matrix<double, 6, 6> H_dev;
        double E_dev = stvk(plate->v_triangularElement[idx].C_dev, kronc,
                            eps, eps, fFFderiv, fFFderiv, fFFhess, fFFhess,
                            &f_dev, needHess ? &H_dev : nullptr);

        double coeff = 1.0 / 4.0;

        if (trS > 0.0)
        {
            // Expansion: all energy degrades
            Eel_plus  = coeff * dA * (E_vol + E_dev);
            fel_plus  = coeff * dA * (f_vol + f_dev);
            if (Hel_plus) *Hel_plus  = coeff * dA * (H_vol + H_dev);
            Eel_minus = 0.0;
            fel_minus.setZero(6);
            if (Hel_minus) Hel_minus->setZero(6, 6);
        }
        else
        {
            // Compression: deviatoric degrades, volumetric resists
            Eel_plus = coeff * dA * E_dev;
            fel_plus = coeff * dA * f_dev;
            if (Hel_plus) *Hel_plus = coeff * dA * H_dev;

            Eel_minus = coeff * dA * E_vol;
            fel_minus = coeff * dA * f_vol;
            if (Hel_minus) *Hel_minus = coeff * dA * H_vol;
        }
    }

void elasticForce::triElementElastic(
        int idx,
        double& Eel_total,
        VectorXd& fel_total,
        MatrixXd& Hel_total)
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

        double dA = plate->v_triangularElement[idx].area;

        Matrix2d s = A - abar;
        Vector4d eps = Eigen::Map<Vector4d>(s.data());

        // Full (unsplit) projection: no spectral decomposition
        Matrix4d kronc = Geometry::kron( Minv.transpose(), M * abarinv);

        // Total material stiffness = C1 (lambda) + C2 (mu)
        Matrix4d C_total = plate->v_triangularElement[idx].C1
                         + plate->v_triangularElement[idx].C2;

        Matrix<double, 6, 1> derivative;
        Matrix<double, 6, 6> hessian;
        double energy = stvk(C_total, kronc, eps, eps,
                             fFFderiv, fFFderiv, fFFhess, fFFhess,
                             &derivative, &hessian);

        double coeff = 1.0 / 4.0;

        Eel_total = coeff * dA * energy;
        fel_total = coeff * dA * derivative;
        Hel_total = coeff * dA * hessian;
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