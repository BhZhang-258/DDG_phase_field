#include "elasticForce.h"
#include <iostream>

elasticForce::elasticForce(elasticPlate &m_plate, timeStepper &m_stepper, double m_Gc, double m_ell, double m_eta)
{
	plate = &m_plate;
	stepper = &m_stepper;

    alpha_input = sqrt(plate->YoungM / (1 - plate->Possion));
    beta_input = sqrt(plate->YoungM / (1 + plate->Possion));

    Gc_input = m_Gc;
    ell_input = m_ell;
    eta_input = m_eta;
}

elasticForce::~elasticForce()
{
    ;
}

void elasticForce::computeFe()
{
    reForce = VectorXd::Zero(plate->ndof);

    for (int i = 0; i < plate->ne; i++)
    {
        MatrixXd undeform = plate->v_triangularElement[i].X;
        MatrixXd deform;
        deform = MatrixXd::Zero(2,3);

        deform.col(0) = plate->getVertex(plate->v_triangularElement[i].nv_1);
        deform.col(1) = plate->getVertex(plate->v_triangularElement[i].nv_2);
        deform.col(2) = plate->getVertex(plate->v_triangularElement[i].nv_3);



        Vector3d Phi;
        Phi(0) = plate->getPhi(plate->v_triangularElement[i].nv_1);
        Phi(1) = plate->getPhi(plate->v_triangularElement[i].nv_2);
        Phi(2) = plate->getPhi(plate->v_triangularElement[i].nv_3);



        double E1;
        VectorXd grad1;
        MatrixXd hession1;
        triElementElasticOnly(deform, undeform, alpha_input, beta_input, E1, grad1, hession1);

    
        double E2;
        Vector3d grad2;
        Matrix3d hession2;
        phaseFieldSurfaceEnergyDerivatives(undeform, Phi, Gc_input, ell_input, E2, grad2, hession2);


        double phiBar = (Phi(0) + Phi(1) + Phi(2)) / 3.0;
        double oneMinusPhiBar = 1.0 - phiBar;

        double g = oneMinusPhiBar * oneMinusPhiBar + eta_input;

        // dg/dPhi
        Vector3d dg;
        dg.setConstant(-2.0 * oneMinusPhiBar / 3.0);

        // d2g/dPhi2
        Matrix3d d2g;
        d2g.setConstant(2.0 / 9.0); 

        double totalEnergy = E1 * g + E2;

        VectorXd grad;
        grad.resize(9);
        grad.setZero();

        // dE/dX = g * grad1
        grad.segment(0, 6) = g * grad1;

        // dE/dPhi = E1 * dg + grad2
        grad.segment(6, 3) = E1 * dg + grad2;


        MatrixXd hession;
        hession.resize(9, 9);
        hession.setZero();

        // Hxx = g * hession1
        hession.block(0, 0, 6, 6) = g * hession1;

        // Hxphi = grad1 * dg^T   (6x3)
        hession.block(0, 6, 6, 3) = grad1 * dg.transpose();

        // Hphix = dg * grad1^T   (3x6)
        hession.block(6, 0, 3, 6) = dg * grad1.transpose();

        // Hphiphi = E1 * d2g + hession2
        hession.block(6, 6, 3, 3) = E1 * d2g + hession2;



        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

        for (int j = 0; j < 9; j++)
        {
            stepper->addForce( arrayIndex(j), grad(j) );

            reForce(arrayIndex(j)) = reForce(arrayIndex(j)) + grad(j);
        }

        for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), hession(j,k) );
            }
        }
        

    }
}

void elasticForce::computeJe()
{
	;
}

void elasticForce::setFirstJacobian()
{
	for (int i = 0; i < plate->ne; i++)
	{
		VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), 1);
            }
        }
	}
}

double elasticForce::cotangent2D(const Vector2d& a, const Vector2d& b)
    {
        const double cross = a.x() * b.y() - a.y() * b.x();
        return a.dot(b) / cross;
    }


void elasticForce::closestRotationAndY2D(
        const Matrix2d& F,
        Matrix2d& R,
        Matrix2d& Y)
    {
        Eigen::JacobiSVD<Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Matrix2d U = svd.matrixU();
        Matrix2d V = svd.matrixV();

        Matrix2d C = Matrix2d::Identity();
        C(1, 1) = (U * V.transpose()).determinant();

        R = U * C * V.transpose();
        Y = R.transpose() * F;
        Y = 0.5 * (Y + Y.transpose());
    }

void elasticForce::setBlock(MatrixXd& H, int bi, int bj, const Matrix2d& A)
    {
        H.block<2, 2>(2 * bi, 2 * bj) = A;
    }


Matrix2d elasticForce::secondPartMat(
        const Vector2d& si,
        const Vector2d& sj,
        double A0,
        double trY)
    {
        return (si * sj.transpose()) / (A0 * trY);
    }


Matrix2d elasticForce::abSplitMat(
        const Vector2d& ri,
        const Vector2d& rj,
        double A0)
    {
        return (ri * rj.transpose()) / A0;
    }

void elasticForce::triElementElasticOnly(
        const Eigen::Matrix<double, 2, 3>& x,
        const Eigen::Matrix<double, 2, 3>& X,
        double alpha,
        double beta,
        double& Eel,
        VectorXd& fel,
        MatrixXd& Hel)
    {
        using std::abs;

        const double alphaSq = alpha * alpha;
        const double betaSq  = beta  * beta;

        const Matrix2d I2 = Matrix2d::Identity();
        Matrix2d J;
        J << 0.0, -1.0,
             1.0,  0.0;

        // Reference and current vertices
        const Vector2d X1 = X.col(0);
        const Vector2d X2 = X.col(1);
        const Vector2d X3 = X.col(2);

        const Vector2d x1 = x.col(0);
        const Vector2d x2 = x.col(1);
        const Vector2d x3 = x.col(2);

        // Reference shape matrix and area
        Matrix2d Dm;
        Dm.col(0) = X2 - X1;
        Dm.col(1) = X3 - X1;

        const double detDm = Dm.determinant();
        if (detDm <= 1e-14)
        {
            throw std::runtime_error("triElementElasticOnly: reference triangle is degenerate or not CCW.");
        }

        const double A0 = 0.5 * detDm;

        // Deformation gradient
        Matrix2d Ds;
        Ds.col(0) = x2 - x1;
        Ds.col(1) = x3 - x1;

        const Matrix2d F = Ds * Dm.inverse();

        // Closest rotation and symmetric stretch
        Matrix2d R, Y;
        closestRotationAndY2D(F, R, Y);

        // Energy
        if (std::abs(alpha - beta) < 1e-14)
        {
            const Matrix2d dF = F - R;
            Eel = 0.5 * alphaSq * A0 * dF.squaredNorm();
        }
        else
        {
            const double trY = Y.trace();
            const Matrix2d devY = Y - 0.5 * trY * I2;
            const double vol = 0.5 * trY - 1.0;
            Eel = 0.5 * A0 * (betaSq * devY.squaredNorm() + 2.0 * alphaSq * vol * vol);
        }

        // Cotangents of reference triangle
        const double cot1 = cotangent2D(X2 - X1, X3 - X1);
        const double cot2 = cotangent2D(X3 - X2, X1 - X2);
        const double cot3 = cotangent2D(X1 - X3, X2 - X3);

        // Area gradients in reference configuration
        const Vector2d a1  = 0.5 * J * (X3 - X2);
        const Vector2d a2v = 0.5 * J * (X1 - X3);
        const Vector2d a3  = 0.5 * J * (X2 - X1);

        // Laplacian part
        const Matrix2d L11 = 0.5 * (cot2 + cot3) * I2;
        const Matrix2d L22 = 0.5 * (cot1 + cot3) * I2;
        const Matrix2d L33 = 0.5 * (cot1 + cot2) * I2;

        const Matrix2d L12 = -0.5 * cot3 * I2;
        const Matrix2d L13 = -0.5 * cot2 * I2;
        const Matrix2d L23 = -0.5 * cot1 * I2;

        const Matrix2d L21 = L12;
        const Matrix2d L31 = L13;
        const Matrix2d L32 = L23;

        // Gradient of elastic energy
        const double trY = Y.trace();
        const double coeff = alphaSq - 0.5 * (alphaSq - betaSq) * trY;

        const Vector2d g1 = betaSq * (L11 * x1 + L12 * x2 + L13 * x3) - coeff * (R * a1);
        const Vector2d g2 = betaSq * (L21 * x1 + L22 * x2 + L23 * x3) - coeff * (R * a2v);
        const Vector2d g3 = betaSq * (L31 * x1 + L32 * x2 + L33 * x3) - coeff * (R * a3);

        // Gradient of elastic energy
        fel.resize(6);
        fel.segment<2>(0) = g1;
        fel.segment<2>(2) = g2;
        fel.segment<2>(4) = g3;

        // Hessian / tangent of elastic part
        double trY_safe = trY;
        if (abs(trY_safe) < 1e-12)
        {
            trY_safe = (trY_safe >= 0.0 ? 1e-12 : -1e-12);
        }

        const Vector2d s1 = R * J * a1;
        const Vector2d s2 = R * J * a2v;
        const Vector2d s3 = R * J * a3;

        const Vector2d r1 = R * a1;
        const Vector2d r2 = R * a2v;
        const Vector2d r3 = R * a3;

        Matrix2d HE11, HE12, HE13;
        Matrix2d HE21, HE22, HE23;
        Matrix2d HE31, HE32, HE33;

        if (std::abs(alpha - beta) < 1e-14)
        {
            HE11 = alphaSq * (L11 - secondPartMat(s1, s1, A0, trY_safe));
            HE12 = alphaSq * (L12 - secondPartMat(s1, s2, A0, trY_safe));
            HE13 = alphaSq * (L13 - secondPartMat(s1, s3, A0, trY_safe));

            HE21 = alphaSq * (L21 - secondPartMat(s2, s1, A0, trY_safe));
            HE22 = alphaSq * (L22 - secondPartMat(s2, s2, A0, trY_safe));
            HE23 = alphaSq * (L23 - secondPartMat(s2, s3, A0, trY_safe));

            HE31 = alphaSq * (L31 - secondPartMat(s3, s1, A0, trY_safe));
            HE32 = alphaSq * (L32 - secondPartMat(s3, s2, A0, trY_safe));
            HE33 = alphaSq * (L33 - secondPartMat(s3, s3, A0, trY_safe));
        }
        else
        {
            const double coeffSecond = alphaSq - 0.5 * (alphaSq - betaSq) * trY;
            const double coeffAB     = 0.5 * (alphaSq - betaSq);

            HE11 = betaSq * L11 - coeffSecond * secondPartMat(s1, s1, A0, trY_safe) + coeffAB * abSplitMat(r1, r1, A0);
            HE12 = betaSq * L12 - coeffSecond * secondPartMat(s1, s2, A0, trY_safe) + coeffAB * abSplitMat(r1, r2, A0);
            HE13 = betaSq * L13 - coeffSecond * secondPartMat(s1, s3, A0, trY_safe) + coeffAB * abSplitMat(r1, r3, A0);

            HE21 = betaSq * L21 - coeffSecond * secondPartMat(s2, s1, A0, trY_safe) + coeffAB * abSplitMat(r2, r1, A0);
            HE22 = betaSq * L22 - coeffSecond * secondPartMat(s2, s2, A0, trY_safe) + coeffAB * abSplitMat(r2, r2, A0);
            HE23 = betaSq * L23 - coeffSecond * secondPartMat(s2, s3, A0, trY_safe) + coeffAB * abSplitMat(r2, r3, A0);

            HE31 = betaSq * L31 - coeffSecond * secondPartMat(s3, s1, A0, trY_safe) + coeffAB * abSplitMat(r3, r1, A0);
            HE32 = betaSq * L32 - coeffSecond * secondPartMat(s3, s2, A0, trY_safe) + coeffAB * abSplitMat(r3, r2, A0);
            HE33 = betaSq * L33 - coeffSecond * secondPartMat(s3, s3, A0, trY_safe) + coeffAB * abSplitMat(r3, r3, A0);
        }

        Hel.resize(6, 6);
        Hel.setZero();

        setBlock(Hel, 0, 0, HE11);
        setBlock(Hel, 0, 1, HE12);
        setBlock(Hel, 0, 2, HE13);

        setBlock(Hel, 1, 0, HE21);
        setBlock(Hel, 1, 1, HE22);
        setBlock(Hel, 1, 2, HE23);

        setBlock(Hel, 2, 0, HE31);
        setBlock(Hel, 2, 1, HE32);
        setBlock(Hel, 2, 2, HE33);
    }

void elasticForce::phaseFieldSurfaceEnergyDerivatives(
    const MatrixXd& referenceShape,   // 2x3
    const Vector3d& phi,              // [phi1, phi2, phi3]
    double Gc,
    double ell,
    double& Efrac,
    Vector3d& gradPhi,                // dE/dphi
    Matrix3d& hessPhi)                // d2E/dphi2
{
    if (referenceShape.rows() != 2 || referenceShape.cols() != 3)
    {
        throw std::runtime_error("referenceShape must be 2x3.");
    }

    const Vector2d X1 = referenceShape.col(0);
    const Vector2d X2 = referenceShape.col(1);
    const Vector2d X3 = referenceShape.col(2);

    // Reference triangle area
    Matrix2d Dm;
    Dm.col(0) = X2 - X1;
    Dm.col(1) = X3 - X1;

    const double detDm = Dm.determinant();
    if (detDm <= 1e-14)
    {
        throw std::runtime_error("Degenerate or non-CCW reference triangle.");
    }

    const double A0 = 0.5 * detDm;

    // --------------------------
    // Mphi
    // --------------------------
    Matrix3d Mphi;
    Mphi << 2.0, 1.0, 1.0,
            1.0, 2.0, 1.0,
            1.0, 1.0, 2.0;
    Mphi *= (A0 / 12.0);

    // --------------------------
    // Bphi
    // --------------------------
    Vector2d gradN1, gradN2, gradN3;
    gradN1 << (X2.y() - X3.y()), (X3.x() - X2.x());
    gradN2 << (X3.y() - X1.y()), (X1.x() - X3.x());
    gradN3 << (X1.y() - X2.y()), (X2.x() - X1.x());

    gradN1 /= (2.0 * A0);
    gradN2 /= (2.0 * A0);
    gradN3 /= (2.0 * A0);

    Eigen::Matrix<double, 2, 3> Bphi;
    Bphi.col(0) = gradN1;
    Bphi.col(1) = gradN2;
    Bphi.col(2) = gradN3;

    // --------------------------
    // Kphi
    // --------------------------
    Matrix3d Kphi = A0 * (Bphi.transpose() * Bphi);

    // --------------------------
    // Aphi = (Gc/ell) Mphi + (Gc*ell) Kphi
    // so that E = 0.5 * phi^T Aphi phi
    // --------------------------
    hessPhi = (Gc / ell) * Mphi + (Gc * ell) * Kphi;

    // Energy
    Efrac = 0.5 * phi.dot(hessPhi * phi);

    // Gradient
    gradPhi = hessPhi * phi;
}