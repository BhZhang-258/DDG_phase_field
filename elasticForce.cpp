#include "elasticForce.h"
#include <iostream>

elasticForce::elasticForce(elasticPlate &m_plate, timeStepper &m_stepper, 
    double m_alpha, double m_beta)
{
	plate = &m_plate;
	stepper = &m_stepper;

    alpha = m_alpha * plate->YoungM;
    beta = m_beta * plate->YoungM;
}

elasticForce::~elasticForce()
{
    ;
}

void elasticForce::computeFe()
{
    for (int i = 0; i < plate->ne; i++)
    {
        MatrixXd undeform = plate->v_triangularElement[i].X;
        MatrixXd deform;
        deform = MatrixXd::Zero(2,3);

        deform.col(0) = plate->getVertex(plate->v_triangularElement[i].nv_1);
        deform.col(1) = plate->getVertex(plate->v_triangularElement[i].nv_2);
        deform.col(2) = plate->getVertex(plate->v_triangularElement[i].nv_3);
        
        double localEnergy;
        VectorXd force;
        MatrixXd jacob;

        triElement(deform, undeform, alpha, beta, localEnergy, force, jacob);

        //cout << force.transpose() << endl;

        VectorXi arrayIndex = plate->v_triangularElement[i].arrayIndex;

        for (int j = 0; j < 6; j++)
        {
            stepper->addForce( arrayIndex(j), -force(j) );
        }

        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), -jacob(j,k) );
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

		for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
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

void elasticForce::triElement(
    const Eigen::Matrix<double, 2, 3>& x,
    const Eigen::Matrix<double, 2, 3>& X,
    double materialParameter1,
    double materialParameter2,
    double& E,
    VectorXd& f,
    MatrixXd& H)
{
    using std::abs;

    const double alpha = materialParameter1;
    const double beta  = materialParameter2;
    const double a2 = alpha * alpha;
    const double b2 = beta  * beta;

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
        throw std::runtime_error("triElement: reference triangle is degenerate or not CCW.");
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
        E = 0.5 * a2 * A0 * dF.squaredNorm();
    }
    else
    {
        const double trY = Y.trace();
        const Matrix2d devY = Y - 0.5 * trY * I2;
        const double vol = 0.5 * trY - 1.0;
        E = 0.5 * A0 * (b2 * devY.squaredNorm() + 2.0 * a2 * vol * vol);
    }

    // Cotangents of reference triangle
    const double cot1 = cotangent2D(X2 - X1, X3 - X1);
    const double cot2 = cotangent2D(X3 - X2, X1 - X2);
    const double cot3 = cotangent2D(X1 - X3, X2 - X3);

    // Area gradients in reference configuration
    const Vector2d a1 = 0.5 * J * (X3 - X2);
    const Vector2d a2v = 0.5 * J * (X1 - X3);
    const Vector2d a3 = 0.5 * J * (X2 - X1);

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

    // Gradient of energy
    const double trY = Y.trace();
    const double coeff = a2 - 0.5 * (a2 - b2) * trY;

    const Vector2d g1 = b2 * (L11 * x1 + L12 * x2 + L13 * x3) - coeff * (R * a1);
    const Vector2d g2 = b2 * (L21 * x1 + L22 * x2 + L23 * x3) - coeff * (R * a2v);
    const Vector2d g3 = b2 * (L31 * x1 + L32 * x2 + L33 * x3) - coeff * (R * a3);

    // Internal force = - gradient
    f.resize(6);
    f.segment<2>(0) = -g1;
    f.segment<2>(2) = -g2;
    f.segment<2>(4) = -g3;

    // Hessian of energy
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
        HE11 = a2 * (L11 - secondPartMat(s1, s1, A0, trY_safe));
        HE12 = a2 * (L12 - secondPartMat(s1, s2, A0, trY_safe));
        HE13 = a2 * (L13 - secondPartMat(s1, s3, A0, trY_safe));

        HE21 = a2 * (L21 - secondPartMat(s2, s1, A0, trY_safe));
        HE22 = a2 * (L22 - secondPartMat(s2, s2, A0, trY_safe));
        HE23 = a2 * (L23 - secondPartMat(s2, s3, A0, trY_safe));

        HE31 = a2 * (L31 - secondPartMat(s3, s1, A0, trY_safe));
        HE32 = a2 * (L32 - secondPartMat(s3, s2, A0, trY_safe));
        HE33 = a2 * (L33 - secondPartMat(s3, s3, A0, trY_safe));
    }
    else
    {
        const double coeffSecond = a2 - 0.5 * (a2 - b2) * trY;
        const double coeffAB     = 0.5 * (a2 - b2);

        HE11 = b2 * L11 - coeffSecond * secondPartMat(s1, s1, A0, trY_safe) + coeffAB * abSplitMat(r1, r1, A0);
        HE12 = b2 * L12 - coeffSecond * secondPartMat(s1, s2, A0, trY_safe) + coeffAB * abSplitMat(r1, r2, A0);
        HE13 = b2 * L13 - coeffSecond * secondPartMat(s1, s3, A0, trY_safe) + coeffAB * abSplitMat(r1, r3, A0);

        HE21 = b2 * L21 - coeffSecond * secondPartMat(s2, s1, A0, trY_safe) + coeffAB * abSplitMat(r2, r1, A0);
        HE22 = b2 * L22 - coeffSecond * secondPartMat(s2, s2, A0, trY_safe) + coeffAB * abSplitMat(r2, r2, A0);
        HE23 = b2 * L23 - coeffSecond * secondPartMat(s2, s3, A0, trY_safe) + coeffAB * abSplitMat(r2, r3, A0);

        HE31 = b2 * L31 - coeffSecond * secondPartMat(s3, s1, A0, trY_safe) + coeffAB * abSplitMat(r3, r1, A0);
        HE32 = b2 * L32 - coeffSecond * secondPartMat(s3, s2, A0, trY_safe) + coeffAB * abSplitMat(r3, r2, A0);
        HE33 = b2 * L33 - coeffSecond * secondPartMat(s3, s3, A0, trY_safe) + coeffAB * abSplitMat(r3, r3, A0);
    }

    // Tangent of force: H = df/dq = -d2E/dq2
    H.resize(6, 6);
    H.setZero();

    setBlock(H, 0, 0, -HE11);
    setBlock(H, 0, 1, -HE12);
    setBlock(H, 0, 2, -HE13);

    setBlock(H, 1, 0, -HE21);
    setBlock(H, 1, 1, -HE22);
    setBlock(H, 1, 2, -HE23);

    setBlock(H, 2, 0, -HE31);
    setBlock(H, 2, 1, -HE32);
    setBlock(H, 2, 2, -HE33);
}