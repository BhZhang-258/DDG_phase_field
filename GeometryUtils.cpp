#include "GeometryUtils.h"

namespace Geometry {

void place3(MatrixXd& M, int bi, int bj, Matrix3d B) 
{
  M.block<3,3>(3*bi, 3*bj) = B;
}

Vector3d faceNormal(
        const Vector3d& q0,
        const Vector3d& q1,
        const Vector3d& q2,
        Matrix<double, 3, 9>* derivative,
        std::vector<Matrix<double, 9, 9> >* hessian)
{
    if (derivative)
    {
        derivative->setZero();
    }

    if (hessian)
    {
        hessian->resize(3);
        for (int i = 0; i < 3; i++) 
        {
            (*hessian)[i].setZero();
        }
    }

    Vector3d n = (q1 - q0).cross(q2 - q0);

    if (derivative)
    {
        derivative->block(0, 0, 3, 3) += crossMat(q2 - q1);
        derivative->block(0, 3, 3, 3) += crossMat(q0 - q2);
        derivative->block(0, 6, 3, 3) += crossMat(q1 - q0);
    }

    if (hessian)
    {
        for (int j = 0; j < 3; j++)
        {
            Vector3d ej(0, 0, 0);
            ej[j] = 1.0;
            Matrix3d ejc = crossMat(ej);
            (*hessian)[j].block(0, 3, 3, 3) -= ejc;
            (*hessian)[j].block(0, 6, 3, 3) += ejc;
            (*hessian)[j].block(3, 6, 3, 3) -= ejc;
            (*hessian)[j].block(3, 0, 3, 3) += ejc;
            (*hessian)[j].block(6, 0, 3, 3) -= ejc;
            (*hessian)[j].block(6, 3, 3, 3) += ejc;
        }
    }

    return n;
}

Matrix4d kron(const Matrix2d &A, const Matrix2d &B)
{
    Matrix4d K;

    K.block<2,2>(0,0) = A(0,0) * B;
    K.block<2,2>(0,2) = A(0,1) * B;
    K.block<2,2>(2,0) = A(1,0) * B;
    K.block<2,2>(2,2) = A(1,1) * B;

    return K;
}

Matrix3d crossMat(Vector3d a)
{
    Matrix3d b;

    b<<0,-a(2),a(1),
    a(2),0,-a(0),
    -a(1),a(0),0;

    return b;
}

Eigen::Matrix3d crossMatrix(Eigen::Vector3d v)
{
    Eigen::Matrix3d ret;
    ret << 0, -v[2], v[1],
        v[2], 0, -v[0],
        -v[1], v[0], 0;
    return ret;
}

Eigen::Matrix4d vecLeftMultiplyOperator(const Eigen::Matrix2d &A)
{
    // This function computes the 4x4 matrix T such that T * vec(M) = vec(A*M) for any 2x2 matrix M . 
    // where vec(M) is the vectorization of M by stacking its columns.
    Eigen::Matrix4d T;

    T << A(0,0), 0,      A(0,1),  0,
         0,      A(0,0), 0,       A(0,1),
         A(1,0), 0,      A(1,1),  0,
         0,      A(1,0), 0,       A(1,1);

    return T;
}

Eigen::Matrix2d  firstFundamentalForm(const Vector3d &vi, const Vector3d &vj, const Vector3d &vk,
                                      Matrix<double, 4, 9>* derivative,
                                      std::vector<Matrix<double, 9, 9> >* hessian)
{
    Matrix3d I3;
    Matrix3d Z3;
    
    I3 << 1,0,0,
          0,1,0,
          0,0,1;

	Z3 << 0,0,0,
          0,0,0,
          0,0,0;

    Vector3d e_1 = vj - vi;  
	Vector3d e_2 = vk - vi; 

	Matrix2d A;

    A(0,0) = e_1.dot(e_1);
    A(0,1) = e_1.dot(e_2);
    A(1,0) = e_2.dot(e_1);
    A(1,1) = e_2.dot(e_2);
    if (derivative)
    {
        derivative->setZero();

        VectorXd gradA1, gradA2, gradA3, gradA4;
        gradA1 = VectorXd::Zero(9);
        gradA2 = VectorXd::Zero(9);
        gradA3 = VectorXd::Zero(9);
        gradA4 = VectorXd::Zero(9);

        gradA1.segment(0,3) = - 2 * e_1; 
        gradA1.segment(3,3) =   2 * e_1; 

        gradA2.segment(0,3) = - e_1 - e_2;
        gradA2.segment(3,3) =         e_2;        
        gradA2.segment(6,3) =   e_1;


        gradA3.segment(0,3) = - e_2 - e_1;
        gradA3.segment(3,3) =   e_2;
        gradA3.segment(6,3) =         e_1;

        gradA4.segment(0,3) = - 2 * e_2;
        gradA4.segment(6,3) =   2 * e_2;

        
        derivative->block<1, 9>(0, 0) = gradA1.transpose();
        derivative->block<1, 9>(1, 0) = gradA2.transpose();
        derivative->block<1, 9>(2, 0) = gradA3.transpose();
        derivative->block<1, 9>(3, 0) = gradA4.transpose();
        
    }
    if (hessian)
    {
        hessian->resize(4);
        for (int i = 0; i < 4; i++)
        {
            (*hessian)[i].setZero();
        }


        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        (*hessian)[0].block<3, 3>(0, 0) += 2.0 * I;
        (*hessian)[0].block<3, 3>(3, 3) += 2.0 * I;
        (*hessian)[0].block<3, 3>(0, 3) -= 2.0 * I;
        (*hessian)[0].block<3, 3>(3, 0) -= 2.0 * I;

        (*hessian)[1].block<3, 3>(3, 6) += I;
        (*hessian)[1].block<3, 3>(6, 3) += I;
        (*hessian)[1].block<3, 3>(0, 3) -= I;
        (*hessian)[1].block<3, 3>(0, 6) -= I;
        (*hessian)[1].block<3, 3>(3, 0) -= I;
        (*hessian)[1].block<3, 3>(6, 0) -= I;
        (*hessian)[1].block<3, 3>(0, 0) += 2.0 * I;

        (*hessian)[2].block<3, 3>(3, 6) += I;
        (*hessian)[2].block<3, 3>(6, 3) += I;
        (*hessian)[2].block<3, 3>(0, 3) -= I;
        (*hessian)[2].block<3, 3>(0, 6) -= I;
        (*hessian)[2].block<3, 3>(3, 0) -= I;
        (*hessian)[2].block<3, 3>(6, 0) -= I;
        (*hessian)[2].block<3, 3>(0, 0) += 2.0 * I;

        (*hessian)[3].block<3, 3>(0, 0) += 2.0 * I;
        (*hessian)[3].block<3, 3>(6, 6) += 2.0 * I;
        (*hessian)[3].block<3, 3>(0, 6) -= 2.0 * I;
        (*hessian)[3].block<3, 3>(6, 0) -= 2.0 * I;
    }


    return A;

}
Eigen::Matrix2d  firstFundamentalForm(const Vector2d &vi, const Vector2d &vj, const Vector2d &vk,
                                      Matrix<double, 4, 6>* derivative,
                                      std::vector<Matrix<double, 6, 6> >* hessian)
{
    using namespace Eigen;

    Vector2d e_1 = vj - vi;
    Vector2d e_2 = vk - vi;

    Matrix2d A;

    A(0,0) = e_1.dot(e_1);
    A(0,1) = e_1.dot(e_2);
    A(1,0) = e_2.dot(e_1);
    A(1,1) = e_2.dot(e_2);

    // ================= derivative =================
    if (derivative)
    {
        derivative->setZero();

        VectorXd gradA1 = VectorXd::Zero(6);
        VectorXd gradA2 = VectorXd::Zero(6);
        VectorXd gradA3 = VectorXd::Zero(6);
        VectorXd gradA4 = VectorXd::Zero(6);

        // A11
        gradA1.segment<2>(0) = -2 * e_1;
        gradA1.segment<2>(2) =  2 * e_1;

        // A12
        gradA2.segment<2>(0) = -e_1 - e_2;
        gradA2.segment<2>(2) =  e_2;
        gradA2.segment<2>(4) =  e_1;

        // A21（和 A12 一样）
        gradA3.segment<2>(0) = -e_2 - e_1;
        gradA3.segment<2>(2) =  e_2;
        gradA3.segment<2>(4) =  e_1;

        // A22
        gradA4.segment<2>(0) = -2 * e_2;
        gradA4.segment<2>(4) =  2 * e_2;

        derivative->block<1,6>(0,0) = gradA1.transpose();
        derivative->block<1,6>(1,0) = gradA2.transpose();
        derivative->block<1,6>(2,0) = gradA3.transpose();
        derivative->block<1,6>(3,0) = gradA4.transpose();
    }

    // ================= hessian =================
    if (hessian)
    {
        hessian->resize(4);
        for (int i = 0; i < 4; i++)
        {
            (*hessian)[i].setZero();
        }

        Matrix2d I = Matrix2d::Identity();

        // A11
        (*hessian)[0].block<2,2>(0,0) += 2.0 * I;
        (*hessian)[0].block<2,2>(2,2) += 2.0 * I;
        (*hessian)[0].block<2,2>(0,2) -= 2.0 * I;
        (*hessian)[0].block<2,2>(2,0) -= 2.0 * I;

        // A12
        (*hessian)[1].block<2,2>(2,4) += I;
        (*hessian)[1].block<2,2>(4,2) += I;
        (*hessian)[1].block<2,2>(0,2) -= I;
        (*hessian)[1].block<2,2>(0,4) -= I;
        (*hessian)[1].block<2,2>(2,0) -= I;
        (*hessian)[1].block<2,2>(4,0) -= I;
        (*hessian)[1].block<2,2>(0,0) += 2.0 * I;

        // A21（同 A12）
        (*hessian)[2].block<2,2>(2,4) += I;
        (*hessian)[2].block<2,2>(4,2) += I;
        (*hessian)[2].block<2,2>(0,2) -= I;
        (*hessian)[2].block<2,2>(0,4) -= I;
        (*hessian)[2].block<2,2>(2,0) -= I;
        (*hessian)[2].block<2,2>(4,0) -= I;
        (*hessian)[2].block<2,2>(0,0) += 2.0 * I;

        // A22
        (*hessian)[3].block<2,2>(0,0) += 2.0 * I;
        (*hessian)[3].block<2,2>(4,4) += 2.0 * I;
        (*hessian)[3].block<2,2>(0,4) -= 2.0 * I;
        (*hessian)[3].block<2,2>(4,0) -= 2.0 * I;
    }

    return A;
}

Matrix2d secondFundamentalForm(
        const sFFinformation &sFFinf,
        Matrix<double, 4, 18>* derivative,
        std::vector<Matrix<double, 18, 18> >* hessian)
    {
        if (derivative)
        {
            derivative->resize(4, 18);
            derivative->setZero();
        }

        if (hessian)
        {
            hessian->resize(4);
            for (int i = 0; i < 4; i++)
            {
                (*hessian)[i].resize(18, 18);
                (*hessian)[i].setZero();
            }
        }

        Matrix<double, 3, 18> IIderiv;
        std::vector < Matrix<double, 18, 18> > IIhess;

        Vector3d II = secondFundamentalFormEntries(sFFinf,  &IIderiv , &IIhess );

        Matrix2d result;
        result << II[0] + II[1], II[0], II[0], II[0] + II[2];

        if (derivative)
        {
            derivative->row(0) += IIderiv.row(0);
            derivative->row(0) += IIderiv.row(1);

            derivative->row(1) += IIderiv.row(0);
            derivative->row(2) += IIderiv.row(0);

            derivative->row(3) += IIderiv.row(0);
            derivative->row(3) += IIderiv.row(2);
        }
        if (hessian)
        {
            (*hessian)[0] += IIhess[0];
            (*hessian)[0] += IIhess[1];

            (*hessian)[1] += IIhess[0];
            (*hessian)[2] += IIhess[0];

            (*hessian)[3] += IIhess[0];
            (*hessian)[3] += IIhess[2];
        }

        return result;
    }

Vector3d secondFundamentalFormEntries(
        const sFFinformation &sFFinf,
        Matrix<double, 3, 18>* derivative,
        std::vector<Matrix<double, 18, 18> >* hessian)
    {
        if (derivative)
        {
            derivative->setZero();
        }

        if (hessian)
        {
            hessian->resize(3);
            for (int i = 0; i < 3; i++)
            {
                (*hessian)[i].setZero();
            }
        }

        Vector3d II;

        Vector3d oppNormals[3];
        oppNormals[0].setZero();
        oppNormals[1].setZero();
        oppNormals[2].setZero();

        Matrix<double, 3, 9> dn[3];
        std::vector<Matrix<double, 9, 9> > hn[3];

        dn[0].setZero();
        dn[1].setZero();
        dn[2].setZero();

        
        for (int i = 0; i < 3; ++i) 
        {
            hn[i].resize(3);

            for (int j = 0; j < 3; ++j) 
            {
                hn[i][j] = Matrix<double, 9, 9>::Zero(); 
            }
        }

        Matrix<double, 3, 9> dcn;
        std::vector<Matrix<double, 9, 9> > hcn;

        
        Vector3d qi = sFFinf.getVertex(0);
        Vector3d qj = sFFinf.getVertex(1);
        Vector3d qk = sFFinf.getVertex(2);
        
        Matrix3d qc;
        qc.col(0) = qi;
        qc.col(1) = qj;
        qc.col(2) = qk;

        Vector3d cNormal = Geometry::faceNormal(qi, qj, qk,   &dcn , &hcn );

        for (int ii = 0; ii < 3; ii++)
        {
            
            int idx = ii;


            if (sFFinf.isInner(idx)) 
            {
                Vector3d v1 = sFFinf.getAdjVertex(idx); 
                
                Vector3d v2 = sFFinf.getVertex((idx + 2) % 3); 

                Vector3d v3 = sFFinf.getVertex((idx + 1) % 3); 

                oppNormals[ii] = Geometry::faceNormal(v1,v2,v3,  &dn[ii] , &hn[ii]);
            }
            else if (sFFinf.isVirtualNormal(idx)) 
            {
                
                Vector3d oppNormal = sFFinf.getVirtualNormal(idx);
                oppNormals[ii] = oppNormal;
            }
        }

        Vector3d qs[3];
        Vector3d mvec[3];
        Vector3d qvec[3];
        double mnorms[3];

        qs[0] = qi;
        qs[1] = qj;
        qs[2] = qk;

        for (int i = 0; i < 3; i++)
        {
            mvec[i] = oppNormals[i] + cNormal;
            mnorms[i] = mvec[i].norm();
        }


        for (int i = 0; i < 3; i++)
        {
            int ip1 = (i + 1) % 3;
            int ip2 = (i + 2) % 3;
            qvec[i] = qs[ip1] + qs[ip2] - 2.0 * qs[i];
        }

        for (int i = 0; i < 3; i++)
        {
            int ip1 = (i + 1) % 3;
            int ip2 = (i + 2) % 3;
            II[i] = (qs[ip1] + qs[ip2] - 2.0 * qs[i]).dot(oppNormals[i]) / mnorms[i];
            if (derivative)
            {
                derivative->block<1, 3>(i, 3 * i) += -2.0 * oppNormals[i].transpose() / mnorms[i];
                derivative->block<1, 3>(i, 3 * ip1) += 1.0 * oppNormals[i].transpose() / mnorms[i];
                derivative->block<1, 3>(i, 3 * ip2) += 1.0 * oppNormals[i].transpose() / mnorms[i];

                derivative->block<1, 3>(i, 9 + 3 * i) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 0);
                derivative->block<1, 3>(i, 3 * ip2) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 3);
                derivative->block<1, 3>(i, 3 * ip1) += qvec[i].transpose() / mnorms[i] * dn[i].block<3, 3>(0, 6);

                derivative->block<1, 3>(i, 9 + 3 * i) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 0); // 9+3i是i的对顶点，也就是qi_bar
                derivative->block<1, 3>(i, 3 * ip2) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 3);   // ip2是i在三角内部左侧的点，表明dn中第二个索引应该组装到qk
                derivative->block<1, 3>(i, 3 * ip1) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dn[i].block<3, 3>(0, 6);

                derivative->block<1, 3>(i, 0) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 0);
                derivative->block<1, 3>(i, 3) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 3);
                derivative->block<1, 3>(i, 6) += -qvec[i].dot(oppNormals[i]) / mnorms[i] / mnorms[i] / mnorms[i] * mvec[i].transpose() * dcn.block<3, 3>(0, 6);
            }
            if (hessian)
            {
                int ip1 = (i + 1) % 3;
                int ip2 = (i + 2) % 3;

                int miidx[3];
                miidx[0] = 9 + 3 * i;
                miidx[1] = 3 * ip2;
                miidx[2] = 3 * ip1;

                Matrix3d dnij[3];
                for (int j = 0; j < 3; j++)
                    dnij[j] = dn[i].block<3, 3>(0, 3 * j);

                for (int j = 0; j < 3; j++)
                {

                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip1) += (1.0 / mnorms[i]) * dnij[j].transpose();
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip2) += (1.0 / mnorms[i]) * dnij[j].transpose();
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * i) += (-2.0 / mnorms[i]) * dnij[j].transpose();

                    Vector3d dnijTm = dnij[j].transpose() * mvec[i];
                    Vector3d dcnjTm = (dcn.block<3, 3>(0, 3 * j).transpose() * mvec[i]);
                    Matrix3d term3 = dnijTm * oppNormals[i].transpose();

                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip1) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * ip2) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;
                    (*hessian)[i].block<3, 3>(miidx[j], 3 * i) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term3;

                    Matrix3d term4 = dcnjTm * oppNormals[i].transpose();

                    (*hessian)[i].block<3, 3>(3 * j, 3 * ip1) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;
                    (*hessian)[i].block<3, 3>(3 * j, 3 * ip2) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;
                    (*hessian)[i].block<3, 3>(3 * j, 3 * i) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term4;

                    (*hessian)[i].block<3, 3>(3 * ip1, miidx[j]) += (1.0 / mnorms[i]) * dnij[j];
                    (*hessian)[i].block<3, 3>(3 * ip2, miidx[j]) += (1.0 / mnorms[i]) * dnij[j];
                    (*hessian)[i].block<3, 3>(3 * i, miidx[j]) += (-2.0 / mnorms[i]) * dnij[j];

                    for (int k = 0; k < 3; k++)
                    {
                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * (dnij[j].transpose() * mvec[i]) * (qvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(3 * j, miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * (dcn.block<3, 3>(0, 3 * j).transpose() * mvec[i]) * (qvec[i].transpose() * dnij[k]);
                    }

                    Matrix3d term1 = oppNormals[i] * dnijTm.transpose();
                    (*hessian)[i].block<3, 3>(3 * ip1, miidx[j]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;
                    (*hessian)[i].block<3, 3>(3 * ip2, miidx[j]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;
                    (*hessian)[i].block<3, 3>(3 * i, miidx[j]) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term1;

                    Matrix3d term2 = oppNormals[i] * dcnjTm.transpose();
                    (*hessian)[i].block<3, 3>(3 * ip1, 3 * j) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;
                    (*hessian)[i].block<3, 3>(3 * ip2, 3 * j) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;
                    (*hessian)[i].block<3, 3>(3 * i, 3 * j) += (2.0 / mnorms[i] / mnorms[i] / mnorms[i]) * term2;

                    Vector3d dnijTq = dnij[j].transpose() * qvec[i];

                    for (int k = 0; k < 3; k++)
                    {
                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * dnijTq * (mvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(miidx[j], 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * dnijTq * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));
                    }

                    double qdoto = qvec[i].dot(oppNormals[i]);

                    for (int k = 0; k < 3; k++)
                    {
                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnij[j].transpose() * dnij[k];
                        (*hessian)[i].block<3, 3>(miidx[j], 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnij[j].transpose() * dcn.block<3, 3>(0, 3 * k);
                        (*hessian)[i].block<3, 3>(3 * j, miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcn.block<3, 3>(0, 3 * j).transpose() * dnij[k];
                        (*hessian)[i].block<3, 3>(3 * j, 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcn.block<3, 3>(0, 3 * j).transpose() * dcn.block<3, 3>(0, 3 * k);

                        (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnijTm * (mvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(miidx[j], 3 * k) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dnijTm * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));
                        (*hessian)[i].block<3, 3>(3 * j, miidx[k]) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcnjTm * (mvec[i].transpose() * dnij[k]);
                        (*hessian)[i].block<3, 3>(3 * j, 3 * k) += (3.0 / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * dcnjTm * (mvec[i].transpose() * dcn.block<3, 3>(0, 3 * k));

                        for (int l = 0; l < 3; l++)
                        {
                            (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * mvec[i][l] * hn[i][l].block<3, 3>(3 * j, 3 * k);
                            (*hessian)[i].block<3, 3>(3 * j, 3 * k) += (-1.0 / mnorms[i] / mnorms[i] / mnorms[i]) * qdoto * mvec[i][l] * hcn[l].block<3, 3>(3 * j, 3 * k);
                            (*hessian)[i].block<3, 3>(miidx[j], miidx[k]) += (1.0 / mnorms[i]) * qvec[i][l] * hn[i][l].block<3, 3>(3 * j, 3 * k);
                        }
                    }
                }
            }
        }

        return II;
    }

}