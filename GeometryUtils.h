#ifndef GEOMETRYUTILS_H
#define GEOMETRYUTILS_H

#include "eigenIncludes.h"
#include <optional>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

namespace Geometry 
{
    constexpr double PI = 3.14159265358979323846;
    // Input information for calculating the second fundamental form
    class sFFinformation
    {
    public:
        
        Eigen::Vector3d vertices[3];
        // 
        std::optional<Eigen::Vector3d> adjVertices[3];

        // 
        std::optional<Eigen::Vector3d> virtual_normal[3];

        // input the three vertices of the triangle
        sFFinformation(
            const Eigen::Vector3d& v0,
            const Eigen::Vector3d& v1,
            const Eigen::Vector3d& v2)
        {
            vertices[0] = v0;
            vertices[1] = v1;
            vertices[2] = v2;
        }
        // input the three adjacent vertices of the triangle
        void setupAdjVertex(int idx , Vector3d adjVertex) {
            adjVertices[idx] = adjVertex;
        }
        // input the Virtual vertices of the triangle
        void setupVirtualNormal(int idx , Vector3d normal) {
            virtual_normal[idx] = normal;
        }

        // get central vertex, idx = 0,1,2
        const Eigen::Vector3d& getVertex(int idx) const {
            return vertices[idx];
        }

        // get adjacent vertex, idx = 0,1,2
        const Eigen::Vector3d& getAdjVertex(int idx) const {
            return adjVertices[idx].value();
        }

        const Eigen::Vector3d& getVirtualNormal(int idx) const {
            return virtual_normal[idx].value();
        }


        // check if the adjacent vertex exists, idx = 0,1,2
        bool isInner(int idx) const {
            return adjVertices[idx].has_value();
        }

        bool isVirtualNormal(int idx) const {
            return virtual_normal[idx].has_value();
        }
    };
    
    void place3(Eigen::MatrixXd& M, int bi, int bj, Eigen::Matrix3d B);

    Vector3d faceNormal(
            const Vector3d& qi0,
            const Vector3d& qi1,
            const Vector3d& qi2,
            Matrix<double, 3, 9>* derivative = nullptr,
            std::vector<Matrix<double, 9, 9> >* hessian = nullptr);

    Matrix4d kron(const Matrix2d& A, const Matrix2d& B);
    
    Matrix3d crossMat(Vector3d a);
    Eigen::Matrix3d crossMatrix(Eigen::Vector3d v);

    Eigen::Matrix4d vecLeftMultiplyOperator(const Eigen::Matrix2d& A);

    Eigen::Matrix2d  firstFundamentalForm(const Vector2d &vi, const Vector2d &vj, const Vector2d &vk,
                    Matrix<double, 4, 6>* derivative,
                    std::vector<Matrix<double, 6, 6> >* hessian);

    Eigen::Matrix2d  firstFundamentalForm(const Vector3d &vi, const Vector3d &vj, const Vector3d &vk,
                    Matrix<double, 4, 9>* derivative,
                    std::vector<Matrix<double, 9, 9> >* hessian);

    Matrix2d secondFundamentalForm(
            const sFFinformation &sFFinf,
            Matrix<double, 4, 18>* derivative,
            std::vector<Matrix<double, 18, 18> >* hessian);

    Vector3d secondFundamentalFormEntries(
            const sFFinformation &sFFinf,
            Matrix<double, 3, 18>* derivative,
            std::vector<Matrix<double, 18, 18> >* hessian);


    //
    template <int N>
    void reorderMatrix(std::array<int, N> idx,
                    Eigen::Matrix<double, 1, N> *derivative, 
                    Eigen::Matrix<double, N, N> *hessian)
    {
        if (derivative)
        {
            Eigen::Matrix<double, 1, N> oldDeriv = *derivative;
            for (int i = 0; i < N; i++)
            {
                (*derivative)(0, idx[i]) = oldDeriv(0, i);
            }
        }

        if (hessian)
        {
            Eigen::Matrix<double, N, N> oldHess = *hessian;
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    (*hessian)(idx[i], idx[j]) = oldHess(i, j);
                }
            }
        }

    }


class UnionFind {
public:
    std::vector<int> parent;

    UnionFind(int n) { 
        parent.resize(n);
        for(int i = 0; i < n; i++) parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(int x, int y) {
        int rx = find(x);
        int ry = find(y);
        if (rx != ry) {
            
            if (rx < ry) parent[ry] = rx;
            else parent[rx] = ry;
        }
    }
};
    

}

#endif