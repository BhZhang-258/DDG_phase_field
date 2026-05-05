#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "elasticPlate.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <unordered_set>


#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <GL/glut.h>

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

struct SolverConfig {
    int total_dof;                     // 全局 DOF 数
    int unconstrained_dof;             // 自由 DOF 数

    VectorXi fullToUncons;     // 映射（-1 表示约束）
    VectorXi isConstrained;   // 或 bool

    int mtype = 11;                   // PARDISO 类型（实对称/非对称）
};

class timeStepper
{
public:

    timeStepper(const SolverConfig& cfg);
    ~timeStepper();

    VectorXd GlobalForceVec;
    MatrixXd GlobalJacobianMax;
    VectorXd GlobalMotionVec;
    VectorXd ReForceVec;

    SparseMatrix<double> sm1;

    void setZero();
    void addForce(int ind, double p);
    void addFirstJacobian(int ind1, int ind2);
    void finalizePattern();
    void addJacobian(int ind1, int ind2, double p);
    void integrator();

    void first_time_PARDISO_setup();

private:
    SolverConfig config;
    int total_dof;
    int uncon_dof;

    int mappedInd, mappedInd1, mappedInd2;
    std::vector<std::unordered_set<int>> rowPattern;

    int n;
    int *ia;
    int *ja;
    double *a;

    int *ia1;
    int *ja1;

    MKL_INT mtype;

    double *b;
    double *x;

    MKL_INT nrhs;

    void *pt[64]; 

    MKL_INT iparm[64];
    double dparm[64];

    MKL_INT maxfct, mnum, phase, error, msglvl;

    double ddum;

    MKL_INT idum;

    int total_num;

    void pardisoSolver();

    VectorXi num1;
    VectorXi num2;

    VectorXi num11;
    VectorXi num22;
};

#endif
