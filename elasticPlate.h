#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include "GeometryUtils.h"
#include <fstream>

struct basicElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	Matrix4d C;
	Matrix4d C1; // Volume term
	Matrix4d C2; // Shear energy term
	Matrix2d M;

	Matrix2d abar;
	Matrix2d abarinv;

	double area; // Initial Geometry area

	Matrix3d Mphi;
	Matrix<double, 2, 3> Bphi;
	Matrix3d Kphi;
	double H_history;

	VectorXi arrayIndex;
	Vector3i arrayIndex_phi;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_Possion, double m_dt);
	~elasticPlate();

	double YoungM;
	double Possion;
	double dt;
	double density;

	double alpha;
	double beta;

	Vector2d getVertex(int i);
	Vector2d getVertexOld(int i);
	Vector2d getVertexStart(int i);
	Vector2d getVelocity(int i);

	double getPhi(int i);

	void setOneVertexBoundaryCondition(double position, int i, int k);
	void setPhiBoundaryCondition(double position, int i);

	std::vector<Vector2d> v_nodes;
	std::vector<Vector3i> v_element;

	std::vector<basicElement> v_triangularElement;

	int nv;
	int ne;

	int ndof_u;
	int uncons_u;
	int ncons_u;

	int ndof_phi;
	int uncons_phi;
	int ncons_phi;

	VectorXd x;
	VectorXd x0;
	VectorXd u;
	VectorXd xStart;

	VectorXd phi;       
	VectorXd phi_old;   
	VectorXd phi_prev;  
	// VectorXd dphi;      

	void setupGeometry();

	void setVertexBoundaryCondition(Vector2d position, int k);
	void setConstraint(double position, int k);

	// boundary conditions
	VectorXi isConstrained_u;
	VectorXi isConstrained_phi;
	int getIfConstrained(int k);
    int getIfConstrained_phi(int k);
    VectorXi unconstrainedMap_u;
    VectorXi unconstrainedMap_phi;
	VectorXi fullToUnconsMap_u;
	VectorXi fullToUnconsMap_phi;

	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
    void updateNewtonMethod_phi(VectorXd m_motion);
    void prepareForIteration();

    VectorXd massArray;
	void setupMass();

	void buildElemet();

    Matrix2d getFFF(int idx, Matrix<double, 4, 9> *derivative, std::vector<Matrix<double, 9, 9>> *hessian);

private:
};

#endif
