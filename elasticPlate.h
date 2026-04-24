#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct basicElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	MatrixXd X;

	VectorXi arrayIndex;
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

	void setOneVertexBoundaryCondition(double position, int i, int k);

	VectorXd x;
	VectorXd x0;
	VectorXd u;
	VectorXd xStart;

	std::vector<Vector2d> v_nodes;
	std::vector<Vector3i> v_element;

	std::vector<basicElement> v_triangularElement;

	int nv;
	int ne;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector2d position, int k);
	void setConstraint(double position, int k);

	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	VectorXd massArray;
	void setupMass();

	void buildElemet();

	private:
};

#endif
