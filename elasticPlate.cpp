#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_Possion, double m_dt)
{
	YoungM = m_YoungM;
	density = m_density;
	Possion = m_Possion;
	dt = m_dt;

	setupGeometry();

	ndof_u = 2 * nv;
	ndof_phi = nv;

	x = VectorXd::Zero(ndof_u);
	x0 = VectorXd::Zero(ndof_u);
	u = VectorXd::Zero(ndof_u);

	phi = VectorXd::Zero(ndof_phi);
	phi_old = VectorXd::Zero(ndof_phi);
	phi_prev = VectorXd::Zero(ndof_phi);

	for (int i = 0; i < nv; i++)
	{
		x(2 * i + 0) = v_nodes[i](0);
		x(2 * i + 1) = v_nodes[i](1);
	}
	x0 = x;
	xStart = x;

	buildElemet();

	setupMass();

	//set up constraint map
	isConstrained_u = VectorXi::Zero(ndof_u);
	isConstrained_phi = VectorXi::Zero(ndof_phi);
}

elasticPlate::~elasticPlate()
{
	;
}

void elasticPlate::setup()
{
	ncons_u = 0;
	ncons_phi = 0;
    for (int i=0; i < ndof_u; i++)
    {
		if (isConstrained_u[i] > 0)
		{
			ncons_u++;
		}
	}
	for (int i=0; i < ndof_phi; i++)
    {
		if (isConstrained_phi[i] > 0)
		{
			ncons_phi++;
		}
	}
	uncons_u = ndof_u - ncons_u;
	uncons_phi = ndof_phi - ncons_phi;
	
	fullToUnconsMap_u = VectorXi::Constant(ndof_u, -1);
	unconstrainedMap_u = VectorXi::Constant(uncons_u, -1);
	
	fullToUnconsMap_phi = VectorXi::Constant(ndof_phi, -1);
	unconstrainedMap_phi = VectorXi::Constant(uncons_phi, -1);

	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof_u; i++)
	{
		if (isConstrained_u[i] == 0)
		{
			unconstrainedMap_u[c] = i;
			fullToUnconsMap_u[i] = c;
			c++;
		}
	}
	c = 0;
	for (int i=0; i < ndof_phi; i++)
	{
		if (isConstrained_phi[i] == 0)
		{
			unconstrainedMap_phi[c] = i;
			fullToUnconsMap_phi[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{

	massArray = VectorXd::Zero(ndof_u);

	for (int i = 0; i < ne; i++)
	{
		int nv1 = v_triangularElement[i].nv_1;
		int nv2 = v_triangularElement[i].nv_2;
		int nv3 = v_triangularElement[i].nv_3;

		Vector2d x1 = getVertex(nv1);
		Vector2d x2 = getVertex(nv2);
		Vector2d x3 = getVertex(nv3);

		double cross = (x2 - x1).x() * (x3 - x1).y() - (x2 - x1).y() * (x3 - x1).x();
		double area = 0.5 * std::abs(cross);
		double deltaMass = area * density / 3;

		massArray(2 * nv1 + 0) = massArray(2 * nv1 + 0) + deltaMass;
		massArray(2 * nv1 + 1) = massArray(2 * nv1 + 1) + deltaMass;

		massArray(2 * nv2 + 0) = massArray(2 * nv2 + 0) + deltaMass;
		massArray(2 * nv2 + 1) = massArray(2 * nv2 + 1) + deltaMass;

		massArray(2 * nv3 + 0) = massArray(2 * nv3 + 0) + deltaMass;
		massArray(2 * nv3 + 1) = massArray(2 * nv3 + 1) + deltaMass;
	}

}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained_u[k];
}
int elasticPlate::getIfConstrained_phi(int k)
{
	return isConstrained_phi[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector2d position, int k)
{
	isConstrained_u[2 * k + 0] = 1;
	isConstrained_u[2 * k + 1] = 1;

	// Store in the constrained dof vector
	x(2 * k + 0) = position(0);
	x(2 * k + 1) = position(1);
}

void elasticPlate::setConstraint(double position, int k)
{
	isConstrained_u[k] = 1;

	// Store in the constrained dof vector
	x(k) = position;
}

void elasticPlate::setOneVertexBoundaryCondition(double position, int i, int k)
{
	isConstrained_u[2 * i + k] = 1;

	// Store in the constrained dof vector
	x(2 * i + k) = position;
}

void elasticPlate::setPhiBoundaryCondition(double position, int i)
{
	isConstrained_u[i] = 1;

	// Store in the constrained dof vector
	x(i) = position;
}

Vector2d elasticPlate::getVertex(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x(2 * i + 0);
	xCurrent(1) = x(2 * i + 1);

	return xCurrent;
}

double elasticPlate::getPhi(int i)
{
	return phi(i);
}

Vector2d elasticPlate::getVertexStart(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = xStart(2 * i + 0);
	xCurrent(1) = xStart(2 * i + 1);

	return xCurrent;

}

Vector2d elasticPlate::getVertexOld(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x0(2 * i + 0);
	xCurrent(1) = x0(2 * i + 1);

	return xCurrent;
}

Vector2d elasticPlate::getVelocity(int i)
{
	Vector2d uCurrent;

	uCurrent(0) = ( x(2 * i + 0) - x0(2 * i + 0) ) / dt;
	uCurrent(1) = ( x(2 * i + 1) - x0(2 * i + 1) ) / dt;

	return uCurrent;
}

void elasticPlate::prepareForIteration()
{
	;
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;


	// make sure th phi can only increase
	// for (int i = 0; i < nv; i++)
	// {
	// 	double localPhi_now = x(2 * nv + i);
	// 	double localPhi_before = x0(2 * nv + i);

	// 	if (localPhi_before > localPhi_now)
	// 	{
	// 		x(2 * nv + i) = localPhi_before;
	// 		x0(2 * nv + i) = localPhi_now;
	// 	}

	// }

	// update x
	x0 = x;

	// make sure phi cannot be larger than 1.0
	// for (int i = 0; i < nv; i++)
	// {
	// 	double localPhi = getPhi(i);

	// 	if (localPhi > 1.0)
	// 	{
	// 		x(2 * nv + i) = 1.0;
	// 		x0(2 * nv + i) = 1.0;
	// 	}

	// }
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons_u; c++)
	{
		x[unconstrainedMap_u[c]] = x[unconstrainedMap_u[c]] + u[unconstrainedMap_u[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons_u; c++)
	{
		x[unconstrainedMap_u[c]] -= m_motion[c];
	}
}
void elasticPlate::updateNewtonMethod_phi(VectorXd m_motion)
{
	for (int c=0; c < uncons_phi; c++)
	{
		phi[unconstrainedMap_phi[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	ifstream inFile1;
	inFile1.open("inputdata/nodesInput.txt");
	v_nodes.clear();
	double a, b;
	while(inFile1 >> a >> b)
	{
		Vector2d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;
		v_nodes.push_back(xCurrent);
	}
	nv = v_nodes.size();
	inFile1.close();

	ifstream inFile2;
	inFile2.open("inputdata/elementInput.txt");
	v_element.clear();
	int aa, bb, cc;
	while(inFile2 >> aa >> bb >> cc)
	{
		aa = aa;
		bb = bb;
		cc = cc;

		Vector3i elementCurrent;
		elementCurrent(0) = aa - 1;
		elementCurrent(1) = bb - 1;
		elementCurrent(2) = cc - 1;

		v_element.push_back(elementCurrent);
	}
	ne = v_element.size();
	inFile2.close();
}

void elasticPlate::buildElemet()
{

	v_triangularElement.clear();

	for (int i = 0; i < ne; i++)
	{
		Vector3i elementCurrent = v_element[i];

		basicElement m_basicElement;

		m_basicElement.nv_1 = elementCurrent(0);
		m_basicElement.nv_2 = elementCurrent(1);
		m_basicElement.nv_3 = elementCurrent(2);

		Vector2d x1 = getVertex(m_basicElement.nv_1);
		Vector2d x2 = getVertex(m_basicElement.nv_2);
		Vector2d x3 = getVertex(m_basicElement.nv_3);

		Vector2d e1 = x2 - x1;
		Vector2d e2 = x3 - x1;
		double A0 = 0.5 * std::abs(e1.x() * e2.y() - e1.y() * e2.x());
		m_basicElement.area = A0;
		m_basicElement.abar = Geometry::firstFundamentalForm( x1 , x2 , x3 , nullptr,  nullptr );
		m_basicElement.abarinv = m_basicElement.abar.inverse();
		
		double fib_theta = 0.0;  // Local material coordinate system
		Vector2d u; // fiber direction 1
		u << cos(fib_theta), sin(fib_theta);

		Vector2d v; // fiber direction 2
		v << -sin(fib_theta), cos(fib_theta);

		// The rotation matrix of edge coordinate system(e1,e2) to material (u,v)
		m_basicElement.M << e1.dot(u), e2.dot(u),
     		e1.dot(v), e2.dot(v);

		
		double d = 0;  // fiber stiffness ratio


		double E_1 = YoungM * (1+d);
		double E_2 = YoungM;
		double nu12 = Possion;
		double G12 = E_2 / ( 2 * (1 + Possion) );
		double nu21 = nu12 * (E_2 / E_1);
		double denom = 1 - nu12 * nu21;
		double C11 = E_1 / denom;
		double C22 = E_2 / denom;
		double C12 = (nu12 * E_2) / denom;  // or (nu21 * E_1) / denom
		double mu_shear = G12;

		m_basicElement.C << 
		C11, 0,           0,           C12,
		0,   0,           2*mu_shear,  0,
		0,   2*mu_shear,  0,           0,
		C12, 0,           0,           C22;

		double lambda = YoungM * Possion / ( 1 - Possion * Possion );
		double mu = YoungM / ( 2 * (1 + Possion) );

		m_basicElement.C1 << 
		lambda, 0,      0,      lambda,
		0,   	0,      0,  	0,
		0,   	0,  	0,      0,
		lambda, 0,      0,      lambda;
		m_basicElement.C2 << 
		2*mu,	0, 		0,      0,
		0,   	0,   	2*mu,  	0,
		0,   	2*mu,  	0,   	0,
		0, 		0,      0,      2*mu;
		//
		// Mphi
		// --------------------------
		Matrix3d Mphi;
		Mphi << 2.0, 1.0, 1.0,
				1.0, 2.0, 1.0,
				1.0, 1.0, 2.0;
		Mphi *= (A0 / 12.0);
		m_basicElement.Mphi = Mphi;

		// --------------------------
		// Bphi
		// --------------------------
		Vector2d gradN1, gradN2, gradN3;
		gradN1 << (x2.y() - x3.y()), (x3.x() - x2.x());
		gradN2 << (x3.y() - x1.y()), (x1.x() - x3.x());
		gradN3 << (x1.y() - x2.y()), (x2.x() - x1.x());

		gradN1 /= (2.0 * A0);
		gradN2 /= (2.0 * A0);
		gradN3 /= (2.0 * A0);

		Eigen::Matrix<double, 2, 3> Bphi;
		Bphi.col(0) = gradN1;
		Bphi.col(1) = gradN2;
		Bphi.col(2) = gradN3;

		m_basicElement.Bphi = Bphi;

		// --------------------------
		// Kphi
		// --------------------------
		Matrix3d Kphi = A0 * (Bphi.transpose() * Bphi);
		m_basicElement.Kphi = Kphi;

	//
		
		m_basicElement.arrayIndex = VectorXi::Zero(9);

		m_basicElement.arrayIndex(0) = 2 * m_basicElement.nv_1 + 0;
		m_basicElement.arrayIndex(1) = 2 * m_basicElement.nv_1 + 1;

		m_basicElement.arrayIndex(2) = 2 * m_basicElement.nv_2 + 0;
		m_basicElement.arrayIndex(3) = 2 * m_basicElement.nv_2 + 1;

		m_basicElement.arrayIndex(4) = 2 * m_basicElement.nv_3 + 0;
		m_basicElement.arrayIndex(5) = 2 * m_basicElement.nv_3 + 1;

		// m_basicElement.arrayIndex(6) = 2 * nv + m_basicElement.nv_1;
		// m_basicElement.arrayIndex(7) = 2 * nv + m_basicElement.nv_2;
		// m_basicElement.arrayIndex(8) = 2 * nv + m_basicElement.nv_3;

		v_triangularElement.push_back(m_basicElement);	
	}


}

