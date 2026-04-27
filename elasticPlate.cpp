#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_Possion, double m_dt)
{
	YoungM = m_YoungM;
	density = m_density;
	Possion = m_Possion;
	dt = m_dt;

	setupGeometry();

	ndof = 3 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

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
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{

	massArray = VectorXd::Zero(ndof);

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

	for (int i = 0; i < nv; i++)
	{
		massArray(2 * nv + i) = 1e-5;
	}


}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector2d position, int k)
{
	isConstrained[2 * k + 0] = 1;
	isConstrained[2 * k + 1] = 1;
	//isConstrained[3 * k + 2] = 1;

	// Store in the constrained dof vector
	x(2 * k + 0) = position(0);
	x(2 * k + 1) = position(1);
	//x(3 * k + 2) = position(2);
}

void elasticPlate::setConstraint(double position, int k)
{
	isConstrained[k] = 1;

	// Store in the constrained dof vector
	x(k) = position;
}

void elasticPlate::setOneVertexBoundaryCondition(double position, int i, int k)
{
	isConstrained[2 * i + k] = 1;

	// Store in the constrained dof vector
	x(2 * i + k) = position;
}

void elasticPlate::setPhiBoundaryCondition(double position, int i)
{
	isConstrained[2 * nv + i] = 1;

	// Store in the constrained dof vector
	x(2 * nv + i) = position;
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
	return x(2 * nv + i);
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
	for (int i = 0; i < nv; i++)
	{
		double localPhi_now = x(2 * nv + i);
		double localPhi_before = x0(2 * nv + i);

		if (localPhi_before > localPhi_now)
		{
			x(2 * nv + i) = localPhi_before;
			x0(2 * nv + i) = localPhi_now;
		}

	}

	// update x
	x0 = x;

	// make sure phi cannot be larger than 1.0
	for (int i = 0; i < nv; i++)
	{
		double localPhi = getPhi(i);

		if (localPhi > 1.0)
		{
			x(2 * nv + i) = 1.0;
			x0(2 * nv + i) = 1.0;
		}

	}
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
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

		Matrix2d Dm;
    	Dm.col(0) = x2 - x1;
    	Dm.col(1) = x3 - x1;

    	const double detDm = Dm.determinant();
    	if (detDm <= 0.0)
    	{
        	m_basicElement.nv_1 = elementCurrent(0);
			m_basicElement.nv_2 = elementCurrent(2);
			m_basicElement.nv_3 = elementCurrent(1);
    	}

		m_basicElement.X = MatrixXd::Zero(2,3);

		m_basicElement.X.col(0) = getVertex(m_basicElement.nv_1);
		m_basicElement.X.col(1) = getVertex(m_basicElement.nv_2);
		m_basicElement.X.col(2) = getVertex(m_basicElement.nv_3);
		
		m_basicElement.arrayIndex = VectorXi::Zero(9);

		m_basicElement.arrayIndex(0) = 2 * m_basicElement.nv_1 + 0;
		m_basicElement.arrayIndex(1) = 2 * m_basicElement.nv_1 + 1;

		m_basicElement.arrayIndex(2) = 2 * m_basicElement.nv_2 + 0;
		m_basicElement.arrayIndex(3) = 2 * m_basicElement.nv_2 + 1;

		m_basicElement.arrayIndex(4) = 2 * m_basicElement.nv_3 + 0;
		m_basicElement.arrayIndex(5) = 2 * m_basicElement.nv_3 + 1;

		m_basicElement.arrayIndex(6) = 2 * nv + m_basicElement.nv_1;
		m_basicElement.arrayIndex(7) = 2 * nv + m_basicElement.nv_2;
		m_basicElement.arrayIndex(8) = 2 * nv + m_basicElement.nv_3;

		v_triangularElement.push_back(m_basicElement);	
	}

}