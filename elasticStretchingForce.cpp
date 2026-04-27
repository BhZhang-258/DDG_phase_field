#include "elasticStretchingForce.h"
#include <iostream>

elasticStretchingForce::elasticStretchingForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

	I3 << 1,0,0,
          0,1,0,
          0,0,1;

	Z3 << 0,0,0,
          0,0,0,
          0,0,0;
}

elasticStretchingForce::~elasticStretchingForce()
{
	;
}

void elasticStretchingForce::computeFs()
{
	for (int i = 0; i < plate->v_triangularElement.size(); i++)
	{
		
		
		A.setZero();
		A = plate->getFFF(i, &fFFderiv, &fFFhess);
		for (int l = 0; l < 2; l++)
		{
			computeFs_layer(i, l);
		}
		
	}

	
}
void elasticStretchingForce::computeFs_layer(int idx, int l)
{
		VectorXi arrayIndex = plate->v_triangularElement[idx].arrayIndex;


		Layer layer = plate->v_triangularElement[idx].layers[l];

    	double thickness = layer.thickness;
    	double coeff = thickness / 4;


		
		Matrix2d Fg = layer.Fg;

		Matrix2d G = Fg.transpose() * Fg;
		Matrix2d Ginv = G.inverse();
		
		
		Matrix2d M = layer.M;
		Matrix2d Minv = layer.Minv;
		
		Matrix2d abar_g = M.transpose() * G * M;
		Matrix2d abarinv = abar_g.inverse();

		Matrix4d C = layer.C;

		Matrix4d kronc = Geometry::kron( Minv.transpose(), M * abarinv);
		// Matrix2d e;
		// e << 1, 2,
		// 	3, 4;

		// std::cout << M * abarinv * e * Minv <<  std::endl;
		// std::cout << kronc * Eigen::Map<Vector4d>(e.data())  <<  std::endl;
		// std::cout << kronc.transpose() * Eigen::Map<Vector4d>(e.data())  <<  std::endl;

		// std::cout << "L:\n" << M * abarinv << std::endl;
		// std::cout << "R:\n" << Minv << std::endl;
		// std::cout << "kronc:\n" << kronc << std::endl;



    	double dA = 0.5 * sqrt( abar_g.determinant() );
	
		Matrix<double, 9, 1> derivative;
		Matrix<double, 9, 9> hessian;
		

		Matrix2d Epsilion =  A - abar_g;
		Vector4d eps = Eigen::Map<Vector4d>(Epsilion.data());
		
		
		// stvk energy/grad/hessian  for epsilon inner product <eps,eps> with C_1 as the stiffness matrix
		double StVK = Constitutive::stvk<9>(C, kronc, 
			eps, eps,
			fFFderiv,fFFderiv, 
			fFFhess, fFFhess, 
			&derivative, &hessian);			

		// std::cout << "derivative:\n" << derivative.transpose() << std::endl;
		// std::cout << "hessian:\n" << hessian << std::endl;	
		derivative = coeff * dA * derivative;
		hessian = coeff * dA * hessian;



		for (int j = 0; j < 9; j++)
        {
            stepper->addForce(arrayIndex(j), derivative(j) );
        }
		

		for (int j = 0; j < 9; j++)
        {
            for (int k = 0; k < 9; k++)
            {
                stepper->addJacobian(arrayIndex(j), arrayIndex(k), hessian(j, k)  );
            }
        }
}

void elasticStretchingForce::computeJs()
{
}

void elasticStretchingForce::setFirstJacobian()
{
	for (int i = 0; i < plate->v_triangularElement.size(); i++)
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

