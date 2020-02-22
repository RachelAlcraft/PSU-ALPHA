//
// This file contains some example code, which may be used to call 
// the matrix, vector and statistics algorithms.
//
// This code is not necessarily meant to be useful in itself, but is provided as 
// an example of how the class may be used.
//
// Copyright (C) 1999 Software Engineering Group, Crystallography Department,
// Birkbeck College, Malet Street, London WC1E 7HX, U.K.
// (d.moss@mail.cryst.bbk.ac.uk or m.williams@biochemistry.ucl.ac.uk)
// 
// This library is free software; you can redistribute it and/or modify it 
// under the terms of the GNU Library General Public License as published by 
// the Free Software Foundation; either version 2 of the License, or (at your
// Handle) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
/////////////////////////////////////////////////////////////////////////////////////
//
// Author: Mark Williams 
// 
/////////////////////////////////////////////////////////////////////////////////////
// 
// Brief Description of Code:
// 
// Carry out simple numerical tasks.
//
// For further details see the documentation for the matrix, vector classes 
// and the algorithms header files.
/////////////////////////////////////////////////////////////////////////////////////
// 
#include <iostream>
using namespace std;

#include "btl_numeric_vector.h"
#include "btl_matrix.h"
#include "btl_vector_algorithms.h"
#include "btl_matrix_algorithms.h"
#include "btl_statistics.h"
using namespace btl;

int main()
{

	matrix<double> A(3,4,2.3);
	matrix<double> B(4,3,2.6);
	matrix<double> C(3,3,1.2);
	matrix<double> D(3,4);

	// member function - A,B,C are instances of the matrix class
	C = A * B;

// Vector stuff

	double temp = 0.0;
	temp = scalar_product(A.begin(),A.end(),B.begin(),temp);
	cout << "temp " << temp << endl;

	temp = 0.0;
	temp = scalar_product(A.begin(),B.begin(),temp);
	cout << "temp triple " << temp << endl;

	vector_product(A.begin(),B.begin(),C.begin());
	cout << "vector product " << C(1,1) << endl;

	triple_vector_product(A.begin(),B.begin(),C.begin(),D.begin());
	cout << "triple vector product " << D(1,1) << " " << D(1,2) << 
	" " << D(1,3) << endl;

	direct_product(A.begin(),A.end(),B.begin(),D.begin());
	cout << "direct " << D(1,1) << endl;

	D += temp;
	cout << "direct+temp " << D(1,1) << endl;

	temp = 0.0;
	temp = triple_scalar_product(A.begin(),B.begin(),D.begin(),temp);
	cout << "triple_scalar " << temp << endl;

	temp = 0.0;
	temp = separation_squared(A.begin(),A.end(),B.begin(),temp);
	cout << "separation squared " << temp << endl;

	temp = 0.0;
	temp = separation(A.begin(),A.end(),B.begin(),temp);
	cout << "separation " << temp << endl;

	temp = 0.0;
	temp = sum(A.begin(),A.end(),temp);
	cout << "sum "<< temp << endl;
	
	temp = 0.0;
	temp = sum_precise(A.begin(),A.end(),temp);
	cout << "sum precise "<< temp << endl;

	temp = 0.0;
	temp = sum_of_squares(A.begin(),A.end(),temp);
	cout << "sum squares "<< temp << endl;

	temp = 0.0;
	temp = sum_of_squares_precise(A.begin(),A.end(),temp);
	cout << "sum squares precise "<< temp << endl;

	temp = 0.0;
	temp = sum(A.begin(),temp);
	cout << "sum triple "<< temp << endl;

	temp = 0.0;
	temp = sum_of_squares(A.begin(),temp);
	cout << "sum squares triple "<< temp << endl;

	temp = 0.0;
	temp = magnitude(A.begin(),A.end(),temp);
	cout << "magnitude "<< temp << endl;

	temp = 0.0;
	temp = magnitude_precise(A.begin(),A.end(),temp);
	cout << "magnitude precise "<< temp << endl;

	temp = 0.0;
	temp = magnitude(A.begin(),temp);
	cout << "magnitude triple "<< temp << endl;

// Statistics stuff

	double mean1 = 0.0;
	mean1 = mean(A.begin(),A.end(),mean1);
	cout << "mean "<< mean1 << endl;

	temp = 0.0;
	temp = mean_absolute_deviation(A.begin(),A.end(),temp,mean1);
	cout << "mean abs dev "<< temp << endl;

	temp = 0.0;
	temp = variance(A.begin(),A.end(),temp,mean1);
	cout << "variance "<< temp << endl;

	double temp1=0.0;
	double temp2=0.0;
	double temp3=0.0;
	double temp4=0.0;
	double temp5=0.0;
	normal_statistics(A.begin(),A.end(),temp1,temp2,temp3,temp4,temp5);
	cout << "mean "<< temp1 << endl;
	cout << "mad "<< temp2 << endl;
	cout << "variance "<< temp3 << endl;
	cout << "skew "<< temp4 << endl;
	cout << "kurtosis "<< temp << endl;

// Matrix stuff

	vector<double> E(5);
	cout << "A " << A(1,1) << endl;
	column_means(A.begin(),A.end(),A.rows(),E.begin());
	cout << "A column mean " << E[0] << endl;

	copy_column(A.begin(),A.end(),A.rows(),E.begin(),1);		
	cout << "A copied to E " << A(1,1) << " " << E[1] << endl;
	
	temp=0.0;
	temp=determinant(C.begin(),C.end(),C.rows(),temp);		
	cout << "det C "<< temp << endl;
	cout << "matrix C "<< C << endl;

	matrix<double> F(4,4);
	matrix<double> G(4,4);
	matrix<double> H(4,4);
	matrix<double> I(4,4);
	matrix<double> J(4,4);
	matrix<double> K(3,3);
	int i,j;
	temp=0.0;
	for(i=1;i<5;i++)
	{
	    for(j=1;j<5;j++)
	    {
		temp += 1.0;
		F(i,j)=temp;
	    }
	}
	cout << "matrix F "<< F << endl;

	transpose(F.begin(),F.end(),F.rows(),G.begin());
	cout << "matrix G=FT "<< G << endl;

	matrix_matrixtranspose_product(F.begin(),F.end(),F.rows(),
					G.begin(),G.end(),G.rows(),
					H.begin());
	cout << "matrix FGT "<< H << endl;

	matrixtranspose_matrix_product(F.begin(),F.end(),F.rows(),
					G.begin(),G.end(),G.rows(),
					I.begin());
	cout << "matrix FTG "<< I << endl;

	matrix_product(G.begin(),G.end(),G.rows(),
			G.begin(),G.end(),G.rows(),
			J.begin());
	cout << "matrix GG "<< J << endl;

	matrix_minor(J.begin(),J.end(),J.rows(),K.begin(),1,3);
	cout << "matrix minor of GG "<< K << endl;

	matrix<double> L(4,4,1.3);
	temp=0.0;
	K(1,1)= 1.0;
	K(1,2)= 3.0;
	K(1,3)= 7.0;
	K(2,1)= 3.0;
	K(2,2)= 4.0;
	K(2,3)= 5.0;
	K(3,1)= 3.0;
	K(3,2)= 2.0;
	K(3,3)= 2.0;
	cout << "K "<< K << endl;
	temp=determinant(K.begin(),K.end(),K.rows(),temp);		
	cout << "det K "<< temp << endl;

	matrix<double> M(3,3,1.3);
	adjoint(K.begin(),K.end(),K.rows(),M.begin());		
	cout << "adjoint K "<< M << endl;
	inverse_square(K.begin(),K.end(),K.rows(),M.begin());		
	cout << "inverse K "<< M << endl;

	matrix<double> N(3,3,0.0);
	matrix_product(M.begin(),M.end(),M.rows(),
			K.begin(),K.end(),K.rows(),
			N.begin());
	cout << "invK.K "<< N << endl;

	matrix<double> U;
	numeric_vector<double> X;
	matrix<double> V;

	K(1,1)= 4.0; K(1,2)= 2.0; K(1,3)= 14.0;
	K(2,1)= 2.0; K(2,2)= 17.0; K(2,3)= -5.0;
	K(3,1)= 14.0; K(3,2)= -5.0; K(3,3)= 83.0;

	inverse_cholesky(K.begin(),K.end(),K.rows(),M.begin());		
	cout << "inverse K "<< M << endl;

	N *= 0.0;
	matrix_product(M.begin(),M.end(),M.rows(),
			K.begin(),K.end(),K.rows(),
			N.begin());
	cout << "invK.K "<< N << endl;

	K(1,1)= 0.0;
	K(1,2)= 1.0;
	K(1,3)= 0.0;
	K(2,1)= 1.0;
	K(2,2)= 0.0;
	K(2,3)= 0.0;
	K(3,1)= 0.0;
	K(3,2)= 0.0;
	K(3,3)= 0.0;
	numeric_vector<double> P(3);
	eigen_solution(K.begin(),K.end(),K.rows(),M.begin(),P.begin());
	cout << "eigenvectors " << M << endl;
	cout << "eigenvalues " << P << endl;


	

	assert(C(1,2)==C(2,1));
	cout << "multiply " << C(1,2) << endl;

    	return 0;
}
