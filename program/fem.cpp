/*
 * =====================================================================================
 *
 *       Filename:  fem.cpp
 *
 *    Description:  实现构造线性方程组并求解
 *
 *        Version:  1.0
 *        Created:  2015年11月30日 22时59分17秒
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include "fem.h"
#include "space.h"

void fem::createLinEqu_triangle(func f, func u0, func beta, func g) {
	EN en;
	CD cd;
	triangle_en(N, en);
	triangle_cd(N, cd);

	MAT<double> tmp(3, 2);
	MAT<double> Ke(3, 3);
	vvector<double> fe(3);
	std::vector< T > List;
	cor x0;
	x0[0] = 0; x0[1] = 0;
	for(int e = 0; e < N*N*2; ++e) {
		tmp[0][0] = cd[1][en[1][e]] - cd[1][en[2][e]];
		tmp[0][1] = cd[0][en[2][e]] - cd[0][en[1][e]];
		tmp[1][0] = cd[1][en[2][e]] - cd[1][en[0][e]];
		tmp[1][1] = cd[0][en[0][e]] - cd[0][en[2][e]];
		tmp[2][0] = cd[1][en[0][e]] - cd[1][en[1][e]];
		tmp[2][1] = cd[0][en[1][e]] - cd[0][en[0][e]];
		//此处下两行不加时，f为常数
		x0[0] = 1./3*(cd[0][en[0][e]] + cd[0][en[1][e]] + cd[0][en[2][e]]);
		x0[1] = 1./3*(cd[1][en[0][e]] + cd[1][en[1][e]] + cd[1][en[2][e]]);
		fe[0] = f(x0)/6./N/N;
		fe[1] = f(x0)/6./N/N;
		fe[2] = f(x0)/6./N/N;
		for(int i = 0; i < 3; ++i) {
			for(int j = 0; j < 3; ++j) {
				Ke[i][j] = 0.5*N*N*(tmp[i][0]*tmp[j][0] + tmp[i][1]*tmp[j][1]);
				List.push_back( T(en[i][e], en[j][e], Ke[i][j]) );
			}
			b(en[i][e]) += fe[i];
		}
	}

	A.setFromTriplets(List.begin(), List.end());
	A.makeCompressed();

	//右边界
	for(int j = 0; j < N+1; ++j) {
		int m = N + j*(N+1);
		A.coeffRef(m, m) += 2./3/N;
		if(j < N) {
			A.coeffRef(m, m+N+1) += 1./6/N;
			A.coeffRef(m+N+1, m) += 1./6/N;
		}
		cor x;
		if(j == 0) {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			//中点积分
			b(m) += g(x)*0.5/N;
		}
		else if(j == N) {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			b(m) += g(x)*0.5/N;
		}
		else {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			b(m) += 2*g(x)*0.5/N;
		}
	}
	//下边界
	for(int j = 1; j < N+1; ++j) {
		if(j < N) A.coeffRef(j, j) += 2./3/N;
		A.coeffRef(j, j-1) += 1./6/N;
		A.coeffRef(j-1, j) += 1./6/N;
		cor x;
		if(j == N) {
			x[0] = cd[0][j] - 0.5/N;
			x[1] = cd[1][j];
			b(j) += g(x)*0.5/N;
		}
		else {
			x[0] = cd[0][j];
			x[1] = cd[1][j];
			b(j) += 2*g(x)*0.5/N;
		}

	}
	//上边界
	for(int j = 1; j < N+1; ++j) {
		int m = j + N*(N+1);
		if(j < N) A.coeffRef(m, m) += 2./3/N;
		A.coeffRef(m, m-1) += 1./6/N;
		A.coeffRef(m-1, m) += 1./6/N;
		cor x;
		if(j == N) {
			x[0] = cd[0][m] - 0.5/N;
			x[1] = cd[1][m];
			b(m) += g(x)*0.5/N;
		}
		else {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			b(m) += 2*g(x)*0.5/N;
		}
	}

	for(int j = 0; j < N+1; ++j) {
		int m = j*(N+1);
		cor x0;
		x0[0] = cd[0][m]; x0[1] = cd[1][m];
		for(int i = 0; i < A.rows(); ++i) {
			if(i != m && A.coeff(i,m) != 0) {
				b(i) -= A.coeff(i,m)*u0(x0);
				A.coeffRef(i, m) = 0;
			}
			if(i != m && A.coeff(m,i) != 0) {
				A.coeffRef(m, i) = 0;
			}
		}
		A.coeffRef(m, m) = 1;
		b(m) = u0(x0);
	}
	A.prune(0.);

}

void fem::createLinEqu_rectangle(func f, func u0, func beta, func g) {
	EN en;
	CD cd;
	rectangle_en(N, en);
	rectangle_cd(N, cd);

	MAT<double> Ke(4, 4);
	vvector<double> fe(4);
	std::vector< T > List;
	cor x0;
	x0[0] = 0; x0[1] = 0;
	Ke[0][0] = 2./3; Ke[1][1] = 2./3; Ke[2][2] = 2./3; Ke[3][3] = 2./3;
	Ke[1][0] = -1./6; Ke[0][1] = -1./6; Ke[2][1] = -1./6; Ke[1][2] = -1./6; Ke[3][2] = -1./6; Ke[2][3] = -1./6;
	Ke[2][0] = -1./3; Ke[0][2] = -1./3; Ke[3][1] = -1./3; Ke[1][3] = -1./3;
	Ke[3][0] = -1./6; Ke[0][3] = -1./6;
	for(int e = 0; e < N*N; ++e) {
		//此处f为常数
		x0[0] = 1./4*(cd[0][en[0][e]] + cd[0][en[1][e]] + cd[0][en[2][e]]+ cd[0][en[3][e]]);
		x0[1] = 1./4*(cd[1][en[0][e]] + cd[1][en[1][e]] + cd[1][en[2][e]]+ cd[1][en[3][e]]);
		fe[0] = f(x0)/4./N/N; fe[1] = f(x0)/4./N/N; fe[2] = f(x0)/4./N/N; fe[3] = f(x0)/4./N/N;
		for(int i = 0; i < 4; ++i) {
			for(int j = 0; j < 4; ++j) {
				List.push_back( T(en[i][e], en[j][e], Ke[i][j]) );
			}
			b(en[i][e]) += fe[i];
		}
	}

	A.setFromTriplets(List.begin(), List.end());
	A.makeCompressed();

	//右边界
	for(int j = 0; j < N+1; ++j) {
		int m = N + j*(N+1);
		A.coeffRef(m, m) += 2./3/N;
		if(j < N) {
			A.coeffRef(m, m+N+1) += 1./6/N;
			A.coeffRef(m+N+1, m) += 1./6/N;
		}
		cor x;
		if(j == 0) {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			//中点积分
			b(m) += g(x)*0.5/N;
		}
		else if(j == N) {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			b(m) += g(x)*0.5/N;
		}
		else {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			b(m) += 2*g(x)*0.5/N;
		}
	}
	//下边界
	for(int j = 1; j < N+1; ++j) {
		if(j < N) A.coeffRef(j, j) += 2./3/N;
		A.coeffRef(j, j-1) += 1./6/N;
		A.coeffRef(j-1, j) += 1./6/N;
		cor x;
		if(j == N) {
			x[0] = cd[0][j] - 0.5/N;
			x[1] = cd[1][j];
			b(j) += g(x)*0.5/N;
		}
		else {
			x[0] = cd[0][j];
			x[1] = cd[1][j];
			b(j) += 2*g(x)*0.5/N;
		}

	}
	//上边界
	for(int j = 1; j < N+1; ++j) {
		int m = j + N*(N+1);
		if(j < N) A.coeffRef(m, m) += 2./3/N;
		A.coeffRef(m, m-1) += 1./6/N;
		A.coeffRef(m-1, m) += 1./6/N;
		cor x;
		if(j == N) {
			x[0] = cd[0][m] - 0.5/N;
			x[1] = cd[1][m];
			b(m) += g(x)*0.5/N;
		}
		else {
			x[0] = cd[0][m];
			x[1] = cd[1][m];
			b(m) += 2*g(x)*0.5/N;
		}
	}

	for(int j = 0; j < N+1; ++j) {
		int m = j*(N+1);
		cor x0;
		x0[0] = cd[0][m]; x0[1] = cd[1][m];
		for(int i = 0; i < A.rows(); ++i) {
			if(i != m && A.coeff(i,m) != 0) {
				b(i) -= A.coeff(i,m)*u0(x0);
				A.coeffRef(i, m) = 0;
			}
			if(i != m && A.coeff(m,i) != 0) {
				A.coeffRef(m, i) = 0;
			}
		}
		A.coeffRef(m, m) = 1;
		b(m) = u0(x0);
	}
	A.prune(0.);

}

void fem::solver(SOL& sol) {
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	sol = solver.solve(b);
	//std::cout<<sol<<std::endl;
}
