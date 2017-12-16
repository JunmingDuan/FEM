/*
 * =====================================================================================
 *
 *       Filename:  fem.h
 *
 *    Description:  生成有限元线性方程组
 *
 *        Version:  1.0
 *        Created:  2015年11月30日 17时26分14秒
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _FEM_H
#define _FEM_H

#include "mvector.h"
#include "vvector.h"
#include <Eigen/SparseLU>
#include <Eigen/SparseCore>

typedef mvector<double,2> cor;
typedef Eigen::VectorXd SOL;

class fem {
	private:
		int N;
		typedef Eigen::Triplet<double> T;
		typedef double (*func)(cor&);
		Eigen::SparseMatrix<double> A;
		Eigen::VectorXd b;


	public:
		fem(int N) : N(N) { 
			A.resize((N+1)*(N+1), (N+1)*(N+1));
			b.resize((N+1)*(N+1));
		}

		void createLinEqu_triangle(func f, func u0, func beta, func g);

		void createLinEqu_rectangle(func f, func u0, func beta, func g);

		void solver(SOL& sol);
};

#endif //_FEM_H
