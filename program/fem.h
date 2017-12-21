#ifndef _FEM_H
#define _FEM_H

#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "MAT.h"
#include "mvector.h"

const double lam = 3.65e10;
const double mu = 1.525e10;
//const double lam = 3.65;
//const double mu = 1.525;
typedef mvector<double,2> cor;
typedef Eigen::VectorXd SOL;
typedef Eigen::VectorXd EVEC;
typedef Eigen::MatrixXd EMAT;
typedef MAT<double> EN;
typedef MAT<double> CD;

class fem {
	private:
		int Nx, N;
    double h;
		typedef Eigen::Triplet<double> T;
		typedef double (*func)(cor&, double);
    EN en;
    CD cd;
		Eigen::SparseMatrix<double> A;
		Eigen::VectorXd b;
    SOL NumSol, Exact;
    EMAT eps, sig;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
    //Eigen::GMRES<Eigen::SparseMatrix<double>> solver;
    Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;

	public:
		fem(int Nx) : Nx(Nx) {
      N = 2*pow(Nx+1,2);
      h = 0.1/Nx;
			A.resize(N, N);
			b.resize(N);
      b.setZero();
			NumSol.resize(N);
      Exact.resize(N);
      eps.resize(2*pow(Nx,2),3);
      sig.resize(2*pow(Nx,2),3);
      triangle_en(Nx, en);
      triangle_cd(Nx, cd);
		}

    void triangle_en(int N, EN& en);
    void triangle_cd(int N, CD& cd);
		void createLinEqu_triangle(func, func , func , func,
        func, func, func , func , func, func);

    EVEC external_force(int , func, func);
    void NeumannBD(func, func, func, func, func, func);
    void DiricletBD(func, func);
    void DiricletBD_full(func u0, func v0);
    void fix_entry(const int m, const double tmp);
		void run();

    void cal_eps_sig();
    void print_geometry();
    //void print(const std::string&, const std::string&, func, func,
    void print(const std::string&, const std::string&, func, func);
};

#endif //_FEM_H

