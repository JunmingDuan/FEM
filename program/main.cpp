/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  主函数
 *
 *        Version:  1.0
 *        Created:  2015年11月30日 23时01分58秒
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <string>
#include "fem.h"

double u0(cor& p){
	//return pow(p[0],3) + pow(p[1],3);
	return exp(p[0])*sin(2*M_PI*p[1]);
}

double f(cor& p){
	//return -6.*(p[0]+p[1]);
	return -(1.-4*M_PI*M_PI)*u0(p);
}

double beta(cor& p){
	return 1.;
}

double g(cor& p){
	//if(std::abs(p[0]-1)<1e-8) return 3*p[0]*p[0] + beta(p)*u0(p);
	//else if(std::abs(p[1]-1)<1e-8) return 3*p[1]*p[1] + beta(p)*u0(p);
	//else return -3*p[1]*p[1] + beta(p)*u0(p);
	if(std::abs(p[0]-1)<1e-8) return u0(p) + beta(p)*u0(p);
	else if(std::abs(p[1]-1)<1e-8) return  2.*M_PI*exp(p[0])*cos(2*M_PI*p[1]) + beta(p)*u0(p);
	else return -2.*M_PI*exp(p[0])*cos(2*M_PI*p[1]) + beta(p)*u0(p);
}

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0]
			<< " <N> <filename>" << std::endl;
		return 1;
	}
	int N = atoi(argv[1]);
	std::string filename(argv[2]);
	SOL NumSol((N+1)*(N+1));
	SOL Exact((N+1)*(N+1));
	for(int j = 0; j < N+1; ++j) {
		for(int i = 0; i < N+1; ++i) {
			cor x;
			x[0] = 1./N*i; x[1] = 1./N*j;
			Exact[j*(N+1)+i] = u0(x);
		}
	}

	clock_t t1,t2;
	t1=clock();
	fem X1(N);
	//X1.createLinEqu_triangle(f, u0, beta, g);
	X1.createLinEqu_rectangle(f, u0, beta, g);
	X1.solver(NumSol);
	t2=clock();

	std::ofstream outfile(filename.c_str(),std::ios::out);
	if(!outfile){
		std::cerr<<"Wrong!"<<std::endl;
	}
	outfile.precision(4);
	//outfile<<std::fixed;
	outfile<<std::showpos;
	outfile.setf(std::ios::scientific);
	outfile<< NumSol <<std::endl;
	outfile<< Exact <<std::endl;
	vvector<double> u1((N+1)*(N+1));
	vvector<double> u2((N+1)*(N+1));
	for(size_t i = 0; i < u1.size(); ++i) {
		u1[i] = NumSol[i];
		u2[i] = Exact[i];
	}
	std::cout.setf(std::ios::scientific);
	std::cout<<"norm1:"<<(u1-u2).Norm1()/u1.size()<<"\n";
	std::cout<<"norm2:"<<sqrt(pow((u1-u2).Norm2(),2)/u1.size())<<"\n";
	std::cout<<"norminf:"<<(u1-u2).NormInfinite()<<"\n";
	//std::cout<<(t2-t1)*pow(10,-6)<<std::endl;
	outfile.close();

	return 0;
}

