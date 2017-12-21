/**
 * @file main.cpp
 * @brief solve plain elastic equations
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-12-17
 */

#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <string>
#include "fem.h"

//body force
double f1(cor& p, double h){
  return 0;
}
double f2(cor& p, double h){
  return 0;
}

//Diriclet BD
double u0(cor& p, double h){
  return 0;
}
double v0(cor& p, double h){
  return 0;
}

double F;

//Neumann BD
double gl1(cor& p, double h){
  if(fabs(p[0]) < 1e-8 && fabs(0.1-p[1]) < h)
    return F;
  else return 0;
  //return 0;
  //return F;
}
double gl2(cor& p, double h){
  return 0;
}
double gr1(cor& p, double h){
  return 0;
}
double gr2(cor& p, double h){
  return 0;
}
double gu1(cor& p, double h){
  return 0;
}
double gu2(cor& p, double h){
  if(fabs(0.1-p[1]) < 1e-8 && fabs(p[0]) < h)
    return -F;
  else return 0;
  //return 0;
  //return -F;
}

int main(int argc, char* argv[])
{
	if(argc != 5)
	{
		std::cerr << "Usage: " << argv[0]
			<< " <N> <uv_filename> <eps_sig_filename> <force> " << std::endl;
		return 1;
	}
	int N = atoi(argv[1]);
	std::string uv_filename(argv[2]);
	std::string eps_sig_filename(argv[3]);
	F = atof(argv[4]);

  std::cout << "Construct the problem..." << std::endl;
	fem X1(N);
  std::cout << "Create the linear equation..." << std::endl;
  X1.createLinEqu_triangle(f1, f2, u0, v0, gl1, gl2, gr1, gr2, gu1, gu2);
  std::cout << "Solve the linear equation..." << std::endl;
  X1.run();
  //std::cout << "Calculate the strain and the stress..." << std::endl;
  X1.cal_eps_sig();
  std::cout << "Print the solution to " << uv_filename << " and " << eps_sig_filename << "..." << std::endl;
  X1.print(uv_filename, eps_sig_filename, u0, v0);

  X1.print_geometry();

	return 0;
}

