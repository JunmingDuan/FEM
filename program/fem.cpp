/**
 * @file fem.cpp
 * @brief construct linear equation
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-12-17
 */

#include "fem.h"

void fem::triangle_en(int N, EN& en) {
	en.resize(6, N*N*2);
	for(int j = 0; j < N; ++j) {
		for(int i = 0; i < N; ++i) {
			en[0][j*N*2+i*2] = 2*(N+1)*j + 2*i;
			en[1][j*N*2+i*2] = 2*(N+1)*j + 2*i + 1;
			en[2][j*N*2+i*2] = 2*(N+1)*j + 2*i + 2*(N+2);
			en[3][j*N*2+i*2] = 2*(N+1)*j + 2*i + 2*(N+2) + 1;
			en[4][j*N*2+i*2] = 2*(N+1)*j + 2*i + 2*(N+1);
			en[5][j*N*2+i*2] = 2*(N+1)*j + 2*i + 2*(N+1) + 1;
			en[0][j*N*2+i*2+1] = 2*(N+1)*j + 2*i;
			en[1][j*N*2+i*2+1] = 2*(N+1)*j + 2*i + 1;
			en[2][j*N*2+i*2+1] = 2*(N+1)*j + 2*i + 2;
			en[3][j*N*2+i*2+1] = 2*(N+1)*j + 2*i + 3;
			en[4][j*N*2+i*2+1] = 2*(N+1)*j + 2*i + 2*(N+2);
			en[5][j*N*2+i*2+1] = 2*(N+1)*j + 2*i + 2*(N+2) + 1;
		}
	}
}

void fem::triangle_cd(int N, CD& cd) {
	cd.resize(2, 2*(N+1)*(N+1));
	for(int j = 0; j <= N; ++j) {
		for(int i = 0; i <= N; ++i) {
			cd[0][2*(N+1)*j+2*i] = i*h;
			cd[1][2*(N+1)*j+2*i] = j*h;
			cd[0][2*(N+1)*j+2*i+1] = i*h;
			cd[1][2*(N+1)*j+2*i+1] = j*h;
		}
	}
}


void fem::createLinEqu_triangle(func f1, func f2, func u0, func v0,
    func gl1, func gl2, func gr1, func gr2, func gu1, func gu2) {
  double b1, b2, b3, c1, c2, c3;
  EMAT B(3,6), L(3,3), Ke(6,6);
  B.setZero();
  L.setZero();
  Ke.setZero();
  L(0,0) = lam+2*mu; L(0,1) = lam;
  L(1,0) = lam; L(1,1) = lam+2*mu; L(2,2) = mu;

  vvector<double> fe(3);
  std::vector< T > List;
  cor x0;
  x0[0] = 0; x0[1] = 0;
  for(int e = 0; e < Nx*Nx*2; ++e) {
    b1 = cd[1][en[2][e]] - cd[1][en[4][e]];
    c1 = cd[0][en[2][e]] - cd[0][en[4][e]];
    b2 = cd[1][en[4][e]] - cd[1][en[0][e]];
    c2 = cd[0][en[4][e]] - cd[0][en[0][e]];
    b3 = cd[1][en[0][e]] - cd[1][en[2][e]];
    c3 = cd[0][en[0][e]] - cd[0][en[2][e]];
    B(0,0) = b1; B(0,2) = b2; B(0,4) = b3;
    B(1,1) = -c1; B(1,3) = -c2; B(1,5) = -c3;
    B(2,0) = -c1; B(2,1) = b1; B(2,2) = -c2;
    B(2,3) = b2; B(2,4) = -c3; B(2,5) = b3;
    double D = pow(h,2);
    Ke = 0.5/D*B.transpose()*L*B;

    EVEC fe = external_force(e, f1, f2);

    for(int i = 0; i < 6; ++i) {
      for(int j = 0; j < 6; ++j) {
        List.push_back( T(en[i][e], en[j][e], Ke(i,j)) );
      }
      b(en[i][e]) += fe[i];
    }
  }

  A.setFromTriplets(List.begin(), List.end());

  //先把边界积分加到右端项,再做固定边界条件的系数修正
  NeumannBD(gl1, gl2, gr1, gr2, gu1, gu2);
  DiricletBD(u0, v0);

  //四条边Diriclet边界(debug)
  //DiricletBD_full(u0, v0);

  A.prune(1e-13);
  A.makeCompressed();

  //std::cout << A << std::endl;
  //std::cout << b << std::endl;
}

EVEC fem::external_force(int e, func f1, func f2) {
  EVEC fe(6);
  fe.setZero();
  //2阶代数精度,3个点积分
  std::vector<cor> vertex(3);
  std::vector<cor> p(3);
  std::vector<double> weight(3,1./3);
  vertex[0][0] = cd[0][en[0][e]];
  vertex[0][1] = cd[1][en[0][e]];
  vertex[1][0] = cd[0][en[2][e]];
  vertex[1][1] = cd[1][en[2][e]];
  vertex[2][0] = cd[0][en[4][e]];
  vertex[2][1] = cd[1][en[4][e]];
  p[0][0] = (vertex[0][0]+vertex[1][0])/2.;
  p[0][1] = (vertex[0][1]+vertex[1][1])/2.;
  p[1][0] = (vertex[1][0]+vertex[2][0])/2.;
  p[1][1] = (vertex[1][1]+vertex[2][1])/2.;
  p[2][0] = (vertex[2][0]+vertex[0][0])/2.;
  p[2][1] = (vertex[2][1]+vertex[0][1])/2.;
  //3个点上函数f的两个分量f1,f2的值
  std::vector<double> f1p(3), f2p(3);
  for(int i = 0; i < 3; ++i) {
    f1p[i] = f1(p[i], h);
    f2p[i] = f2(p[i], h);
  }
  //3个点上3个基函数phi的值
  std::vector<double> phi1p(3), phi2p(3), phi3p(3);
  phi1p[0] = 0.5; phi1p[1] = 0; phi1p[2] = 0.5;
  phi2p[0] = 0.5; phi2p[1] = 0.5; phi2p[2] = 0;
  phi3p[0] = 0; phi3p[1] = 0.5; phi3p[2] = 0.5;
  //数值积分
  for(int j = 0; j < 3; ++j) {
    fe[0] += weight[j]*f1p[j]*phi1p[j];
    fe[1] += weight[j]*f2p[j]*phi1p[j];
    fe[2] += weight[j]*f1p[j]*phi2p[j];
    fe[3] += weight[j]*f2p[j]*phi2p[j];
    fe[4] += weight[j]*f1p[j]*phi3p[j];
    fe[5] += weight[j]*f2p[j]*phi3p[j];
  }
  fe *= pow(h,2)/2.;

  return fe;
}

/*EVEC fem::external_force(int e, func f1, func f2) {*/
  //EVEC fe(6);
  ////2阶代数精度,6个点积分
  //std::vector<cor> p[6];
  //p[0][0] = cd[0][en[0][e]];
  //p[0][1] = cd[1][en[0][e]];
  //p[1][0] = cd[0][en[2][e]];
  //p[1][1] = cd[1][en[2][e]];
  //p[2][0] = cd[0][en[4][e]];
  //p[2][1] = cd[1][en[4][e]];
  //p[3][0] = (p[0][0]+p[1][0])/2.;
  //p[3][1] = (p[0][1]+p[1][1])/2.;
  //p[4][0] = (p[1][0]+p[2][0])/2.;
  //p[4][1] = (p[1][1]+p[2][1])/2.;
  //p[5][0] = (p[2][0]+p[0][0])/2.;
  //p[5][1] = (p[2][1]+p[0][1])/2.;
  ////6个点上函数f的两个分量f1,f2的值
  //std::vector<double> f1p[6], f2p[6];
  //for(int i = 0; i < 6; ++i) {
    //f1p[i] = f1(p[i]);
    //f2p[i] = f2(p[i]);
  //}
  ////6个点上3个基函数phi的值
  //std::vector<double> phi1p[6], phi2p[6], phi3p[6];
  //phi1p[0] = 1; phi1p[1] = 0; phi1p[2] = 0;
  //phi1p[3] = 0.5; phi1p[4] = 0; phi1p[5] = 0.5;
  //phi2p[0] = 0; phi2p[1] = 1; phi2p[2] = 0;
  //phi2p[3] = 0.5; phi2p[4] = 0.5; phi2p[5] = 0;
  //phi3p[0] = 0; phi3p[1] = 0; phi3p[2] = 1;
  //phi3p[3] = 0; phi3p[4] = 0.5; phi3p[5] = 0.5;
  ////数值积分

  //return fe;
//}

void fem::NeumannBD(func gl1, func gl2, func gr1, func gr2, func gu1, func gu2) {
  //线积分,2个Gauss积分点
  cor x1, x2, g1, g2;
  std::vector<double> phi1p(2), phi2p(2), phi3p(2);
  //左边界
  for(int j = 0; j < Nx; ++j) {
    phi1p[0] = (1+1/sqrt(3))/2; phi1p[1] = (1-1/sqrt(3))/2;
    phi3p[0] = (1-1/sqrt(3))/2; phi3p[1] = (1+1/sqrt(3))/2;
    x1[0] = cd[0][en[0][2*Nx*j]];
    x1[1] = cd[1][en[0][2*Nx*j]];
    x2[0] = cd[0][en[4][2*Nx*j]];
    x2[1] = cd[1][en[4][2*Nx*j]];
    g1 = -(x2-x1)/2./sqrt(3) + (x2+x1)/2.;
    g2 =  (x2-x1)/2./sqrt(3) + (x2+x1)/2.;
    int m = j*2*(Nx+1);
    b(m) += (gl1(g1, h)*phi1p[0]+gl1(g2, h)*phi1p[1])*h/2;
    m = j*2*(Nx+1)+1;
    b(m) += (gl2(g1, h)*phi1p[0]+gl2(g2, h)*phi1p[1])*h/2;
    m = (j+1)*2*(Nx+1);
    b(m) += (gl1(g1, h)*phi3p[0]+gl1(g2, h)*phi3p[1])*h/2;
    m = (j+1)*2*(Nx+1)+1;
    b(m) += (gl2(g1, h)*phi3p[0]+gl2(g2, h)*phi3p[1])*h/2;
  }
  //右边界
  for(int j = 0; j < Nx; ++j) {
    phi2p[0] = (1+1/sqrt(3))/2; phi2p[1] = (1-1/sqrt(3))/2;
    phi3p[0] = (1-1/sqrt(3))/2; phi3p[1] = (1+1/sqrt(3))/2;
    x1[0] = cd[0][en[2][2*Nx*(j+1)-1]];
    x1[1] = cd[1][en[2][2*Nx*(j+1)-1]];
    x2[0] = cd[0][en[4][2*Nx*(j+1)-1]];
    x2[1] = cd[1][en[4][2*Nx*(j+1)-1]];
    g1 = -(x2-x1)/2./sqrt(3) + (x2+x1)/2.;
    g2 = (x2-x1)/2./sqrt(3) + (x2+x1)/2.;
    int m = 2*Nx + j*2*(Nx+1);
    b(m) += (gr1(g1, h)*phi2p[0]+gr1(g2, h)*phi2p[1])*h/2;
    m = 2*Nx + j*2*(Nx+1) + 1;
    b(m) += (gr2(g1, h)*phi2p[0]+gr2(g2, h)*phi2p[1])*h/2;
    m = 2*Nx + (j+1)*2*(Nx+1);
    b(m) += (gr1(g1, h)*phi3p[0]+gr1(g2, h)*phi3p[1])*h/2;
    m = 2*Nx + (j+1)*2*(Nx+1) + 1;
    b(m) += (gr2(g1, h)*phi3p[0]+gr2(g2, h)*phi3p[1])*h/2;
  }
  //上边界
  for(int i = 0; i < Nx; ++i) {
    phi2p[0] = (1-1/sqrt(3))/2; phi2p[1] = (1+1/sqrt(3))/2;
    phi3p[0] = (1+1/sqrt(3))/2; phi3p[1] = (1-1/sqrt(3))/2;
    x1[0] = cd[0][en[4][2*Nx*(Nx-1)+2*i]];
    x1[1] = cd[1][en[4][2*Nx*(Nx-1)+2*i]];
    x2[0] = cd[0][en[2][2*Nx*(Nx-1)+2*i]];
    x2[1] = cd[1][en[2][2*Nx*(Nx-1)+2*i]];
    g1 = -(x2-x1)/2./sqrt(3) + (x2+x1)/2.;
    g2 = (x2-x1)/2./sqrt(3) + (x2+x1)/2.;
    int m = 2*i + Nx*2*(Nx+1);
    b(m) += (gu1(g1, h)*phi3p[0]+gu1(g2, h)*phi3p[1])*h/2;
    m = 2*i + Nx*2*(Nx+1) + 1;
    b(m) += (gu2(g1, h)*phi3p[0]+gu2(g2, h)*phi3p[1])*h/2;
    m = 2*i + Nx*2*(Nx+1) + 2;
    b(m) += (gu1(g1, h)*phi2p[0]+gu1(g2, h)*phi2p[1])*h/2;
    m = 2*i + Nx*2*(Nx+1) + 3;
    b(m) += (gu2(g1, h)*phi2p[0]+gu2(g2, h)*phi2p[1])*h/2;
  }

}

void fem::DiricletBD(func u0, func v0) {
  //下边界,Diriclet, zero BD
  cor x0;
  for(int i = 0; i < Nx+1; ++i) {
    x0[0] = cd[0][i]; x0[1] = cd[1][i];
    int m = 2*i;
    double tmp = u0(x0, h);
    fix_entry(m, tmp);

    m = 2*i+1;
    tmp = v0(x0, h);
    fix_entry(m, tmp);
  }
}

void fem::fix_entry(const int m, const double tmp) {
  for(int j = 0; j < A.rows(); ++j) {
    if(j != m && A.coeff(j,m) != 0) {
      b(j) -= A.coeff(j,m)*tmp;
      A.coeffRef(j, m) = 0;
    }
    if(j != m && A.coeff(m,j) != 0) {
      A.coeffRef(m, j) = 0;
    }
  }
  //A.coeffRef(m, m) = 1;
  b(m) = tmp*A.coeffRef(m, m);
}

void fem::DiricletBD_full(func u0, func v0) {
  //四条边,Diriclet, zero BD
  cor x0;
  for(int i = 0; i < Nx+1; ++i) {
    int m = 2*i;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    double tmp = u0(x0, h);
    fix_entry(m, tmp);

    m = 2*i+1;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    tmp = v0(x0, h);
    fix_entry(m, tmp);

    m = 2*i+2*(Nx+1)*Nx;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    tmp = u0(x0, h);
    fix_entry(m, tmp);

    m = 2*i+2*(Nx+1)*Nx + 1;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    tmp = v0(x0, h);
    fix_entry(m, tmp);
  }
  for(int i = 1; i < Nx; ++i) {
    int m = 2*(Nx+1)*i;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    double tmp = u0(x0, h);
    fix_entry(m, tmp);

    m = 2*(Nx+1)*i + 1;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    tmp = v0(x0, h);
    fix_entry(m, tmp);

    m = 2*(Nx+1)*(i+1)-2;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    tmp = u0(x0, h);
    fix_entry(m, tmp);

    m = 2*(Nx+1)*(i+1)-1;
    x0[0] = cd[0][m]; x0[1] = cd[1][m];
    tmp = v0(x0, h);
    fix_entry(m, tmp);
  }

}


void fem::run() {
  //std::cout << "======solve_leqn by iterative solver======" << std::endl;
  solver.setMaxIterations(1e3);
  solver.setTolerance(1e-12);
  solver.compute(A);
  NumSol = solver.solve(b);
  std::cout << "iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;
  //std::cout << "======================" << std::endl;
}

void fem::cal_eps_sig() {
  double b1, b2, b3, c1, c2, c3;
  double D = pow(h,2);
  EMAT B(3,6), L(3,3);
  EVEC EPS(3), SIG(3), PHI(6);
  L(0,0) = lam+2*mu; L(0,1) = lam;
  L(1,0) = lam; L(1,1) = lam+2*mu; L(2,2) = 2*mu;
  for(int e = 0; e < Nx*Nx*2; ++e) {
    b1 = cd[1][en[2][e]] - cd[1][en[4][e]];
    c1 = cd[0][en[2][e]] - cd[0][en[4][e]];
    b2 = cd[1][en[4][e]] - cd[1][en[0][e]];
    c2 = cd[0][en[4][e]] - cd[0][en[0][e]];
    b3 = cd[1][en[0][e]] - cd[1][en[2][e]];
    c3 = cd[0][en[0][e]] - cd[0][en[2][e]];
    B(0,0) = b1; B(0,2) = b2; B(0,4) = b3;
    B(1,1) = -c1; B(1,3) = -c2; B(1,5) = -c3;
    B(2,0) = -c1/2; B(2,1) = b1/2; B(2,2) = -c2/2;
    B(2,3) = b2/2; B(2,4) = -c3/2; B(2,5) = b3/2;
    PHI[0] = NumSol(en[0][e]);
    PHI[1] = NumSol(en[1][e]);
    PHI[2] = NumSol(en[2][e]);
    PHI[3] = NumSol(en[3][e]);
    PHI[4] = NumSol(en[4][e]);
    PHI[5] = NumSol(en[5][e]);
    EPS = B*PHI/D;
    SIG = L*EPS;
    eps(e, 0) = EPS(0);
    eps(e, 1) = EPS(1);
    eps(e, 2) = EPS(2);
    sig(e, 0) = SIG(0);
    sig(e, 1) = SIG(1);
    sig(e, 2) = SIG(2);
  }
}

void fem::print(const std::string& filename, const std::string& filename1,
    func u0, func v0) {
  cor x;

  std::ofstream outfile(filename.c_str(),std::ios::out);
  if(!outfile){
    std::cerr<<"Wrong!"<<std::endl;
  }
  outfile.precision(16);
  outfile << std::showpos;
  outfile.setf(std::ios::scientific);
  std::ofstream outfile1(filename1.c_str(),std::ios::out);
  if(!outfile1){
    std::cerr<<"Wrong!"<<std::endl;
  }
  outfile1.precision(16);
  outfile1 << std::showpos;
  outfile1.setf(std::ios::scientific);

  for(int j = 0; j < Nx+1; ++j) {
    for(int i = 0; i < Nx+1; ++i) {
      x[0] = h*i;
      x[1] = h*j;
      outfile
        << NumSol[j*2*(Nx+1)+2*i] << " "
        << NumSol[j*2*(Nx+1)+2*i+1] << " "
        << u0(x, h) << " "
        << v0(x, h) << "\n";
    }
  }
  outfile << std::endl;
  outfile.close();

  for(int e = 0; e < Nx*Nx*2; ++e) {
    x[0] = (cd[0][en[0][e]]+cd[0][en[2][e]]+cd[0][en[4][e]])/3.;
    x[1] = (cd[1][en[0][e]]+cd[1][en[2][e]]+cd[1][en[4][e]])/3.;
    outfile1
      << eps(e, 0) << " "
      << eps(e, 1) << " "
      << eps(e, 2) << " "
      << sig(e, 0) << " "
      << sig(e, 1) << " "
      << sig(e, 2) << "\n";
  }
  outfile1 << std::endl;
  outfile1.close();

}

//void fem::print(const std::string& filename, const std::string& filename1,
    //func u0, func v0,
    //func ex, func ey, func exy, func sx, func sy, func sxy) {
  //cor x;

  //std::ofstream outfile(filename.c_str(),std::ios::out);
  //if(!outfile){
    //std::cerr<<"Wrong!"<<std::endl;
  //}
  //outfile.precision(16);
  //outfile << std::showpos;
  //outfile.setf(std::ios::scientific);
  //std::ofstream outfile1(filename1.c_str(),std::ios::out);
  //if(!outfile1){
    //std::cerr<<"Wrong!"<<std::endl;
  //}
  //outfile1.precision(16);
  //outfile1 << std::showpos;
  //outfile1.setf(std::ios::scientific);

  //for(int j = 0; j < Nx+1; ++j) {
    //for(int i = 0; i < Nx+1; ++i) {
      //x[0] = h*i;
      //x[1] = h*j;
      //outfile
        //<< NumSol[j*2*(Nx+1)+2*i] << " "
        //<< NumSol[j*2*(Nx+1)+2*i+1] << " "
        //<< u0(x) << " "
        //<< v0(x) << "\n";
    //}
  //}
  //outfile << std::endl;
  //outfile.close();

  //for(int e = 0; e < Nx*Nx*2; ++e) {
    //x[0] = (cd[0][en[0][e]]+cd[0][en[2][e]]+cd[0][en[4][e]])/3.;
    //x[1] = (cd[1][en[0][e]]+cd[1][en[2][e]]+cd[1][en[4][e]])/3.;
    //outfile1
      //<< eps(e, 0) << " "
      //<< eps(e, 1) << " "
      //<< eps(e, 2) << " "
      //<< sig(e, 0) << " "
      //<< sig(e, 1) << " "
      //<< sig(e, 2) << " "
      //<< ex(x) << " "
      //<< ey(x) << " "
      //<< exy(x) << " "
      //<< sx(x) << " "
      //<< sy(x) << " "
      //<< sxy(x) << "\n";
  //}
  //outfile1 << std::endl;
  //outfile1.close();

//}

void fem::print_geometry() {
  cor x;
  char str[100];
  sprintf(str, "P%d.dat", Nx);
  std::ofstream hehe(str,std::ios::out);
  if(!hehe){
    std::cerr<<"Wrong!"<<std::endl;
  }
  for(int j = 0; j < Nx+1; ++j) {
    for(int i = 0; i < Nx+1; ++i) {
      x[0] = h*i;
      x[1] = h*j;
      hehe << x[0] << " " << x[1] << "\n";
    }
  }
  hehe << std::endl;
  hehe.close();

  sprintf(str, "TRI%d.dat", Nx);
  std::ofstream haha(str,std::ios::out);
  if(!haha){
    std::cerr<<"Wrong!"<<std::endl;
  }
  for(int e = 0; e < Nx*Nx*2; ++e) {
    haha << en[0][e]/2+1 << " "
      << en[2][e]/2+1 << " "
      << en[4][e]/2+1 << "\n";
  }
  haha << std::endl;
  haha.close();
}

