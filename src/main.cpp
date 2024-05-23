#include <iostream>
#include <vector>
#include "writeVTK.hpp"
#include <cmath>
#include <numbers>
#include "Jacobi.hpp"

int main(){ 
    
  using namespace JacobiMethod;

  constexpr double pi = std::numbers::pi;
  int NX, NY;
  NX = 100;
  NY = 100;

  double hx, hy;
  hx = 1.0/ (NX+1);
  hy = 1.0/ (NY+1);

  std::vector<std::vector<double>> exacSol;

  double x,y;
  auto  fun = [](const double & x, const double & y){ return std::sin(2*pi*x)*std::sin(2*pi*y);};

  for (int i = 0; i < NX+1; i++) {
    std::vector<double> row;
    x = i*hx;
    for (int j = 0; j < NY+1; j++) {
      y = j*hy;
      row.push_back(fun(x, y));
    }
    exacSol.push_back(row);
  }
  std::cout << "Generating VTK file..."<<std::flush ;
  generateVTKFile("results/exact.vtk", exacSol, NX,NY, hx, hy);
  std::cout << "Done!" << std::endl;

  ElemType<double> sol;
  auto rhs = [](const double & x, const double & y){ return 8*pi*pi*std::sin(2*pi*x)*std::sin(2*pi*y);};
  Jacobi<double> solver(rhs, NX, hx, 1e-4, 1000);
  std::cout << "Solving the problem..."<<std::flush ;
  solver.solve(sol);
  std::cout << "Done!" << std::endl;

  std::vector<std::vector<double>> numSol;
  for (int i = 0; i < NX+1; i++) {
    std::vector<double> row;
    x = i*hx;
    for (int j = 0; j < NY+1; j++) {
      y = j*hy;
      row.push_back(sol[{i,j}]);
    }
    numSol.push_back(row);
  }
  std::cout << "Generating VTK file..."<<std::flush ;
  generateVTKFile("results/numerical.vtk", numSol, NX,NY, hx, hy);
  std::cout << "Done!" << std::endl;
  return 0;
}
