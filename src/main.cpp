#include <iostream>
#include <vector>
#include "writeVTK.hpp"
#include <cmath>
#include <numbers>
#include <mpi.h>
#include <omp.h>

#include "Jacobi.hpp"

int main(int argc, char** argv){ 
    
  using namespace JacobiMethod;
  if (argc != 2) {
    std::cerr << "Usage of" << argv[0]
              << " number-of-threads\n";
    std::exit(EXIT_FAILURE);
  }
  const auto nthreads = std::atoi(argv[1]);

  // initialize MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0) {
    std::cout << "Provided parallelization lvl: " << provided << std::endl;
    std::cout << "# OpenMP threads: " << nthreads << std::endl;
    std::cout << "# MPI procs: " << size << std::endl;
  }
 
  constexpr double pi = std::numbers::pi;
  int N;
  N = 100;


  double h;
  h = 1.0/ (N+1);
  //Computing the exact solution by master processor
  if( rank==0){
    std::vector<std::vector<double>> exacSol;

    double x,y;
    auto  fun = [](const double & x, const double & y){ return std::sin(2*pi*x)*std::sin(2*pi*y);};

    for (int i = 0; i < N+1; i++) {
      std::vector<double> row;
      x = i*h;
      for (int j = 0; j < N+1; j++) {
        y = j*h;
        row.push_back(fun(x, y));
      }
      exacSol.push_back(row);
    }
    std::cout << "Generating VTK file..."<<std::flush ;
    generateVTKFile("results/exact.vtk", exacSol, N,N, h, h);
    std::cout << "Done!" << std::endl;
  }
  //wait for the master processor to finish generating the exact solution
  MPI_Barrier(MPI_COMM_WORLD);
  auto rhs = [](const double & x, const double & y){ return 8*pi*pi*std::sin(2*pi*x)*std::sin(2*pi*y);};
  
  ElemType<double> sol, local_sol;
  //const auto N_local= N/size ;// assume no remainder for simplicity
  Jacobi<double> solver(rhs, N,N, h, 1e-4, 1000);
  std::cout << "Solving the problem..."<<std::flush ;
  solver.solve(sol, nthreads);
  std::cout << "Done!" << std::endl;

  std::vector<std::vector<double>> numSol;
  for (int i = 0; i < N+1; i++) {
    std::vector<double> row;
    for (int j = 0; j < N+1; j++) {
      row.push_back(sol[{i,j}]);
    }
    numSol.push_back(row);
  }
  std::cout << "Generating VTK file..."<<std::flush ;
  generateVTKFile("results/numerical.vtk", numSol, N,N, h, h);
  std::cout << "Done!" << std::endl;
  MPI_Finalize();
  return 0;
}
