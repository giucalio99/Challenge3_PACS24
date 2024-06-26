#ifndef HH_JACOBI_IMPL_HH
#define HH_JACOBI_IMPL_HH
#include <omp.h>
#include "Jacobi.hpp"

using namespace JacobiMethod;

template <class T>
Jacobi<T>::Jacobi(const std::function<T(T,T)>& fun, int nx, int ny, T h,T tol, int maxIter): 
    m_fun(fun), 
    m_nx(nx),
    m_ny(ny), 
    m_h(h),
    m_tol(tol),
    m_maxIter(maxIter){
    //initialize the solution with zeros
    //account for the boundary conditions of homogeous Dirichlet type
    for (int i = 0; i < nx+1; i++) {
        for (int j = 0; j < ny+1; j++) {
            m_U[{i,j}] = 0.0;
        }
    }
}

template <class T>
void Jacobi<T>::solve(ElemType<T> &sol, unsigned int num_threads){
    int iter = 0;  //set the iteration counter to zero
    T error = 1.0; //initialize the error to a value greater than the tolerance
    ElemType<T> U_old = m_U; //initialize the old solution to the initial guess
    while (iter < m_maxIter && error > m_tol){
        error = 0.0; //reset the error
        #pragma omp parallel for num_threads(num_threads) shared(m_nx,m_ny, iter, U_old)   \
        reduction(+:error)
        for (auto i = 1; i < m_nx; i++) {
            for (auto j = 1; j < m_ny; j++) {
                //update the solution with the Jacobi iteration
                m_U[{i,j}] = 0.25*(U_old[{i-1,j}] + U_old[{i+1,j}] + U_old[{i,j-1}] + U_old[{i,j+1}] + m_h*m_h*m_fun(i*m_h,j*m_h));
                //compute the error
                error +=(m_U[{i,j}]-U_old[{i,j}]) * (m_U[{i,j}]-U_old[{i,j}]);
                }
        }
        error= std::sqrt(m_h*error);//compute the norm of the error
        U_old = m_U;    //update the old solution
        iter++;//increment the iteration counter
        
    }
    sol = m_U;//return the solution
}
#endif//HH_JACOBI_IMPL_HH