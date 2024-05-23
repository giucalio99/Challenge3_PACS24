#ifndef HH_JACOBI_HH
#define HH_JACOBI_HH
#include <map>
#include <array>
#include <functional>
namespace JacobiMethod{

template <class T>
    using ElemType = std::map<std::array<std::size_t,2>,T>;

template <class T>
class Jacobi {
    private: 
        std::function<T(T,T)> m_fun;
        int m_n;
        T m_h;
        T m_tol;
        int m_maxIter;
        ElemType<T> m_U;
    public:
        Jacobi(const std::function<T(T,T)>& fun, int n, T h,T tol, int maxIter);
        void solve(ElemType<T>& sol);
};

#include "Jacobi_impl.hpp"
}//namespace JacobiMethod
#endif//HH_JACOBI_HH

