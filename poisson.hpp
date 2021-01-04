///
/// \file include/fista/poisson.hpp
/// \brief FISTA (Fast Iterative Shrinkage Tresholding Algorithm) solver for Poisson distributed noise.
/// \author Hatef Monajemi <monajemi@stanford.edu> 2012-2014
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_FISTA_POISSON_HPP
#define ASTROQUT_FISTA_POISSON_HPP

#include "utils/linearop/matrix.hpp"
#include "utils/linearop/operator/matmult.hpp"
#include "utils/linearop/operator/convolution.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>

namespace alias
{
namespace fista
{
namespace poisson
{

template<class T = double>
struct Parameters
{
    /** Default constructor
     *  Parameters with default values
     */
    Parameters() noexcept
        : tol(1e-6)
        , iter_max(1000)
        , init_value{}
        , indices(Matrix<size_t>(0,1,1))
        , log(true)
        , log_period(20)
    {}

    T tol; //!< Member variable "tol"
    size_t iter_max; //!< Member variable "iter_max"
    Matrix<T> init_value; //!< Member variable "init_value"
    Matrix<size_t> indices; //!< Member variable "indices"
    bool log; //!< Member variable "log"
    unsigned int log_period; //!< Member variable "log_period"
};

/** Poisson distributed noise solver
 *  \param A Explicit regression matrix
 *  \param u Background shift
 *  \param b Response data to the regression matrix
 *  \param lambda Regularization parameter
 *  \param options Parameters that defines various value for FISTA to work
 */
template<class T>
T Func(const Matrix<T>& Axu,
       const Matrix<T>& b )
{
    // sum(A*x+u - b.*log(A*x+u))
    return (Axu - (b & (Axu.Log()))).Sum();
}
template<class T>
Matrix<T> FuncGrad(const Matrix<T>& Axu,
                   const Operator<T>& At,
                   const Matrix<T>& b )
{
    // A' * ((A*x+u - b) ./ (A*x+u))
    return At * ((Axu - b) / Axu);
}
template<class T>
T FLasso(const Matrix<T>& Axu,
         const Matrix<T>& x_woi,
         const Matrix<T>& b,
         T lambda )
{
    // sum(A*x+u - b.*log(A*x+u)) + lambda * norm(x[-0],1)
    return Func(Axu, b) + lambda*x_woi.Norm(one);
}
template<class T>
T FLassoApprox(const Matrix<T>& Ayu,
               const Operator<T>& At,
               const Matrix<T>& x,
               const Matrix<T>& x_woi,
               const Matrix<T>& y,
               const Matrix<T>& b,
               T lambda,
               T L )
{
    // sum(A*x+u - b.*log(A*x+u)) + <(x-y), gradfunc(A,y,u,b,w)> + 0.5*L*||x-y||^2  + lambda * norm(x[-0],1)
    return Func(Ayu, b) + Inner( x-y, FuncGrad(Ayu, At, b) ) + 0.5*L*(x-y).Norm(two_squared) + lambda*x_woi.Norm(one);
}
template<class T>
Matrix<T> Solve(const Operator<T>& A,
                const Matrix<T>& u,
                const Matrix<T>& b,
                T lambda,
                const Parameters<T>& options )
{
    std::cout << std::defaultfloat;
    std::cout << std::string(37, '*') << " FISTA " << std::string(36, '*') << std::endl;
    std::cout << "A: " << A.Height() << "x" << A.Width() << " matrix";
    std::cout << ", u: " << u.Length() << " vector";
    std::cout << ", b: " << b.Length() << " vector" << std::endl;
    std::cout << "lambda:" << lambda << ", tol:" << options.tol << std::endl;
    std::cout << std::string(80, '*') << std::endl << std::endl;
    std::cout << " iter" << " | " << "         tol        " << " | " << "       FLasso       " << " | " << "     Lf      " << " | " << "  NNZ   " << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::scientific;

    // x and y variables
    Matrix<T> x((T)0, A.Width(), 1);
    if( !options.init_value.IsEmpty() )
        x = options.init_value;
    Matrix<T> x_next(x);
    Matrix<T> x_next_woi(x_next.Data()+1, x_next.Height()-1, 1); // points to second element of x_next;
    Matrix<T> y(x);

    // intermediate results
    Matrix<T> Axu = A*x+u;
    Matrix<T> Ax_nextu;
    Matrix<T> Ayu = Axu;
    Operator<T> *A_copy = A.Clone();
    Operator<T> &At = A_copy->Transpose();
    T f_lasso_next = (T)0;
    T f_lasso_previous[10] {};
    f_lasso_previous[0] = FLasso(Axu, x_next_woi, b, lambda);
    Matrix<T> grad_current = FuncGrad(Axu, At, b);

    // FISTA variables
    T tol = std::numeric_limits<T>::infinity();
    T Lf = (T)1;
    T eta = (T)2;
    T L_bar = (T)0;
#ifdef CLASSIC_FISTA
    T t = 1.0L;
    T t_next = (1.0L + std::sqrt(1.0L + 4.0L * t * t)) / 2.0L;
#endif // CLASSIC_FISTA
    size_t k = 0;

    // main loop
    while( std::abs(tol) > options.tol && k < options.iter_max )
    {
        // backtracking loop
        T beta = std::numeric_limits<T>::infinity();
        for( int ik = 0; beta > 0; ++ik )
        {
            L_bar = std::pow(eta, ik) * Lf;
            x_next = y - (grad_current/L_bar);
            x_next_woi.Data(x_next.Data()+1); // points to second element of new x_next
            std::move(x_next_woi).Shrink(lambda/L_bar); //cast to an rvalue to allow in-place shrinkage
            std::move(x_next).RemoveNeg(options.indices);
            Ax_nextu = (A*x_next)+u;
            if( Ax_nextu.ContainsNeg() ) // skip function evaluation if we have negative values
                continue;
            f_lasso_next = FLasso(Ax_nextu, x_next_woi, b, lambda);
            beta = f_lasso_next - FLassoApprox(Ayu, At, x_next, x_next_woi, y, b, lambda, L_bar);
        }

        // FISTA step
#ifdef CLASSIC_FISTA
        y = x_next + (x_next - x) * ((t - 1.0)/t_next);
#else
        y = x_next;
#endif // CLASSIC_FISTA

        // compute tol from previous function value
        T f_lasso_previous_sum = std::accumulate(f_lasso_previous, f_lasso_previous+10, (T)0) / std::min((T) k+1, (T)10);
        tol = std::abs( f_lasso_next - f_lasso_previous_sum ) / f_lasso_previous_sum;

        // actualize values for next iteration
        ++k;
        x = x_next;
        Axu = Ax_nextu;
        Ayu = A*y+u;
        f_lasso_previous[k % 10] = f_lasso_next;
        grad_current = FuncGrad(Axu, At, b);
#ifdef CLASSIC_FISTA
        t = t_next;
        t_next = (1.0L + std::sqrt(1.0L + 4.0L * t * t)) / 2.0L;
#else
        Lf = (k % 100 == 0 ? (T)1 : L_bar / (T)2);
#endif // CLASSIC_FISTA

        if( options.log )
        {
            if( k % (options.log_period * 50) == 0 )
            {
                std::cout << std::endl << " iter" << " | " << "         tol        " << " | " << "       FLasso       " << " | " << "     Lf      " << " | " << "  NNZ  " << std::endl;
                std::cout << std::string(80, '-') << std::endl;
            }
            if( k % options.log_period == 0 )
            {
                std::cout << std::setw(5) << k << " | " << std::scientific << std::setprecision(10) << std::setw(20) << tol << " | " << std::setw(20) << f_lasso_next << " | " << std::defaultfloat << std::setw(13) << Lf << " | " << std::setw(8) << x.NonZeroAmount() << std::endl;
            }
        }
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::setw(5) << k << " | " << std::scientific << std::setprecision(10) << std::setw(20) << std::abs(tol) << " | " << std::setw(20) << f_lasso_next << " | " << std::defaultfloat << std::setw(13) << Lf << " | " << std::setw(8) << x.NonZeroAmount() << std::endl << std::endl << std::endl;

    if(k < options.iter_max)
        std::cout << "FISTA: converged in " << k << " iterations" << std::endl;
    else
        std::cout << "FISTA: did not converge after " << k << " iterations" << std::endl;

    std::cout << "FISTA: Relative error: " << std::abs(tol) << std::endl << std::endl;

    x_next_woi.Data(nullptr); // release pointer
    delete A_copy;

    return x;
}

} // namespace poisson
} // namespace fista
} // namespace alias

#endif // ASTROQUT_FISTA_POISSON_HPP

