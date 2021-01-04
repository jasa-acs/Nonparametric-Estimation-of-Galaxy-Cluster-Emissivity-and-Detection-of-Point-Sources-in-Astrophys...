///
/// \file include/WS/astroQUT.hpp
/// \brief AstroQUT solver header
/// \author Jairo Diaz <jairo.diaz@unige.ch> MATLAB version 2016-2017
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_WS_ASTROQUT_HPP
#define ASTROQUT_WS_ASTROQUT_HPP

#include "fista/poisson.hpp"
#include "utils/linearop/matrix.hpp"

namespace alias
{
namespace WS
{

template<class T = double>
struct Parameters
{
    /** Default constructor
     *  Parameters with default values
     */
    Parameters()
        : pic_size(0)
        , model_size(0)
        , blurring_filter(std::string("data/blurring.data"))
        , bootstrap_max(1)
        , wavelet{3,8}
        , resample_windows_size(4)
        , MC_max(1000)
        , MC_quantile_PF(800)
        , MC_quantile_PS(999)
        , beta0(std::numeric_limits<double>::infinity())
        , lambda(1.0)
        , standardize{}
        , fista_params{}
    {
#ifdef DEBUG
        std::cerr << "WS::Parameters Default constructor called." << std::endl;
#endif // DEBUG
    }

    /** Default destructor
     */
    ~Parameters()
    {
#ifdef DEBUG
        std::cerr << "WS::Parameters Default destructor called." << std::endl;
#endif // DEBUG
    }

    size_t pic_size; //!< Member variable "pic_size" side of picture in pixel
    size_t model_size; //!< Member variable "model_size" wavelet + spline + point sources
    std::string blurring_filter; //!< Member variable "blurring_filter" path to the blurring filter data file
    size_t bootstrap_max; //!< Member variable "bootstrap_max" total amount of bootstraps computations, 0 for no bootstrapping
    size_t wavelet[2]; //!< Member variable "wavelet" wavelet type and wavelet parameter
    size_t resample_windows_size; //!< Member variable "resample_windows_size" side size of the resampling square
    size_t MC_max; //!< Member variable "MC_max" amount of Monte Carlo simulations to perform
    size_t MC_quantile_PF; //!< Member variable "MC_quantile_PF" quantile for lambda
    size_t MC_quantile_PS; //!< Member variable "MC_quantile_PS" quantile for lambdaI
    T beta0; //!< Member variable "beta0" intercept
    T lambda; //!< Member variable "lambda" first regularization parameter
    Matrix<T> standardize; //!< Member variable "standardize" standardisation matrix
    fista::poisson::Parameters<T> fista_params; //!< Member variable "fista_params" parameters to be given to the FISTA solver
};

Matrix<double> Solve(std::string picture_path,
                     std::string sensitivity_path,
                     std::string background_path,
                     std::string solution_path,
                     Parameters<double>& options );

} // namespace WS
} // namespace alias

#endif // ASTROQUT_WS_ASTROQUT_HPP

