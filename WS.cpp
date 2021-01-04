///
/// \file src/test/WS.cpp
/// \brief Test suite to test the WS method.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#include "test/WS.hpp"

namespace alias
{
namespace test
{
namespace astro
{

void Chandra()
{
    std::cout << std::endl << std::endl << "WS test with 512x512 Chandra data." << std::endl;

    Matrix<double> expected_result(std::string("data/512_chandra/MATLAB/2000/fhat.data"), 512, 1);

    WS::Parameters<double> options;
    options.resample_windows_size = 4;
    options.pic_size = 512;
    options.MC_max = 1000;
    options.blurring_filter = "data/512_chandra/blurring.data";
    options.fista_params.iter_max = 1000;
    options.fista_params.log_period = 1;

    Matrix<double> solution = WS::Solve("data/512_chandra/F.data",
                                        "data/512_chandra/E.data",
                                        "data/512_chandra/O.data",
                                        "data/512_chandra/AstroQUT/" + std::to_string(options.MC_max) + "/solution.data",
                                        options);

    std::cout << std::endl;
    std::cout << "MATLAB result : " << std::endl;
    std::cout << "FISTA: did not converge in 1000 iterations (nnz(x)= 2329)" << std::endl;
    std::cout << "FISTA: Relative error: -1.93e-06" << std::endl;
    std::cout << "FISTA: Final FLasso value: -1.387380e+06" << std::endl;
    std::cout << "Elapsed time is 9293.529344 seconds." << std::endl << std::endl;
}

void Sym256()
{
    std::cout << std::endl << std::endl << "Test with 256x256 simulated data." << std::endl;

    WS::Parameters<double> options;
    options.bootstrap_max = 1;
    options.resample_windows_size = 4;
    options.pic_size = 256;
    options.MC_max = 5000;
    options.blurring_filter = "data/512_chandra/blurring.data";
    options.fista_params.iter_max = 1000;

    for(size_t i = 1; i <=24; ++i)
    {
        WS::Solve("data/256/M256blockssymPS1B1W1S1Erealnk24forcepos2/result" + std::to_string(i) + ".data",
                  "data/256/E.data",
                  "data/256/O.data",
                  "data/256/M256blockssymPS1B1W1S1Erealnk24forcepos2/solution" + std::to_string(i) + ".data",
                  options);
    }

    for(size_t i = 1; i <=24; ++i)
    {
        WS::Solve("data/256/M256f4PS1B1W1S1Erealnk24forcepos2/result" + std::to_string(i) + ".data",
                  "data/256/E.data",
                  "data/256/O.data",
                  "data/256/M256f4PS1B1W1S1Erealnk24forcepos2/solution" + std::to_string(i) + ".data",
                  options);
    }

    for(size_t i = 1; i <=24; ++i)
    {
        WS::Solve("data/256/M256f4symPS1B1W1S1Erealnk24forcepos2/result" + std::to_string(i) + ".data",
                  "data/256/E.data",
                  "data/256/O.data",
                  "data/256/M256f4symPS1B1W1S1Erealnk24forcepos2/solution" + std::to_string(i) + ".data",
                  options);
    }

    for(size_t i = 1; i <=24; ++i)
    {
        WS::Solve("data/256/M256truefPS1B1W1S1Erealnk24forcepos2/result" + std::to_string(i) + ".data",
                  "data/256/E.data",
                  "data/256/O.data",
                  "data/256/M256truefPS1B1W1S1Erealnk24forcepos2/solution" + std::to_string(i) + ".data",
                  options);
    }

}

} // namespace astro
} // namespace test
} // namespace alias
