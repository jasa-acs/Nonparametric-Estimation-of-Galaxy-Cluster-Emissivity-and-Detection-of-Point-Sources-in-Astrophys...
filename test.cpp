///
/// \file src/test.cpp
/// \brief Implementation of the test suites.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#include "test.hpp"

namespace alias
{
namespace test
{

bool OperatorTest()
{
    using namespace oper;

    bool convolution = ConvolutionTest();

//    ConvolutionTime(1024, 5);

    bool abel_build = AbelTestBuild();
    bool abel_apply = AbelTestApply();
    bool abel_apply2 = AbelTestApply2();
    bool abel_transposed = AbelTestTransposed();
    bool abel_transposed2 = AbelTestTransposed2();

//    AbelTime(512);

    bool wavelet = WaveletTest();
    bool wavelet2 = WaveletTest2();
    bool wavelet3 = WaveletTest3();

    bool spline = SplineTest();

    bool blur = BlurTest();

    bool astro = AstroTest();
    bool astro_transposed = AstroTestTransposed();

    return convolution && abel_build && abel_apply && abel_apply2 && abel_transposed && abel_transposed2 && wavelet && wavelet2 && wavelet3 && spline && blur && astro && astro_transposed;
}

bool FISTATest()
{
    bool fista_small = fista::SmallExample();

//    fista::Time(1024);

    return fista_small;
}

} // namespace test
} // namespace alias
