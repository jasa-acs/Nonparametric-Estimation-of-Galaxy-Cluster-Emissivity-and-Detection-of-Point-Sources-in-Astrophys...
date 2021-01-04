///
/// \file src/test/matrix.cpp
/// \brief Test suite specializations to validate the Matrix class with particular types
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#include "test/matrix.hpp"

namespace alias
{
namespace test
{
namespace matrix
{

template <>
const Matrix<std::complex<double>> SquareMatrix()
{
    std::complex<double> data[9] = {{1,2},{3,4},{5,6},{7,8},{9,10},{11,12},{13,14},{15,16},{17,18}};
    return Matrix<std::complex<double>>(data, 9, 3, 3);
}

template <>
const Matrix<std::complex<double>> RectMatrix()
{
    std::complex<double> data[10] = {{1,2},{3,4},{5,6},{7,8},{9,10},{11,12},{13,14},{15,16},{17,18},{19,20}};
    return Matrix<std::complex<double>>(data, 10, 2, 5);
}

template <>
bool TransposeSquare<std::complex<double>>()
{
    std::cout << "Transpose test with complex<double> square matrices : ";

    std::complex<double> data[9] = {{1,2},{7,8},{13,14},{3,4},{9,10},{15,16},{5,6},{11,12},{17,18}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result = (Compare(expected_result, SquareMatrix<std::complex<double>>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool TransposeRect<std::complex<double>>()
{
    std::cout << "Transpose test with complex<double> rect matrices: ";

    std::complex<double> data[10] = {{1,2},{11,12},{3,4},{13,14},{5,6},{15,16},{7,8},{17,18},{9,10},{19,20}};
    const Matrix<std::complex<double>> expected_result(data, 10, 5, 2);
    bool test_result =  (Compare(expected_result, RectMatrix<std::complex<double>>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Add<std::complex<double>>()
{
    std::cout << "Add test with complex<double> : ";

    std::complex<double>data[9] = {{2,4},{6,8},{10,12},{14,16},{18,20},{22,24},{26,28},{30,32},{34,36}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<std::complex<double>>() + SquareMatrix<std::complex<double>>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Sub<std::complex<double>>()
{
    std::cout << "Sub test with complex<double> : ";

    std::complex<double> data[9] = {{-1,-2},{-3,-4},{-5,-6},{-7,-8},{-9,-10},{-11,-12},{-13,-14},{-15,-16},{-17,-18}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<std::complex<double>>() - std::complex(2.0,0.0)*SquareMatrix<std::complex<double>>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultSquare<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> square matrix : ";

    std::complex<double> data[9] = {{-33,204},{-39,246},{-45,288},{-51,474},{-57,588},{-63,702},{-69,744},{-75,930},{-81,1116}};
    const Matrix<std::complex<double>> expected_result(data, 9, 3, 3);
    bool test_result =  (Compare(expected_result, SquareMatrix<std::complex<double>>() * SquareMatrix<std::complex<double>>()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultRect<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> rect matrix: ";

    std::complex<double> data[4] = {{-55,380},{-105,930},{-105,930},{-155,2480}};
    const Matrix<std::complex<double>> expected_result(data, 4, 2, 2);
    bool test_result =  (Compare(expected_result, RectMatrix<std::complex<double>>() * RectMatrix<std::complex<double>>().Transpose()));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultVectMat<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> vector times matrix : ";

    std::complex<double> data_vect[2] = {{1,2},{3,4}};
    const Matrix<std::complex<double>> vect(data_vect, 2, 1, 2);

    std::complex<double> data[5] = {{-18,84},{-22,104},{-26,124},{-30,144},{-34,164}};
    const Matrix<std::complex<double>> expected_result(data, 5, 1, 5);

    bool test_result = Compare(expected_result, vect * RectMatrix<std::complex<double>>());
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool MultMatVect<std::complex<double>>()
{
    std::cout << "Mult test with complex<double> matrix times vector : ";

    std::complex<double> data_vect[5] = {{1,2},{3,4},{5,6},{7,8},{9,10}};
    const Matrix<std::complex<double>> vect(data_vect, 5, 5, 1);

    std::complex<double> data[2] = {{-55,380},{-105,930}};
    const Matrix<std::complex<double>> expected_result(data, 2, 2, 1);

    bool test_result =  Compare(expected_result, RectMatrix<std::complex<double>>() * vect);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool NormOne<std::complex<double>>()
{
    std::cout << "Norm one test with complex<double> : ";

    bool test_result =  IsEqual(SquareMatrix<std::complex<double>>().Norm(one), 121.2044302895985);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool NormTwo<std::complex<double>>()
{
    std::cout << "Norm two test with complex<double> : ";

    bool test_result =  IsEqual(SquareMatrix<std::complex<double>>().Norm(two), 45.923850012820139);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool NormInf<std::complex<double>>()
{
    std::cout << "Norm inf test with complex<double> : ";

    bool test_result = IsEqual(SquareMatrix<std::complex<double>>().Norm(inf), 24.758836806279895);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Sum<std::complex<double>>()
{
    std::cout << "Sum test with complex<double> : ";

    bool test_result = IsEqual(SquareMatrix<std::complex<double>>().Sum(), std::complex(81.0,90.0));
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

template <>
bool Shrink<std::complex<double>>()
{
    std::cout << "Shrinkage test with complex<double> is not implemented" << std::endl;

    return true;
}

bool Input()
{
    std::cout << "Input test : ";

    double data[64] = {0.718326887014344, 0.254697476919073, 0.00927382750852335, 0.577122634981908, 0.633734394010610, 0.421283139598629, 0.460202078925928, 0.674443684458196, 0.986048510786069, 0.748893933718032, 0.374450622100372, 0.120330608524745, 0.840897872853126, 0.563880064958332, 0.948640666270834, 0.0205356939531407, 0.961849177962106, 0.950107798625595, 0.560196923165714, 0.752753955796127, 0.971311449868029, 0.401363407979538, 0.897968197861966, 0.863428619539009, 0.885987805788211, 0.0354932749668265, 0.296377093786706, 0.00436407800590133, 0.511142964607578, 0.865052190242670, 0.0515327248096797, 0.612173298525157, 0.128803848149846, 0.571319022444656, 0.721881326171109, 0.375275208475366, 0.533503174291409, 0.457215880934024, 0.561481964909078, 0.0836644349582446, 0.626739647975978, 0.234699749896591, 0.573726220563258, 0.517413402072890, 0.793261259811798, 0.0624005653201738, 0.241007512948713, 0.521477853178805, 0.819861713926915, 0.292257563646496, 0.0719880830244852, 0.253174979246155, 0.796468203492791, 0.642068697604116, 0.0314088549147913, 0.482489712294326, 0.0977842049266892, 0.795171287382555, 0.880957791911562, 0.218089090918185, 0.549041240475796, 0.347361789599811, 0.293446717464101, 0.470367130245693};
    const Matrix<double> expected_result(data, 64, 8, 8);
#ifdef VERBOSE
    std::cout << std::endl << "Expected result :" << expected_result;
#endif // VERBOSE

    Matrix<double> input_test("data/test/test.data", 8, 8);

#ifdef VERBOSE
    std::cout << std::endl << "Read data :" << input_test;
#endif // VERBOSE

    bool test_result =  Compare(expected_result, input_test);
    std::cout << (test_result ? "Success" : "Failure") << std::endl;

    return test_result;
}

} // namespace matrix
} // namespace test
} // namespace alias
