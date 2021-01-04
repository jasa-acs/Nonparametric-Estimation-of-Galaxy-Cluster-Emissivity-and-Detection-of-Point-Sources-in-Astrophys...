///
/// \file include/utils/linearop/operator/matmult/spline.hpp
/// \brief Spline operator class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_SPLINE_HPP
#define ASTROQUT_UTILS_OPERATOR_SPLINE_HPP

#include "utils/linearop/operator/matmult.hpp"

namespace alias
{

template<class T = double>
class Spline : public MatMult<T>
{
private:

public:

    /** Default constructor
     */
    Spline()
    {
#ifdef DEBUG
        std::cerr << "Spline : Default constructor called." << std::endl;
#endif // DEBUG
    }

    /** Copy constructor
     *  \param other Object to copy from
     */
    Spline(const Spline& other)
        : MatMult<T>(other)
    {
#ifdef DEBUG
        std::cout << "Spline : Copy constructor called" << std::endl;
#endif // DEBUG
    }

    /** Move constructor
     *  \param other Object to move from
     */
    Spline(Spline&& other)
        : Spline()
    {
#ifdef DEBUG
        std::cout << "Spline : Move constructor called" << std::endl;
#endif // DEBUG
        swap(*this, other);
    }

    /** Full member constructor
     *  \param low_pass_filter QMF matrix for low pass filtering
     *  \param low_pass_filter Mirrored QMF matrix for high pass filtering
     *  \param height Height of the full Abel matrix
     *  \param width Width of the full Abel matrix
     */
    explicit Spline(Matrix<T> data, size_t height, size_t width)
        : MatMult<T>(data, height, width)
    {
#ifdef DEBUG
        std::cerr << "Spline : Full member constructor called." << std::endl;
#endif // DEBUG
    }

    /** Build constructor
     *  \brief Builds the Spline operator with the qmf matrix corresponding to type and parameter
     *  \param pic_size Side size of the picture in pixel
     */
    explicit Spline(size_t pic_size)
        : MatMult<T>(Generate(pic_size), pic_size, pic_size)
    {
#ifdef DEBUG
        std::cerr << "Spline : Build constructor called." << std::endl;
#endif // DEBUG
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    Spline* Clone() const override final
    {
        return new Spline(*this);
    }

    /** Default destructor
     */
    virtual ~Spline()
    {
#ifdef DEBUG
        std::cout << "Spline : Destructor called" << std::endl;
#endif // DEBUG
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    Spline& Transpose() override final
    {
        std::swap(this->height_, this->width_);
        this->transposed_ = !this->transposed_;
        this->data_ = std::move(this->data_.Transpose());
        return *this;
    }

    /** Swap function
     *  \param first First object to swap
     *  \param second Second object to swap
     */
    friend void swap(Spline& first, Spline& second) noexcept
    {
        using std::swap;

        swap(static_cast<Operator<T>&>(first), static_cast<Operator<T>&>(second));
    }

    /** Copy assignment operator
     *  \param other Object to assign to current object
     *  \return A reference to this
     */
    Spline& operator=(Spline other)
    {
        swap(*this, other);

        return *this;
    }

    /** Generate a full King spline matrix
     *  \brief Builds an Abel transform matrix with diagonal radius, without duplicating data or inserting zeros
     *  \param pic_side Width of the square picture
     */
    Matrix<T> Generate( size_t pic_side )
    {
        size_t height = 2*pic_side;
        size_t height_sqrt = (size_t) std::sqrt(height);

        Matrix<T> temp((T)0, height, pic_side);
        Matrix<T> result(pic_side, pic_side);
        T* bi = new T[height_sqrt];
        T* ri = new T[height_sqrt];

        for( size_t i = 0; i < height_sqrt; ++i )
            bi[i] = (T)1/(T)6 + (T)i*((T)10 - (T)1 / (T)6) / ((T)height_sqrt - (T)1);

        for( size_t i = 0; i < height_sqrt; ++i )
            ri[i] = (T)1 + (T)i*((T)pic_side / (T)2 - (T)1) / ((T)height_sqrt - (T)1);

        size_t row = 0;
        for( size_t subrow1 = 0; subrow1 < height_sqrt; ++subrow1 )
            for( size_t subrow2 = 0; subrow2 < height_sqrt; ++subrow2 )
            {
                for( size_t col = 0; col < pic_side/2; ++col )
                {
                    T local_result = std::pow((T)1 + std::pow( ((T)col + (T)1) / ri[subrow2], (T)2), (T)-3.0 * bi[subrow1]);
                    temp[row*pic_side + pic_side/2 + col] = local_result;
                    temp[row*pic_side + pic_side/2 - col - 1] = local_result;
                }
                T min_value = *std::min_element(&temp[row*pic_side], &temp[(row+1)*pic_side]);
                T accumulator = (T)0;
                for( size_t i = 0; i < pic_side; ++i )
                {
                    temp[row*pic_side + i] -= min_value;
                    accumulator += temp[row*pic_side + i];
                }
                for( size_t i = 0; i < pic_side; ++i )
                    temp[row*pic_side + i] /= accumulator;
                ++row;
            }

        for( size_t i = 0; i < pic_side/2; ++i )
            for( size_t j = pic_side/2; j < pic_side; ++j )
                temp[i*pic_side + j] = (T)0;

        for( size_t i = pic_side/2; i < height; ++i )
            for( size_t j = 0; j < pic_side/2; ++j )
                temp[i*pic_side + j] = (T)0;

        for( size_t i = 0; i < pic_side; ++i )
        {
            size_t index = std::floor((T)i*((T)height_sqrt*(T)height_sqrt-(T)1)/((T)pic_side-(T)1));
            std::copy( &temp[index*pic_side], &temp[(index+1)*pic_side], &result[i*pic_side]);
        }

        delete[] bi;
        delete[] ri;
        return result.Transpose();
    }
};

} // namespace alias

#endif // ASTROQUT_UTILS_OPERATOR_SPLINE_HPP
