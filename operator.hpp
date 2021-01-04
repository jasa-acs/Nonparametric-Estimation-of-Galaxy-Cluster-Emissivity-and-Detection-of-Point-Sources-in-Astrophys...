///
/// \file include/utils/linearop/operator.hpp
/// \brief Operator class header
/// \details Provide operator container with operator-vector operations
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_HPP
#define ASTROQUT_UTILS_OPERATOR_HPP

#include "utils/linearop/matrix.hpp"

#include <functional>
#include <utility>

namespace alias
{

template<class T = double>
class Operator : public LinearOp
{
protected:
    Matrix<T> data_; //!< Member variable "data_"
    bool transposed_; //!< Member variable "transposed_"

public:
    /** Default constructor
     */
    Operator()
        : data_()
        , transposed_(false)
    {
#ifdef DEBUG
        std::cout << "Operator : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Copy constructor
     *  \param other Object to copy from
     */
    Operator(const Operator& other)
        : LinearOp(other.height_, other.width_)
        , data_(other.data_)
        , transposed_(other.transposed_)
    {
#ifdef DEBUG
        std::cout << "Operator : Copy constructor called" << std::endl;
#endif // DEBUG
    }

    /** Move constructor
     *  \param other Object to copy from
     */
    Operator(Operator&& other)
        : LinearOp()
    {
#ifdef DEBUG
        std::cout << "Operator : Move constructor called" << std::endl;
#endif // DEBUG
        swap(*this, other);
    }

    /** Empty constructor
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    Operator(size_t height, size_t width) noexcept
        : LinearOp(height, width)
        , data_()
        , transposed_(false)
    {
#ifdef DEBUG
        std::cout << "Operator : Empty constructor called with height=" << height << ", width=" << width << std::endl;
#endif // DEBUG
    }

    /** Full member constructor
     *  \param data Array containing data used by the operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     *  \param transposed If the operator is transposed
     */
    Operator(Matrix<T> data, size_t height, size_t width, bool transposed) noexcept
        : LinearOp(height, width)
        , data_(Matrix<T>(data))
        , transposed_(transposed)
    {
#ifdef DEBUG
        std::cout << "Operator : Full member constructor called with data=" << &data << ", height=" << height << ", width=" << width << ", transposed=" << transposed << std::endl;
#endif // DEBUG
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    virtual Operator* Clone() const = 0;

    /** Default destructor
     */
    virtual ~Operator()
    {
#ifdef DEBUG
        std::cout << "Operator : Destructor called" << std::endl;
#endif // DEBUG
    }

    /** Access data_
     * \return The current value of data_
     */
    const Matrix<T>& Data() const noexcept
    {
        return data_;
    }
    /** Set data_
     * \param data New value to set
     */
    void Data(const Matrix<T>& data) noexcept
    {
        data_ = data;
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    virtual bool IsValid() const override
    {
        if( this->height_ != 0 && this->width_ != 0 )
            return true;
        else
            throw std::invalid_argument("Operator dimensions must be non-zero!");
    }

    /** Array subscript setter operator
     *  \param index Array subscript
     *  \return A reference to the array element at index
     */
    T& operator[](size_t index) noexcept
    {
        return data_[index];
    }

    /** Array subscript getter operator
     *  \param index Array subscript
     *  \return A reference to the array element at index
     */
    template <class S = T, typename std::enable_if_t<std::is_arithmetic<S>::value>* = nullptr>
    const S operator[](size_t index) const noexcept
    {
        return data_[index];
    }
    template <class S = T, typename std::enable_if_t<is_complex<S> {}>* = nullptr>
    const S& operator[](size_t index) const noexcept
    {
        return data_[index];
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    virtual Operator& Transpose()
    {
        std::swap(this->height_, this->width_);
        transposed_ = !transposed_;
        return *this;
    }

    /** Access transposed_
     * \return The current value of transposed_
     */
    bool IsTransposed() const noexcept
    {
        return transposed_;
    }

    /** Comparison operator equal
     *  \param other Object to compare to
     *  \return True if both object are the same, False otherwise
     */
    bool operator==(const Operator& other) const noexcept
    {
        if( this->height_ != other.height_ || this->width_ != other.width_ || data_ != other.data_ )
            return false;

        return true;
    }

    /** Comparison operator not-equal
     *  \param other Object to compare to
     *  \return False if both object are the same element-wise, True else
     */
    bool operator!=(const Operator& other) const noexcept
    {
        return !(*this == other);
    }

    /** Swap function
     *  \param first First object to swap
     *  \param second Second object to swap
     */
    friend void swap(Operator& first, Operator& second) noexcept
    {
        using std::swap;

        swap(static_cast<LinearOp&>(first), static_cast<LinearOp&>(second));
        swap(first.data_, second.data_);
        swap(first.transposed_, second.transposed_);
    }

    virtual Matrix<T> operator*(const Matrix<T>& other) const = 0;
};

} // namespace alias

#endif // ASTROQUT_UTILS_OPERATOR_HPP
