///
/// \file include/utils/linearop/operator/matmult.hpp
/// \brief Matrix Multiplication class header
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_OPERATOR_MATMULT_HPP
#define ASTROQUT_UTILS_OPERATOR_MATMULT_HPP

#include "utils/linearop/operator.hpp"

namespace alias
{

template <class T = double>
class MatMult : public Operator<T>
{
public:

    /** Default constructor
     */
    MatMult()
        : Operator<T>()
    {
#ifdef DEBUG
        std::cout << "MatMult : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Copy constructor
     *  \param other Object to copy from
     */
    MatMult(const MatMult& other)
        : Operator<T>(other)
    {
#ifdef DEBUG
        std::cout << "MatMult : Copy constructor called" << std::endl;
#endif // DEBUG
    }

    /** Move constructor
     *  \param other Object to move from
     */
    MatMult(MatMult&& other)
        : MatMult()
    {
#ifdef DEBUG
        std::cout << "MatMult : Move constructor called" << std::endl;
#endif // DEBUG
        swap(*this, other);
    }

    /** Full member constructor
     *  \param data Array containing data used by the operator
     *  \param height Height of the operator
     *  \param width Width of the operator
     */
    explicit MatMult(Matrix<T> data, size_t height, size_t width)
        : Operator<T>(data, height, width, false)
    {
#ifdef DEBUG
        std::cout << "MatMult : Full member constructor called" << std::endl;
#endif // DEBUG
    }

    /** Clone function
     *  \return A copy of the current instance
     */
    virtual MatMult* Clone() const override
    {
        return new MatMult(*this);
    }

    /** Default destructor
     */
    virtual ~MatMult()
    {
#ifdef DEBUG
        std::cout << "MatMult : Destructor called" << std::endl;
#endif // DEBUG
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    bool IsValid() const override final
    {
        if( this->height_ != 0 && this->width_ != 0 && !this->data_.IsEmpty() )
            return true;
        else
            throw std::invalid_argument("Operator dimensions must be non-zero and function shall not be nullptr!");
    }

    /** Swap function
     *  \param first First object to swap
     *  \param second Second object to swap
     */
    friend void swap(MatMult& first, MatMult& second) noexcept
    {
        using std::swap;

        swap(static_cast<Operator<T>&>(first), static_cast<Operator<T>&>(second));
    }

    /** Copy assignment operator
     *  \param other Object to assign to current object
     *  \return A reference to this
     */
    MatMult& operator=(MatMult other)
    {
        swap(*this, other);

        return *this;
    }

    Matrix<T> operator*(const Matrix<T>& other) const override final
    {
#ifdef DEBUG
        std::cerr << "MatMult: operator* called" << std::endl;
#endif // DEBUG
#ifdef DO_ARGCHECKS
        try
        {
            this->ArgTest(other, mult);
        }
        catch (const std::exception&)
        {
            throw;
        }
#endif // DO_ARGCHECKS

        return std::move(this->data_ * other);
    }

    /** Transpose in-place
     *   \return A reference to this
     */
    virtual MatMult& Transpose() override
    {
        std::swap(this->height_, this->width_);
        this->transposed_ = !this->transposed_;
        std::move(this->data_).Transpose();
        return *this;
    }

};

} // namespace alias

#endif // ASTROQUT_UTILS_OPERATOR_MATMULT_HPP
