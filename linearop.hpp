///
/// \file include/utils/linearop.hpp
/// \brief LinearOp class header
/// \details Provide generic linear operator base class.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#ifndef ASTROQUT_UTILS_LINEAROP_HPP
#define ASTROQUT_UTILS_LINEAROP_HPP

#include "const.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>

namespace alias
{

/** Types of argument tests
 */
enum ArgTestType {mult, add, element_wise};

class LinearOp
{
protected:
    size_t height_; //!< Member variable "height_"
    size_t width_; //!< Member variable "width_"
    size_t length_; //!< Member variable "length_"

public:

    /** Default constructor
     *  Create an empty container
     */
    LinearOp() noexcept
        : height_(0)
        , width_(0)
        , length_(0)
    {
#ifdef DEBUG
        std::cout << "LinearOp : Default constructor called" << std::endl;
#endif // DEBUG
    }

    /** Copy constructor
     *  \param other Object to copy from
     */
    LinearOp(const LinearOp& other)
        : height_(other.height_)
        , width_(other.width_)
        , length_(other.length_)
    {
#ifdef DEBUG
        std::cout << "LinearOp : Copy constructor called" << std::endl;
#endif // DEBUG
    }

    /** Move constructor
     *  \param other Object to copy from
     */
    LinearOp(LinearOp&& other)
        : LinearOp()
    {
#ifdef DEBUG
        std::cout << "LinearOp : Move constructor called" << std::endl;
#endif // DEBUG
        swap(*this, other);
    }

    /** Full member constructor
     *  \param height Height of the data
     *  \param width Width of the data
     */
    LinearOp(size_t height, size_t width) noexcept
        : height_(height)
        , width_(width)
        , length_(height*width)
    {
#ifdef DEBUG
        std::cout << "LinearOp : Full member constructor called with height=" << height << ", width=" << width << std::endl;
#endif // DEBUG
    }

    /** Default destructor
     */
    virtual ~LinearOp()
    {
#ifdef DEBUG
        std::cout << "LinearOp : Destructor called" << std::endl;
#endif // DEBUG
    }

    /** Access height_
     * \return The current value of height_
     */
    size_t Height() const noexcept
    {
        return height_;
    }
    /** Set height_
     * \param val New value to set
     */
    void Height(size_t val) noexcept
    {
        height_ = val;
        length_ = height_*width_;
    }

    /** Access width_
     * \return The current value of width_
     */
    size_t Width() const noexcept
    {
        return width_;
    }
    /** Set width_
     * \param val New value to set
     */
    void Width(size_t val) noexcept
    {
        width_ = val;
        length_ = height_*width_;
    }

    /** Access length_
     * \return The current value of length_
     */
    size_t Length() const noexcept
    {
        return length_;
    }

    /** Valid instance test
     *  \return Throws an error message if instance is not valid.
     */
    virtual bool IsValid() const = 0;

    /** Argument test
     *  \param other Other object to test
     *  \return Throws an error message if instance is not valid.
     */
    bool ArgTest(const LinearOp& other, ArgTestType type) const
    {
        bool test_result = false;
        switch(type)
        {
        case mult:
        {
            if( IsValid() && other.IsValid() && this->width_ == other.height_ )
                test_result = true;
            break;
        }
        case add:
        {
            if( IsValid() && other.IsValid() && this->length_ == other.length_ )
                test_result = true;
            break;
        }
        case element_wise:
        {
            if( IsValid() && other.IsValid() && this->length_ == other.length_ )
                test_result = true;
            break;
        }
        default:
        {
            break;
        }
        }
        if( test_result )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument("Matrix dimensions must agree, be non-zero and data shall not be empty!");
        }
    }

    /** Swap function
     *  \param first First object to swap
     *  \param second Second object to swap
     */
    friend void swap(LinearOp& first, LinearOp& second) noexcept
    {
        using std::swap;

        swap(first.height_, second.height_);
        swap(first.width_, second.width_);
        swap(first.length_, second.length_);
    }

};
} // namespace alias

#endif // ASTROQUT_UTILS_LINEAROP_HPP
