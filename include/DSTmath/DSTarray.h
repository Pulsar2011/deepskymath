//
//  DSTarray.h
//  DeepSkyTools
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#ifndef DSM_ARRAY_H
#define DSM_ARRAY_H

#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace DST
{
    namespace Math
    {
#pragma region - Masked array class definition
        /**
         *  @class DST::Math::MaskedArray DSTarray.h "DSTmath/DSTarray.h"
         *  @brief Definition of a masked array.
         *  @details MaskedArray class describe a masked array. It is a wrapper around std::valarray to handle masked array.
         */
        template <class T>
        class MaskedArray
        {
#pragma region -- private members
        protected:
            std::valarray<T> fData;
            std::valarray<bool> fMask;

#pragma endregion
#pragma region -- public members functions
        public:
#pragma region ctor/dtor
            /**
             *  @brief Default constructor.
             *  @details Create an empty masked array.
             */
            MaskedArray():fData(),fMask() {};

            /**
             *  @brief Constructor with size.
             *  @details Create a masked array with a given size.
             *  @param size Size of the masked array.
             */
            MaskedArray(size_t size):fData(size),fMask(size)
            {
                fData *= 0;
                fMask &= false;
            };

            /**
             *  @brief Constructor with size and value.
             *  @details Create a masked array with a given size and a default value.
             *  @param size Size of the masked array.
             *  @param value Default value of the masked array.
             */
            MaskedArray(size_t size, T value):fData(size),fMask(size)
            {
                fData += value;
                fMask &= false;
            };

            /**
             *  @brief Constructor with array and mask.
             *  @details Create a masked array with a given size and a default value.
             *  @param value Default value of the masked array.
             *  @param mask Mask of the masked array.
             */
            MaskedArray(std::valarray<T> value, std::valarray<bool> mask):fData(value.size()),fMask(mask.size())
            {
                if (value.size() != mask.size())
                    throw std::invalid_argument("MaskedArray::MaskedArray: value and mask size mismatch ["+std::to_string(__LINE__)+"]");
                
                fData += value;
                fMask |= mask;
            };

            /**
             *  @brief Constructor with array.
             *  @details Create a masked array with a given size and a default value.
             *  @param value Default value of the masked array.
             */
            MaskedArray(std::valarray<T> value):fData(value.size()),fMask(value.size())
            {
                fData += value;
                fMask &= false;
            };

            /**
             *  @brief Constructor with array and mask.
             *  @details Create a masked array with a given size and a default value.
             *  @param value Default value of the masked array.
             *  @param mask Mask of the masked array.
             */
            MaskedArray(std::vector<T> value, std::vector<bool> mask):fData(value.size()),fMask(mask.size())
            {
                if (value.size() != mask.size())
                    throw std::invalid_argument("MaskedArray::MaskedArray: value and mask size mismatch ["+std::to_string(__LINE__)+"]");
                
                std::copy(value.begin(), value.end(), std::begin(fData)); 
                std::copy(mask.begin(), mask.end(), std::begin(fMask)); 
            };

            /**
             *  @brief Constructor with array.
             *  @details Create a masked array with a given size and a default value.
             *  @param value Default value of the masked array.
             *  @param mask Mask of the masked array.
             */
            MaskedArray(std::vector<T> value):fData(value.size()),fMask(value.size())
            {
                std::copy(value.begin(), value.end(), std::begin(fData)); 
                fMask &= false;
            };

            /**
             *  @brief Copy constructor.
             *  @details Create a masked array with a given size and a default value.
             *  @param value Default value of the masked array.
             *  @param mask Mask of the masked array.
             */
            MaskedArray(const MaskedArray<T>& ma):fData(ma.fData.size()),fMask(ma.fMask.size())
            {
                fData += ma.fData;
                fMask |= ma.fMask;
            }

            /**
             *  @brief Destructor.
             *  @details Destroy the masked array.
             */
            ~MaskedArray()
            {
                fData.resize(0);
                fMask.resize(0);
            };

            /**
             * @brief Assignement operator
             * 
             */
            MaskedArray<T>& operator=(const MaskedArray<T>& ma)
            {
                fData.resize(ma.fData.size());
                fMask.resize(ma.fMask.size());

                fData = ma.fData;
                fMask = ma.fMask;
                return *this;
            }

            /**
             * @brief Assignement operator
             * 
             */
            MaskedArray<T>& operator=(const std::slice_array<T>& ma)
            {
                std::valarray<T> _data=ma;
                fData.resize(_data.size());
                fMask.resize(fData.size());
                
                fData *= 0;
                fData += _data;
                fMask &= false;
                
                return *this;
            }

            /**
             * @brief Assignement operator
             * 
             */
            MaskedArray<T>& operator=(const std::gslice_array<T>& ma)
            {
                std::valarray<T> _data=ma;
                fData.resize(_data.size());
                fMask.resize(fData.size());

                fData *= 0;
                fData += _data;
                fMask &= false;
                return *this;
            }

            /**
             * @brief Assignement operator
             * 
             */
            MaskedArray<T>& operator=(const std::mask_array<T>& ma)
            {
                std::valarray<T> _data=ma;
                fData.resize(_data.size());
                fMask.resize(fData.size());

                fData *= 0;
                fData += _data;
                fMask &= false;
                return *this;
            }

            /**
             * @brief Assignement operator
             * 
             */
            MaskedArray<T>& operator=(const std::indirect_array<T>& ma)
            {
                std::valarray<T> _data=ma;
                fData.resize(_data.size());
                fMask.resize(fData.size());

                fData *= 0;
                fData += _data;
                fMask &= false;
                return *this;
            }

            /**
             * @brief Assignement operator
             * 
             */
            MaskedArray<T>& operator=(std::initializer_list<T> il)
            {
                std::valarray<T> _data=il;
                fData.resize(_data.size());
                fMask.resize(fData.size());

                fData *= 0;
                fData += _data;
                fMask &= false;
                return *this;
            }

#pragma endregion
#pragma region accessors
            /**
             *  @brief Accessor to the data.
             *  @details Return reference to the data array.
             *  @return  reference to the fData array.
             */
            std::valarray<T>& data() { return fData; }
            const std::valarray<T>& data() const { return fData; }

            /**
             *  @brief Accessor to the mask.
             *  @details Return reference to the mask array.
             *  @return  reference to the fMask array.
             */
            std::valarray<bool>& mask() { return fMask; }
            const std::valarray<bool>& mask() const { return fMask; }

            /**
             *  @brief Accessor to the size.
             *  @details Return the size of the masked array.
             *  @return The size of the masked array.
             */
            size_t size() const {return fData.size();}

            /**
             * @brief Count number of unmasked data
             * @return Number of unmasked value
             */
            size_t count() const
            {
                return fMask[!fMask].size();
            }

            /**
             * @brief Access data at index index
             * @param idx Index of the data to be accessed
             * @return constante reference at position idx
             */
            const T& operator[](size_t idx) const
            {
                if( idx >= fData.size() )
                    throw std::out_of_range("MaskedArray::operator[]: index out of range ["+std::to_string(__LINE__)+"]");

                return fData[idx];
            }

            /**
             * @brief Access data at index index
             * @param idx Index of the data to be accessed
             * @return  reference at position idx
             */
            T& operator[](size_t idx)
            {
                if( idx >= fData.size() )
                    throw std::out_of_range("MaskedArray::operator[]: index out of range ["+std::to_string(__LINE__)+"]");

                return fData[idx];
            }

#pragma endregion
#pragma region modifiers

            /**
             * @brief Reset the data object preserving the original size of the array
             * @param data value to be set to the data array 
             */
            void set_data(const std::valarray<T>& data)
            {
                fData.resize(data.size());
                fMask.resize(fData.size());
                
                fData *= 0;
                fData += data;
            }

            /**
             * @brief Change an element of the data array
             * @param data value to be set to the data array 
             */
            void set_data(const std::valarray<T>& data, const std::valarray<bool>& mask)
            {
                if( mask.size() != data.size() )
                    throw std::out_of_range("MaskedArray::set_data: unconsistant data size ["+std::to_string(__LINE__)+"]");

                fData.resize(data.size());
                fMask.resize(fData.size());
                
                fData *= 0;
                fData += data;
                
                fMask &= false;
                fMask |= mask;
            }

            /**
             * @brief Reset the data mask preserving the original size of the array
             * @param mask to be set to the data array 
             */
            void set_mask(const std::valarray<bool>& mask)
            {
                if (mask.size() != fMask.size())
                    throw std::invalid_argument("MaskedArray::set_mask: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask &= false;
                fMask |= mask;
            }

            /**
             * @brief Reset the data mask preserving the original size of the array
             * @param mask to be set to the data array 
             */
            void set_mask(const size_t& idx, const bool& mask)
            {
                if (idx >= fMask.size())
                    throw std::invalid_argument("MaskedArray::set_mask: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask[idx] = mask;
            }

            /**
             * @brief swap content of two MaskedArray
             * @param ma MaskedArray to be swapped whith this 
             */
            void swap(MaskedArray<T>& ma)
            {
                std::swap(fData,ma.fData);
                std::swap(fMask,ma.fMask);
            }

            /**
             * @brief resize the data and mask array
             * @param size new size of the data and mask array
             */
            void resize(size_t size)
            {
                fData.resize(size);
                fMask.resize(size);
            }

            /**
             * @brief evaluate the sum of unmasked data
             * @return the arithmetic sum of unmasked data
             */
            T sum() const
            {
                std::valarray<T> _v = fData[!fMask];
                return _v.sum();
            }

            /**
             * @brief evaluate the minimum of unmasked data
             * @return the minimum of unmasked data
             */
            T min() const
            {
                std::valarray<T> _v = fData[!fMask];
                return _v.min();
            }

            /**
             * @brief evaluate the maximum of unmasked data
             * @return the maximum of unmasked data
             */
            T max() const
            {
                std::valarray<T> _v = fData[!fMask];
                return _v.max();
            }

            /**
             * @brief evaluate the mean of unmasked data
             * @return the arithmetic mean of unmasked data
             */
            double mean() const
            {
                std::valarray<T> _v = fData[!fMask];
                return static_cast<double>(_v.sum())/static_cast<double>(_v.size());
            }

            /**
             * @brief evaluate the sum of square of unmasked data
             * @return the arithmetic sum of square of unmasked data
             */
            T sum_of_square() const
            {
                std::valarray<T> _v = fData[!fMask];
                return (_v*_v).sum();
            }

            /**
             * @brief evaluate the mean of square of unmasked data
             * @return the arithmetic mean of square of unmasked data
             */
            double mean_of_square() const
            {
                std::valarray<T> _v = fData[!fMask];
                return static_cast<double>(sum_of_square())/static_cast<double>(_v.size());
            }

            /**
             * @brief evaluate the quadratic mean of unmasked data
             * @return the quadratic mean of unmasked data
             */
            double quadratic_mean() const
            {
                return sqrt(mean_of_square());
            }

            /**
             * @brief evaluate the variance of unmasked data
             * @return the variance of unmasked data
             */
            double variance() const
            {
                double _m = mean();
                double _n = static_cast<double>(count());
                return ( mean_of_square() - _m*_m )* (_n) / (_n - 1) ;
            }

            /**
             * @brief evaluate the error on the variance of unmasked data
             * @return the error on the variance of unmasked data
             */
            double variance_error() const
            {
                double _v = 0;
                double _n = 0;
                double vv = variance();
                double m = mean();
                
                for (size_t i = 0; i < fData.size(); i++)
                {
                    if (!fMask[i])
                    {
                        _v += (static_cast<double>(fData[i]) - m)*(static_cast<double>(fData[i]) - m)*(static_cast<double>(fData[i]) - m)*(static_cast<double>(fData[i]) - m);
                        _n++;
                    }
                }

                _v /= _n;
                double b = (_n - 3) / (_n - 1) * vv * vv;
                
                return sqrt((_v - b) / _n);
            }

            /**
             * @brief evaluate the covariance of unmasked data between a set of two MaskedArray 
             * @param ma MaskedArray to be used to compute the covariance
             * @return the covariance between the two MaskedArray
             */
            double covariance(const MaskedArray<T>& ma) const
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::covariance: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                std::valarray<bool> msk_array(fMask.size());
                msk_array &= false;
                msk_array |= fMask;
                msk_array |= ma.fMask;

                std::valarray<T> _v1 = fData[!msk_array];
                std::valarray<T> _v2 = ma.fData[!msk_array];

                return static_cast<double>((_v2*_v1).sum())/static_cast<double>(_v1.size());
            }

            /**
             * @brief Identify the first index of the minimum value of unmasked data
             * @return Index of the minimum value which is unmasked
             */
            size_t find_minimum() const
            {   
                std::valarray<T> readings(fData.size());
                readings *= 0;
                readings += fData;

                readings[fMask] = std::numeric_limits<T>::max();

                size_t index = std::distance(std::begin(readings), std::min_element(std::begin(readings), std::end(readings)));

                return index;
            }

            /**
             * @brief Identify the first index of the maximum value of unmasked data
             * @return Index of the maximum value which is unmasked
             */
            size_t find_maximum() const
            {   
                std::valarray<T> readings(fData.size());
                readings *= 0;
                readings += fData;

                readings[fMask] = -1*std::numeric_limits<T>::max();

                size_t index = std::distance(std::begin(readings), std::max_element(std::begin(readings), std::end(readings)));
                
                return index;
            }

            /**
             * @brief evaluate the median of unmasked data
             * @return the median of unmasked data
             */
            T median() const
            {
                std::valarray<T> readings(count());
                readings *= 0;
                readings += (std::valarray<T>) fData[!fMask];

                if (readings.size() < 1)
                    throw std::invalid_argument("MaskedArray::median: no umasked data available do compute median ["+std::to_string(__LINE__)+"]");


                std::sort(std::begin(readings), std::end(readings));

                size_t index = 0;
                if (readings.size() < 2)
                    return readings[index];
                else if (!(readings.size() % 2))
                {
                    index = readings.size()/2;
                    return (readings[index-1]+readings[index])/2;
                }
                else
                {
                    index = readings.size()/2;
                    return readings[index];
                }
            
                return readings[index];
            }

            /**
             * @brief evaluate the median absolute deviation of unmasked data
             * @return the median absolute deviation of unmasked data
             */
            T nmad() const
            {
                T med = median();
                MaskedArray<T> mad_array(std::abs(fData[!fMask]-med));
                return mad_array.median();
            }



#pragma endregion
#pragma region operator

            /**
             * @brief operator to add two MaskedArray
             * @param ma MaskedArray to be added to this
             * @details 'this' is transformed with the addition of ma. The addition of the masked data is evaluated as if unmasked but the fMask of this is updated to mask the masked data of ma
             */
            void operator+=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator+=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fData += ma.fData;
                fMask |= ma.fMask;
            }

            /**
             * @brief operator to substract two MaskedArray
             * @param ma MaskedArray to be substracted to this
             * @details 'this' is transformed with the substraction of ma. The substraction of the masked data is evaluated as if unmasked but the fMask of this is updated to mask the masked data of ma
             */
            void operator-=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator-=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fData -= ma.fData;
                fMask |= ma.fMask;
            }

            /**
             * @brief operator to multiply two MaskedArray
             * @param ma MaskedArray to be multiplied to this
             * @details 'this' is transformed with the multiplication of ma. The multiplication of the masked data is evaluated as if unmasked but the fMask of this is updated to mask the masked data of ma
             */
            void operator*=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator*=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fData *= ma.fData;
                fMask |= ma.fMask;
            }

            /**
             * @brief operator to divide two MaskedArray
             * @param ma MaskedArray to be divided to this
             * @details 'this' is transformed with the division of ma. The division of the masked data is evaluated as if unmasked but the fMask of this is updated to mask the masked data of ma
             */
            void operator/=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator/=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fData /= ma.fData;
                fMask |= ma.fMask;
            }

            /**
             * @brief operator to add a value to the MaskedArray
             * @param val value to be added to the MaskedArray
             * @details 'this' is transformed with the addition of val. The addition of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator+=(const T& val){fData += val;}

            /**
             * @brief operator to substract a value to the MaskedArray
             * @param val value to be substracted to the MaskedArray
             * @details 'this' is transformed with the subsraction of val. The substraction of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator-=(const T& val){fData -= val;}

            /**
             * @brief operator to multiply a value to the MaskedArray
             * @param val value to be multiplied to the MaskedArray
             * @details 'this' is transformed with the multiplication of val. The multiplication of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator*=(const T& val){fData *= val;}

            /**
             * @brief operator to divide a value to the MaskedArray
             * @param val value to be divided to the MaskedArray
             * @details 'this' is transformed with the division of val. The division of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator/=(const T& val){fData /= val;}

            /**
             * @brief operator to add values to the MaskedArray
             * @param val value to be added to the MaskedArray
             * @details 'this' is transformed with the addition of val. The addition of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator+=(const std::valarray<T>& val)
            {
                if (val.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator+=: data size mismatch ["+std::to_string(__LINE__)+"]");
                fData += val;
            }

            /**
             * @brief operator to substract values to the MaskedArray
             * @param val value to be substracted to the MaskedArray
             * @details 'this' is transformed with the substraction of val. The substraction of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator-=(const std::valarray<T>& val)
            {
                if (val.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator-=: data size mismatch ["+std::to_string(__LINE__)+"]");
                fData -= val;
            }

            /**
             * @brief operator to multiply values to the MaskedArray
             * @param val value to be multiplied to the MaskedArray
             * @details 'this' is transformed with the multiplication of val. The multiplication of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator*=(const std::valarray<T>& val)
            {
                if (val.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator*=: data size mismatch ["+std::to_string(__LINE__)+"]");
                fData *= val;
            }

            /**
             * @brief operator to divide values to the MaskedArray
             * @param val value to be divided to the MaskedArray
             * @details 'this' is transformed with the division of val. The division of the masked data is evaluated as if unmasked but mask isn't modified as no mask information is provided as input
             */
            void operator/=(const std::valarray<T>& val){fData /= val;}

            void operator|=(const bool& mask){fMask |= mask;}
            void operator&=(const bool& mask){fMask &= mask;}
            void operator^=(const bool& mask){fMask ^= mask;}

            void operator|=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator|=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask |= ma.fMask;
            }

            void operator&=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator&=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask &= ma.fMask;
            }

            void operator^=(const MaskedArray<T>& ma)
            {
                if (ma.size() != fData.size())
                    throw std::invalid_argument("MaskedArray::operator^=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask ^= ma.fMask;
            }

            void operator|=(const std::valarray<bool>& mask)
            {
                if (mask.size() != fMask.size())
                    throw std::invalid_argument("MaskedArray::operator|=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask |= mask;
            }

            void operator&=(const std::valarray<bool>& mask)
            {
                if (mask.size() != fMask.size())
                    throw std::invalid_argument("MaskedArray::operator&=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask &= mask;
            }

            void operator^=(const std::valarray<bool>& mask)
            {
                if (mask.size() != fMask.size())
                    throw std::invalid_argument("MaskedArray::operator^=: data size mismatch ["+std::to_string(__LINE__)+"]");
                
                fMask ^= mask;
            }
        };

#pragma endregion
#pragma endregion
#pragma endregion
#pragma region - Logical Opperator
        template <class T>
        std::valarray<bool> operator!=(const MaskedArray<T>& a, const MaskedArray<T>& b)
        {
            if (a.size() != b.size())
                throw std::invalid_argument("MaskedArray::operator!=: data size mismatch ["+std::to_string(__LINE__)+"]");

            std::valarray<bool> _ma(a.size());
            _ma = (a.data() != b.data()) && (a.mask() != b.mask());
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator!=(const MaskedArray<T>& a, const T& b)
        {
            std::valarray<bool> _ma(a.size());
            _ma = (a.data() != b) && !a.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator!=(const T& b,const MaskedArray<T>& a)
        {
            return a != b;
        }

        template <class T>
        std::valarray<bool> operator==(const MaskedArray<T>& a, const MaskedArray<T>& b)
        {
            if (a.size() != b.size())
                throw std::invalid_argument("MaskedArray::operator==: data size mismatch ["+std::to_string(__LINE__)+"]");

            std::valarray<bool> _ma(a.size());
            _ma = (a.data() == b.data()) && (a.mask() == b.mask());
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator==(const MaskedArray<T>& a, const T& b)
        {
            std::valarray<bool> _ma(a.size());
            _ma = (a.data() == b)&&(!a.mask());
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator==(const T& b, const MaskedArray<T>& a)
        {
            return a == b;
        }

        template <class T>
        std::valarray<bool> operator>(const MaskedArray<T>& a, const MaskedArray<T>& b)
        {
            if (a.size() != b.size())
                throw std::invalid_argument("MaskedArray::operator>: data size mismatch ["+std::to_string(__LINE__)+"]");

            std::valarray<bool> _ma(a.size());
            _ma = (a.data() > b.data()) && !a.mask() && !b.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator>(const MaskedArray<T>& a, const T& b)
        {
            std::valarray<bool> _ma(a.size());
            _ma = (a.data() > b) && !a.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator>(const T& b, const MaskedArray<T>& a)
        {
            return a < b;
        }

        template <class T>
        std::valarray<bool> operator<(const MaskedArray<T>& a, const MaskedArray<T>& b)
        {
            if (a.size() != b.size())
                throw std::invalid_argument("MaskedArray::operator<: data size mismatch ["+std::to_string(__LINE__)+"]");

            std::valarray<bool> _ma(a.size());
            _ma = (a.data() < b.data()) && !a.mask() && !b.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator<(const MaskedArray<T>& a, const T& b)
        {
            std::valarray<bool> _ma(a.size());
            _ma = (a.data() < b) && !a.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator<(const T& b,const MaskedArray<T>& a)
        {
            return a > b;
        }

        template <class T>
        std::valarray<bool> operator>=(const MaskedArray<T>& a, const MaskedArray<T>& b)
        {
            if (a.size() != b.size())
                throw std::invalid_argument("MaskedArray::operator>=: data size mismatch ["+std::to_string(__LINE__)+"]");

            std::valarray<bool> _ma(a.size());
            _ma = (a.data() >= b.data()) && !a.mask() && !b.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator>=(const MaskedArray<T>& a, const T& b)
        {
            std::valarray<bool> _ma(a.size());
            _ma = (a.data() >= b) && !a.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator>=(const T& b,const MaskedArray<T>& a)
        {
            return a <= b;
        }

        template <class T>
        std::valarray<bool> operator<=(const MaskedArray<T>& a, const MaskedArray<T>& b)
        {
            if (a.size() != b.size())
                throw std::invalid_argument("MaskedArray::operator<=: data size mismatch ["+std::to_string(__LINE__)+"]");

            std::valarray<bool> _ma(a.size());
            _ma = (a.data() <= b.data()) && !a.mask() && !b.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator<=(const MaskedArray<T>& a, const T& b)
        {
            std::valarray<bool> _ma(a.size());
            _ma = (a.data() <= b) && !a.mask();
            return _ma;
        }

        template <class T>
        std::valarray<bool> operator<=(const T& b, const MaskedArray<T>& a)
        {
            return a >= b;
        }
    }
}

#pragma endregion
#pragma region - Helper functions

template <class T>
std::ostream& operator<<(std::ostream& os, const DST::Math::MaskedArray<T>& v)
{
    for (size_t i = 0; i < v.size(); i++)
        os << ' ' << v.data()[i] << '[' << ((v.mask()[i])?"\033[32mtrue\033[0m":"\033[31mfalse\033[0m")<<"]";
    return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& v)
{
    for (size_t i = 0; i < v.size(); i++)
        os << ' ' << v[i];
    return os;
}

#pragma endregion


#endif