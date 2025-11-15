//
//  DSTpoint.h
//  DeepSkyTools
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#ifndef __DST_MATH_POINT__
#define __DST_MATH_POINT__

#include <vector>
#include <limits>
#include <cmath>
#include <cstdarg>
#include <stdexcept>
#include <string>

namespace DST
{
    namespace Math
    {
#pragma mark - point class definition
        /**
         *  @class point DSTpoint.h "DSTmath/DSTgeo.h"
         *  @brief Definition of cartesian coordinate of flating point.
         *  @details point class describe cartesian coordinate of floating point. It can handle multidimentional reference frame but it is crrently optimized for 1D, 2D and 3D geometry.
         * @author William GILLARD
         */

        /*!
        * @class point DSTpoint.h "DSTmath/DSTgeo.h"
        * @brief Definition floating point cartesian coordinates.
        */
        class point
        {
#pragma mark -- protected memeber
        protected:
            std::vector<double>fx;
            
        public:
#pragma mark -- ctor/dto r
            point();
            point(const size_t& n):point()
            {
                while(fx.size() < n)
                    fx.push_back(0);
            }
            
            point(const std::initializer_list<double>& );
            point(const std::vector<double>& );

            point(const point&);
            virtual ~point();
            
#pragma mark -- Modifier
            void SetPoints( const std::initializer_list<double>& );
            void SetPoints( const std::vector<double>& );
            void SetPoints( const std::initializer_list<float>& );
            void SetPoints( const std::vector<float>& );
            void SetPoints( const point& _p );

            inline void SetPoint(const size_t& i ,const double& x){if( i < fx.size()) {*(fx.begin()+i) *=0; *(fx.begin()+i) += x;} else throw std::out_of_range("[DST::Math::points::SetPoints(const size_t&, const double&)] out of range.");}
            
            inline void SetX(double x)    {SetPoint(0,x);}
            inline void SetY(double y)    {SetPoint(1,y);}
            inline void SetZ(double z)    {SetPoint(2,z);}
            
#pragma mark -- Operator
            inline const double operator[](const size_t& i) const {return (i < fx.size()) ? *(fx.cbegin()+i) : throw std::out_of_range("[DST::Math::points[]] out of range.");}
            
            bool operator!=(const DST::Math::point&) const;
            bool operator==(const DST::Math::point&) const;
            bool operator<=(const DST::Math::point&) const;
            bool operator< (const DST::Math::point&) const;
            bool operator>=(const DST::Math::point&) const;
            bool operator> (const DST::Math::point&) const;
            
            point& operator=(const point&);
            
            void operator+=(const DST::Math::point&);
            void operator-=(const DST::Math::point&);
            void operator*=(const DST::Math::point&);
            void operator/=(const DST::Math::point&);
            
            void operator+=(const double&);
            void operator-=(const double&);
            void operator*=(const double&);
            void operator/=(const double&);
            
            void operator+=(const int&);
            void operator-=(const int&);
            void operator*=(const int&);
            void operator/=(const int&);
            
#pragma mark -- Accessor
            inline const std::vector<double>& GetPoint() const {return fx;}
            inline size_t size() const {return fx.size();}
            inline double X() const {return operator[](0);}
            inline double Y() const {return operator[](1);}
            inline double Z() const {return operator[](2);}
            double R()     const;
            double Theta() const;
            double Phi()   const;
            
            static double precision;
            
#pragma mark -- Debug
            virtual std::string Dump() const;
            static bool debug;

            template<typename T,
            std::enable_if_t<(std::is_arithmetic<T>::value && !std::is_same<T,bool>::value), int> = 0>
            static point pointFromVector(const std::vector<T>& lst)
            {
                std::vector<double> il;
                for(const auto& val : lst)
                    il.push_back(static_cast<double>(val));
                return point(il);
            }
        };
        
#pragma mark - vector2D class definition

        /*!
        * @class vector2D DSTpoint.h "DSTmath/DSTpoint.h"
        * @brief Definition of 2D vector in spherical coordinate.
        * @details vector2D class describe a 2D vector in spherical coordinate. It is defined by its length and its azimuth angle (phi).
        * @author William GILLARD
        */
        class vector2D
        {
#pragma mark -- protected memeber
        protected:
            double flength;
            double fphi;

            void _phi();
         
         public:
#pragma mark -- ctor/dtor
            vector2D();
            vector2D(const point&);
            vector2D(const point&, const point&);
            vector2D(double, double);
            vector2D(double);
            vector2D(const vector2D&);
         
            virtual ~vector2D();
         
#pragma mark -- Modifier
            inline void SetPhi   (double _t){fphi    = _t; _phi();}
            inline void SetLength(double _l){flength = fabs(_l);}
         
#pragma mark -- Accessor
            inline double Phi      () const {return fphi;}
            inline double Length   () const {return flength;}
            inline virtual double Theta() const {return acos(-1)/2.;}
            
            inline virtual double X() const {return flength*cos(fphi);}
            inline virtual double Y() const {return flength*sin(fphi);}
            inline virtual double Z() const {return 0;}
            
            virtual DST::Math::point EndPoint() const {return point({X(), Y()});}
         
#pragma mark -- Opperator
            virtual bool operator!=(const vector2D&) const;
            virtual bool operator==(const vector2D&) const;
            virtual bool operator<=(const vector2D&) const;
            virtual bool operator< (const vector2D&) const;
            virtual bool operator>=(const vector2D&) const;
            virtual bool operator> (const vector2D&) const;
            
            vector2D operator=(const vector2D&);

            void operator+=(const vector2D&);
            void operator-=(const vector2D&);
            void operator*=(const double);
            void operator/=(const double);
            
#pragma mark -- Dump
            virtual std::string Dump() const;
            
         };
        
#pragma mark - vector3D class definition

        /*!
        * @brief Definition of 3D vector in spherical coordinate.
        * @details vector3D class describe a 3D vector in spherical coordinate. It is defined by its length, its azimuth angle (phi) and its polar angle (theta).
        * @author William GILLARD
        */
        class vector3D: public vector2D
        {
#pragma mark -- protected memeber
        protected:
            double ftheta;

            void _theta();
            
        public:
#pragma mark -- ctor/dtor
            vector3D();
            vector3D(const point&);
            vector3D(const point&,const point&);
            vector3D(double, double, double);
            vector3D(double, double);
            vector3D(const vector3D&);
            vector3D(const vector2D&);
            
            virtual ~vector3D();
            
#pragma mark -- Modifier
            inline void SetTheta  (double _p){ftheta   = _p;}
            
#pragma mark -- Accessor
            inline double Theta() const {return ftheta;}
            
            inline virtual double X() const {return flength*sin(ftheta)*cos(fphi);}
            inline virtual double Y() const {return flength*sin(ftheta)*sin(fphi);}
            inline virtual double Z() const {return flength*cos(ftheta);}
            
            virtual DST::Math::point EndPoint() const {return point({X(), Y(), Z()});}
            
#pragma mark -- Opperator
#ifdef Linuxx86_64
            using vector2D::operator<=;
            using vector2D::operator<;
            using vector2D::operator>=;
            using vector2D::operator>;
            using vector2D::operator==;
            using vector2D::operator!=;
#endif

            bool operator!=(const vector3D&) const;
            bool operator==(const vector3D&) const;
            bool operator<=(const vector3D&) const;
            bool operator< (const vector3D&) const;
            bool operator>=(const vector3D&) const;
            bool operator> (const vector3D&) const;
            
            vector3D operator=(const vector3D&);
            
            void operator^=(const vector3D&);
            void operator+=(const vector3D&);
            void operator-=(const vector3D&);
            void operator^=(const vector2D&);
            void operator+=(const vector2D&);
            void operator-=(const vector2D&);
            void operator+=(const double&);
            void operator-=(const double&);
            void operator/=(const double&);
            void operator*=(const double&);
            
#pragma mark -- Dump
            virtual std::string Dump() const;
        };

    }
}

#pragma mark operator

DST::Math::point operator+(const DST::Math::point& p1, const DST::Math::point& p2);
DST::Math::point operator-(const DST::Math::point& p1, const DST::Math::point& p2);
DST::Math::point operator*(const DST::Math::point& p1, const DST::Math::point& p2);
DST::Math::point operator/(const DST::Math::point& p1, const DST::Math::point& p2);

DST::Math::point operator+(const DST::Math::point& p1, const double s);
DST::Math::point operator-(const DST::Math::point& p1, const double s);
DST::Math::point operator*(const DST::Math::point& p1, const double s);
DST::Math::point operator/(const DST::Math::point& p1, const double s);

DST::Math::point operator+(const double s, const DST::Math::point& p1);
DST::Math::point operator-(const double s, const DST::Math::point& p1);
DST::Math::point operator*(const double s, const DST::Math::point& p1);
DST::Math::point operator/(const double s, const DST::Math::point& p1);

double operator*(const DST::Math::vector2D& v1, const DST::Math::vector2D& v2);
double operator*(const DST::Math::vector2D& v1, const DST::Math::vector3D& v2);
double operator*(const DST::Math::vector3D& v1, const DST::Math::vector2D& v2);
double operator*(const DST::Math::vector3D& v1, const DST::Math::vector3D& v2);

DST::Math::vector2D operator*(const DST::Math::vector2D& v1, const double s);
DST::Math::vector2D operator*(const double s, const DST::Math::vector2D& v1);

DST::Math::vector2D operator/(const DST::Math::vector2D& v1, const double s);

DST::Math::vector3D operator*(const DST::Math::vector3D& v1, const double s);
DST::Math::vector3D operator*(const double s, const DST::Math::vector3D& v1);

DST::Math::vector3D operator/(const DST::Math::vector3D& v1, const double s);

DST::Math::vector3D operator^(const DST::Math::vector2D& v1, const DST::Math::vector2D& v2);
DST::Math::vector3D operator^(const DST::Math::vector3D& v1, const DST::Math::vector2D& v2);
DST::Math::vector3D operator^(const DST::Math::vector2D& v1, const DST::Math::vector3D& v2);
DST::Math::vector3D operator^(const DST::Math::vector3D& v1, const DST::Math::vector3D& v2);

#endif
