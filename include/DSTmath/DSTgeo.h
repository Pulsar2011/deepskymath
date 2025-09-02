//
//  DSTpoint.h
//  DeepSkyTools
//
//  File created by GILLARD William on 04/11/15.
//  Definition of cartesian coordinates
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#ifndef __DST_MATH_POINT__
#define __DST_MATH_POINT__

#include <vector>
#include <limits>
#include <cmath>
#include <cstdarg>

std::string to_string(const double, const int);
std::string to_string(const float, const int);

namespace DST
{
    namespace Math
    {
#pragma mark - point class definition
        /**
         *  @class point DSTpoint.h "DSTmath/DSTpoint.h"
         *  @brief Definition of cartesian coordinate of flating point.
         *  @details point class describe cartesian coordinate of floating point. It can handle multidimentional reference frame but it is crrently optimized for 1D, 2D and 3D geometry.
         */
        class point
        {
#pragma mark • protected memeber
        protected:
            std::vector<double>fx;
            
        public:
#pragma mark • ctor/dto r
            point();
            point(const size_t& n)
            {
                fx = std::vector<double>();
                while(fx.size() < n)
                    fx.push_back(0);
            }
            
            point(const size_t& n, double, ...);
            point(const point&);
            virtual ~point();
            
#pragma mark • Modifier
            void SetPoints( const std::initializer_list<double>& );
            void SetPoints( const std::vector<double>& );
            void SetPoints( const std::initializer_list<float>& );
            void SetPoints( const std::vector<float>& );
            void SetPoints( const point& _p );


            inline void SetPoint(const size_t& i ,const double& x){if( i < fx.size()) {fx[i]*=0; fx[i] += x;} else throw std::out_of_range("[DST::Math::points::SetPoints(const size_t&, const double&)] out of range.");}
            
            inline void SetX(double x)    {SetPoint(0,x);}
            inline void SetY(double y)    {SetPoint(1,y);}
            inline void SetZ(double z)    {SetPoint(2,z);}
            
#pragma mark • Operator
            inline const double operator[](unsigned int i) const {return (i < fx.size()) ? fx[i] : throw std::out_of_range("[DST::Math::points[]] out of range.");}
            
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
            
#pragma mark • Accessor
            inline const std::vector<double>& GetPoint() const {return fx;}
            inline size_t size() const {return fx.size();}
            inline double X() const {return operator[](0);}
            inline double Y() const {return operator[](1);}
            inline double Z() const {return operator[](2);}
            double R()     const;
            double Theta() const;
            double Phi()   const;
            
            static double precision;
            
#pragma mark • Debug
            virtual std::string Dump() const;
            static bool debug;
        };
        
#pragma mark - vector2D class definition
        class vector2D
        {
#pragma mark • protected memeber
        protected:
            double flength;
            double fphi;

            void _phi();
         
         public:
#pragma mark • ctor/dtor
            vector2D();
            vector2D(const point&);
            vector2D(const point&, const point&);
            vector2D(double, double);
            vector2D(double);
            vector2D(const vector2D&);
         
            virtual ~vector2D();
         
#pragma mark • Modifier
            inline void SetPhi   (double _t){fphi    = _t; _phi();}
            inline void SetLength(double _l){flength = fabs(_l);}
         
#pragma mark • Accessor
            inline double Phi      () const {return fphi;}
            inline double Length   () const {return flength;}
            inline virtual double Theta() const {return acos(-1)/2.;}
            
            inline virtual double X() const {return flength*cos(fphi);}
            inline virtual double Y() const {return flength*sin(fphi);}
            inline virtual double Z() const {return 0;}
            
            virtual DST::Math::point EndPoint() const {return point(2, X(), Y());}
         
#pragma mark • Opperator
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
            
#pragma mark • Dump
            virtual std::string Dump() const;
            
         };
        
#pragma mark - vector3D class definition
        class vector3D: public vector2D
        {
#pragma mark • protected memeber
        protected:
            double ftheta;

            void _theta();
            
        public:
#pragma mark • ctor/dtor
            vector3D();
            vector3D(const point&);
            vector3D(const point&,const point&);
            vector3D(double, double, double);
            vector3D(double, double);
            vector3D(const vector3D&);
            vector3D(const vector2D&);
            
            virtual ~vector3D();
            
#pragma mark • Modifier
            inline void SetTheta  (double _p){ftheta   = _p;}
            
#pragma mark • Accessor
            inline double Theta() const {return ftheta;}
            
            inline virtual double X() const {return flength*sin(ftheta)*cos(fphi);}
            inline virtual double Y() const {return flength*sin(ftheta)*sin(fphi);}
            inline virtual double Z() const {return flength*cos(ftheta);}
            
            virtual DST::Math::point EndPoint() const {return point(3, X(), Y(), Z());}
            
#pragma mark • Opperator
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
            
#pragma mark • Dump
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
