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

#include <iostream>
#include <iomanip>
#include <DSTmath/DSTgeo.h>
#include <DSTmath/DSTmath.h>

std::string to_string(const double a_value, const int n)
{
    std::ostringstream sout;
    sout << std::scientific << std::setprecision(n) << a_value;
    return sout.str();
}

std::string to_string(const float a_value, const int n)
{
    return to_string(static_cast<const double>( a_value ), n);
}

namespace DST
{
    namespace Math
    {
#pragma mark - point class implementation
        double point::precision = std::numeric_limits<double>::min()*1000.;
        bool   point::debug = false;
        
#pragma mark • ctor/dtor
        /**
         *  @brief Default constructor
         *  @details Create floating point coordinates
         */
        point::point()
        {
            fx = std::vector<double>();
        }
        
        /**
         *  @brief Default constructor
         *  @details Create floating point coordinates with n dimension
         *  @param n number of dimension
         *  @param x floating point cartesian coordinates
         */
        point::point(unsigned int n, double x, ... )
        {
            
            double arg = x;
            fx.push_back(arg);
            
            va_list ap;
            va_start(ap, x);
            
            for(unsigned int naxe = 1; naxe < n; naxe++)
            {
                arg = va_arg(ap, double);
                
                fx.push_back(arg);
            }
            va_end(ap);
        }
        
        /**
         *  @brief Copy constructor
         */
        point::point(const point& p)
        {
            for(unsigned int i = 0; i < p.fx.size(); i++)
                fx.push_back(p.fx[i]);
        }
        
        point::~point()
        {
            fx.clear();
        }
        
#pragma mark • Modifier
        /**
         *  @brief Assign coordinates
         *  @details Assign floating point cartesian coordinate to this point
         *
         *  @param x   floating point cartesian coordinate.
         *  @note The number of floating point coordinate given as parameters shall equals the dimension of this point.
         */
        void point::SetPoints( double x, ...)
        {
            double arg = x;
            if(fx.size() < 1)
                fx.push_back(arg);
            else
                fx[0] = arg;
            
            va_list ap;
            va_start(ap, x);
            
            for(unsigned int naxe = 1; naxe < fx.size(); naxe++)
            {
                arg = va_arg(ap, double);
                
                fx[naxe] = arg;
            }
            va_end(ap);
        }
        
        /**
         *  @brief Assign coordinates
         *  @details Assign floating point cartesian coordinate to this point by copying an other point coordinates
         *
         *  @param p   Cartesian coordinate to assign to this.
         *  @note If the dimension of p isn't the same as this, the dimension of this is extended to the dimension of p.
         */
        void point::SetPoints( const point& p )
        {
            for(unsigned int i = 0; i < p.fx.size(); i++)
            {
                if(fx.size() <= i)
                    fx.push_back(p.fx[i]);
                else
                    fx[i] = p.fx[i];
            }
            
        }
        
#pragma mark • Accessor
        
        /**
         *  @brief Get distance from the origine
         *  @return \f$\sqrt{\sum_i x_i^2}\f$
         */
        double point::R() const
        {
            double radius = 0;
            for(unsigned int i = 0; i < fx.size(); i++)
                radius += fx[i]*fx[i];
            
            return sqrt(radius);
        }
        
        /**
         *  @brief Get azimuth angle
         *  @return \f$\atan\left(\frac{y}{x}\right)\f$
         */
        double point::Phi() const
        {
            if(fx.size() >= 2)
                return (atan2(fx[1], fx[0]) > 0 ) ? atan2(fx[1], fx[0]) : 2.*acos(-1) - atan2(fx[1], fx[0]);
            
            return 0;
        }
        
        /**
         *  @brief Get Elevation angle
         *  @return \f$\sqrt{\sum_i x_i^2}\f$
         */
        double point::Theta() const
        {
            if(fx.size() >= 3)
                return acos(fx[2]/R());
            
            return 0;
        }
        
#pragma mark • Operator
        
        /**
         *  @brief Not equal comparator
         *  @param p cartesian coordinate to compare to this
         *  @return true if each floating point coordinate of p differs from this
         */
        bool point::operator!=(const point& p)  const
        {
            return !(operator==(p));
        }
        
        /**
         *  @brief Equal comparator
         *  @param p cartesian coordinate to compare to this
         *  @return true if each floating point coordinate of p equal those of this
         */
        bool point::operator==(const point& p) const
        {
            bool isEqual = true;
            
            if(fx.size() != p.fx.size())
                return !isEqual;
            
            for(unsigned int i = 0; i < fx.size(); i++)
                isEqual &= ( (fx[i] == p.fx[i]) || (fabs(fx[i] - p.fx[i]) <= point::precision) );
    
            return isEqual;
        }
        
        /**
         *  @brief less or equal comparator
         *  @param p cartesian coordinate to compare to this
         *  @return true if this is closer to the cartesian axis than p
         *  @note DST::Math::point is ordered by distance to the main axis, dimension per dimension. Not by distance to the center.
         */
        bool point::operator<=(const point& p)  const
        {
            if(operator==(p))
                return true;
            
            unsigned int imax = (fx.size() <= p.fx.size()) ? fx.size() : p.fx.size();
            
            for(unsigned int i = imax; i != 0; i--)
            {
                double this_radius = 0;
                double p_radius    = 0;
                
                for(unsigned int k = 0; k < imax; k++)
                {
                    this_radius += fx[k]*fx[k];
                    p_radius    += p.fx[k]*p.fx[k];
                }
                
                this_radius = sqrt(this_radius);
                p_radius    = sqrt(p_radius);
                
                if(fabs(this_radius - p_radius) > point::precision)
                    return this_radius - p_radius < - point::precision;
                
                if(i < 1)
                    continue;
                
                if( fabs(fx[i-1] - p.fx[i-1]) > point::precision)
                    return fx[i-1] - p.fx[i-1] < - point::precision;
            }
            
            
            return fx[0] - p.fx[0]   < - point::precision;
        }
        
        /**
         *  @brief less comparator
         *  @param p cartesian coordinate to compare to this
         *  @return true if this is closer to the cartesian axis than p
         *  @note DST::Math::point is ordered by distance to the main axis, dimension per dimension. Not by distance to the center.
         */
        bool point::operator< (const point& p)  const
        {
            return operator<=(p) && !(operator==(p));
        }
        
        /**
         *  @brief greater or equal comparator
         *  @param p cartesian coordinate to compare to this
         *  @return true if this is further away to the cartesian axis than p
         *  @note DST::Math::point is ordered by distance to the main axis, dimension per dimension. Not by distance to the center.
         */
        bool point::operator>=(const point& p)  const
        {
            if(operator==(p))
                return true;
            
            unsigned int imax = (fx.size() <= p.fx.size()) ? fx.size() : p.fx.size();
            
            for(unsigned int i = imax; i != 0; i--)
            {
                double this_radius = 0;
                double p_radius    = 0;
                
                for(unsigned int k = 0; k < imax; k++)
                {
                    this_radius += fx[k]*fx[k];
                    p_radius    += p.fx[k]*p.fx[k];
                }
                
                this_radius = sqrt(this_radius);
                p_radius    = sqrt(p_radius);
                
                if(fabs(this_radius - p_radius) > point::precision)
                    return this_radius - p_radius > point::precision;
                
                if(i < 1)
                    continue;
                
                if( fabs(fx[i-1] - p.fx[i-1]) > point::precision)
                    return fx[i-1] - p.fx[i-1] < - point::precision;
            }
            
            
            return fx[0] - p.fx[0]   > - point::precision;
        }
        
        /**
         *  @brief greater comparator
         *  @param p cartesian coordinate to compare to this
         *  @return true if this is further away to the cartesian axis than p
         *  @note DST::Math::point is ordered by distance to the main axis, dimension per dimension. Not by distance to the center.
         */
        bool point::operator> (const point& p)  const
        {
            return operator>=(p) && !(operator==(p));
        }
        
        /**
         *  Assignement operator
         *  @param p point to copy
         */
        point& point::operator=(const point& p)
        {
            if(this != &p)
            {
                point tmp(p);
                
                fx.swap(tmp.fx);
            }
            return *this;
        }
        
        void point::operator+=(const DST::Math::point& p)
        {
            unsigned int n_elem = (fx.size() <= p.fx.size())? fx.size(): p.fx.size();
            for(unsigned int i = 0; i < n_elem; i++)
                fx[i]+=p.fx[i];
        }
        
        void point::operator-=(const DST::Math::point& p)
        {
            unsigned int n_elem = (fx.size() <= p.fx.size())? fx.size(): p.fx.size();
            for(unsigned int i = 0; i < n_elem; i++)
                fx[i]-=p.fx[i];
        }
        
        void point::operator*=(const DST::Math::point& p)
        {
            unsigned int n_elem = (fx.size() <= p.fx.size())? fx.size(): p.fx.size();
            for(unsigned int i = 0; i < n_elem; i++)
                fx[i]*=p.fx[i];
        }
        
        void point::operator/=(const DST::Math::point& p)
        {
            unsigned int n_elem = (fx.size() <= p.fx.size())? fx.size(): p.fx.size();
            for(unsigned int i = 0; i < n_elem; i++)
                fx[i]/=p.fx[i];
        }
        
        void point::operator+=(const double p)
        {
            for(unsigned int i = 0; i < fx.size(); i++)
                fx[i]+=p;
        }
        
        void point::operator-=(const double p)
        {
            for(unsigned int i = 0; i < fx.size(); i++)
                fx[i]-=p;
        }
        void point::operator*=(const double p)
        {
            for(unsigned int i = 0; i < fx.size(); i++)
                fx[i]*=p;
        }
        void point::operator/=(const double p)
        {
            for(unsigned int i = 0; i < fx.size(); i++)
                fx[i]/=p;
        }
        
        void point::operator+=(const int p)
        {
            operator+=(static_cast<double>( p ));
        }
        
        void point::operator-=(const int p)
        {
            operator-=(static_cast<double>( p ));
        }
        
        void point::operator*=(const int p)
        {
            operator*=(static_cast<double>( p ));
        }
    
        void point::operator/=(const int p)
        {
            operator+=(static_cast<double>( p ));
        }
        

        
#pragma mark • Dump
        
        std::string point::Dump() const
        {
            std::string sdump = std::string();
            
            if(DST::Math::point::debug)
                sdump += std::string("\033[31m");
            
            sdump += std::string("(");
            
            for(unsigned int i = 0 ; i < fx.size()-1; i++)
            {
                if( DST::Math::sgn( fx[i] ) >= 0 )
                    sdump += " ";
                
                sdump += to_string(fx[i],4)+std::string(" , ");
            }
            if( DST::Math::sgn( fx[fx.size()-1] ) >= 0 )
                sdump += " ";
            
            sdump += to_string(fx[fx.size()-1],4)+std::string(") r = ")+to_string(R(),4);
            
            if(DST::Math::point::debug)
                sdump += std::string("\033[0m");
                
            return sdump;
        }
        
#pragma mark - vector2D class implementation
#pragma mark • ctor/dtor
        
        /**
         *  @brief Default constructor
         *  @details Create unitary vector aligned with the \f$x\f$ axis
         *  @note The origin of the vector is at the origin
         */
        vector2D::vector2D():flength(1.),fphi(0.)
        {};
        
        /**
         *  @brief Constructor
         *  @param p The angle with the \f$x\f$ axis
         *  @param l The vector length
         *  @details Create vector with length and orientation
         *  @note The origin of the vector is at the origin
         */
        vector2D::vector2D(double l, double p):flength(fabs(l)), fphi(p)
        {};
        
        /**
         *  @brief Constructor
         *  @param p The angle with the \f$x\f$ axis
         *  @details Create unitary vector with orientation
         *  @note The origin of the vector is at the origin
         */
        vector2D::vector2D(double p):flength(1.), fphi(p)
        {};
        
        /**
         *  @brief Copy constructor
         *  @param v vector to be copied
         *  @note The origin of the vector is at the origin
         */
        vector2D::vector2D(const vector2D& v): flength(v.flength), fphi(v.fphi)
        {};
        
        /**
         *  @brief Constructor
         *  @details Construc vector by coordinate of its end-point
         *  @param p cartesian coordinate of the vector end-point
         *  @note The origin of the vector is at the origin
         */
        vector2D::vector2D(const point& p):flength(1),fphi(0)
        {
            if(p.size() <= 1)
            {
                return;
            }
            
            flength = sqrt( p.X()*p.X() + p.Y()*p.Y() );
            fphi    = atan2(p.Y(), p.X());
            
                
        }
        
        /**
         *  @brief destructor
         */
        vector2D::~vector2D(){};
        
#pragma mark • Opperator
        
        bool vector2D::operator!=(const vector2D& v) const
        {
            return fabs(flength - v.flength) > point::precision || fabs(fphi - v.fphi) > point::precision;
        }
        
        bool vector2D::operator==(const vector2D& v) const
        {
            return fabs(flength - v.flength) < point::precision && fabs(fphi - v.fphi) < point::precision;
        }
        
        bool vector2D::operator<=(const vector2D& v) const
        {
            if(operator==(v))
                return true;
            
            if( fabs(flength - v.flength) < point::precision )
                return fphi - v.fphi < - point::precision;
            
            return flength - v.flength < - point::precision;
        }
        
        bool vector2D::operator< (const vector2D& v) const
        {
            return operator<=(v) && operator!=(v);
        }

        bool vector2D::operator>=(const vector2D& v) const
        {
            if(operator==(v))
                return true;
            
            if( fabs(flength - v.flength) < point::precision )
                return fphi - v.fphi > point::precision;
            
            return (flength - v.flength) > point::precision;
        }

        bool vector2D::operator> (const vector2D& v) const
        {
            return operator>=(v) && operator!=(v);
        }
        
        /**
         *  Assignement operator
         *  @param v vector to copy
         */
        vector2D vector2D::operator=(const vector2D& v)
        {
            if(this != &v)
            {
                vector2D tmp(v);
                
                std::swap(flength, tmp.flength);
                std::swap(fphi, tmp.fphi);
            }
            return *this;
        }
        
        /**
         *  @brief Addition
         */
        void vector2D::operator+=(const vector2D& v)
        {
            double fx = X()+v.X();
            double fy = Y()+v.Y();
            flength = sqrt(fx*fx + fy*fy);
            fphi    = acos(fx/flength);
        }
        
        /**
         *  @brief Substraction
         */
        void vector2D::operator-=(const vector2D& v)
        {
            double fx = X() - v.X();
            double fy = Y() - v.Y();
            flength = sqrt(fx*fx + fy*fy);
            fphi    = acos(fx/flength);
        }
        
        /**
         *  @brief Scaling
         */
        void vector2D::operator*=(const double s)
        {
            flength *= s;
        }
        
        /**
         *  @brief Scaling
         */
        void vector2D::operator/=(const double s)
        {
            flength /= s;
        }
        
#pragma mark • Dump
        
        std::string vector2D::Dump() const
        {
            std::string sdump = std::string();
            
            if(DST::Math::point::debug)
                sdump += std::string("\033[31m");
            
            sdump += std::string("r = ");
            sdump += to_string(flength,4);
            
            sdump += std::string(" | phi = ");
            
            if( DST::Math::sgn( fphi ) >= 0 )
                    sdump += " ";
            
            sdump += to_string(fphi,4);
            
            if(DST::Math::point::debug)
                sdump += std::string("\033[0m");
            
            return sdump;
        }
        
#pragma mark - vector3D class implementation
#pragma mark • ctor/dtor
        /**
         *  @brief Default constructor
         *  @details Create unitary vector aligned with the \f$x\f$ axis
         *  @note The origin of the vector is at the origin
         */
        vector3D::vector3D(): vector2D(), ftheta(0.)
        {};
        
        /**
         *  @brief Constructor
         *  @details Construc vector by coordinate of its end-point
         *  @param p cartesian coordinate of the vector end-point
         *  @note The origin of the vector is at the origin
         */
        vector3D::vector3D(const point& p)
        {
            if(p.size() < 3)
                return;
            
            flength = sqrt( p.X()*p.X() + p.Y()*p.Y() );
            fphi    = atan2(p.Y(), p.X());
            ftheta  = acos(p.Z());
        }
        
        /**
         *  @brief Constructor
         *  @details Construc vector with length and orientation
         *  @param l vector length
         *  @param t vector elevation angle
         *  @param p vector azimuth angle
         *  @note The origin of the vector is at the origin
         */
        vector3D::vector3D(double l, double p, double t):vector2D(l,p), ftheta(t)
        {};
        
        vector3D::vector3D(double p, double t):vector2D(p), ftheta(t)
        {};
        
        vector3D::vector3D(const vector3D& v):vector2D()
        {
            flength = v.flength;
            fphi    = v.fphi;
            ftheta  = v.ftheta;
        }
        
        vector3D::vector3D(const vector2D& v):vector2D(v.Length(), v.Phi())
        {
            ftheta  = v.Theta();
        }
            
        vector3D::~vector3D(){};
        
#pragma mark • Opperator
        
        bool vector3D::operator!=(const vector3D& v) const
        {
            return vector2D::operator!=( static_cast<vector2D>(v) ) || fabs(ftheta - v.ftheta) > point::precision;
        }
        
        bool vector3D::operator==(const vector3D& v) const
        {
            return vector2D::operator==( static_cast<vector2D>(v) ) && fabs(ftheta - v.ftheta) <= point::precision;
        }
        
        bool vector3D::operator<=(const vector3D& v) const
        {
            if(operator==(v))
                return true;
            
            if(fabs(flength - v.flength) > point::precision)
                return static_cast<int>((flength - v.flength)/fabs(flength - v.flength)) <  0 ;
            
            double this_ad = acos(sin(ftheta)  * cos(fphi));
            double v_ad    = acos(sin(v.ftheta)* cos(v.fphi));
            
            if(fabs(this_ad - v_ad) > point::precision)
                return static_cast<int>((this_ad-v_ad)/fabs(this_ad-v_ad)) < 0;
            
            if(fabs(fphi - v.fphi) > point::precision)
                return static_cast<int>((fphi - v.fphi)/fabs(fphi - v.fphi)) < 0 ;
            
            return fabs(ftheta - v.ftheta) > point::precision && ftheta*180/acos(-1.) < v.ftheta*180/acos(-1.);
                
        }
        
        bool vector3D::operator< (const vector3D& v) const
        {
            return operator<=(v) && operator!=(v);
        }
        
        bool vector3D::operator>=(const vector3D& v) const
        {
            if(operator==(v))
                return true;
            
            if(fabs(flength - v.flength) > point::precision)
                return flength - v.flength > point::precision;
            
            double this_ad = acos(sin(ftheta)  * cos(fphi));
            double v_ad    = acos(sin(v.ftheta)* cos(v.fphi));
            
            if(fabs(this_ad - v_ad) > point::precision)
                return this_ad - v_ad > point::precision;
            
            if(fabs(fphi - v.fphi) > point::precision)
                return fphi - v.fphi > point::precision;
            
            return ftheta - v.ftheta > point::precision;
            
        }
        
        bool vector3D::operator> (const vector3D& v) const
        {
            return operator>=(v) && operator!=(v);
        }
        
        /**
         *  Assignement operator
         *  @param v vector to copy
         */
        vector3D vector3D::operator=(const vector3D& v)
        {
            if(this != &v)
            {
                vector3D tmp(v);
                
                std::swap(flength, tmp.flength);
                std::swap(fphi, tmp.fphi);
                std::swap(ftheta, tmp.ftheta);
            }
            return *this;
        }
        
        /**
         *  @brief Addition
         */
        void vector3D::operator+=(const vector3D& v)
        {
            double fx = X()+v.X();
            double fy = Y()+v.Y();
            double fz = Z()+v.Z();
            flength = sqrt(fx*fx + fy*fy + fz*fz);
            ftheta  = acos(fz/flength);
            fphi    = atan2(fy,fx);
        }
        
        /**
         *  @brief Substraction
         */
        void vector3D::operator-=(const vector3D& v)
        {
            double fx = X()-v.X();
            double fy = Y()-v.Y();
            double fz = Z()-v.Z();
            flength = sqrt(fx*fx + fy*fy + fz*fz);
            ftheta  = acos(fz/flength);
            fphi    = atan2(fy,fx);
        }
        
        /**
         *  @brief Substraction
         */
        void vector3D::operator^=(const vector3D& v)
        {
            double fx = Y()*v.Z() - Z()*v.Y();
            double fy = Z()*v.X() - X()*v.Z();
            double fz = X()*v.Y() - Y()*v.X();
            
            flength = sqrt(fx*fx + fy*fy + fz*fz);
            ftheta  = acos(fz/flength);
            fphi    = atan2(fy,fx);
        }
        
#pragma mark • Dump
        
        std::string vector3D::Dump() const
        {
            std::string sdump = vector2D::Dump();
            
            if(DST::Math::point::debug)
                sdump += std::string("\033[31m");
            
            sdump += " | theta = ";
            
            if( DST::Math::sgn( ftheta ) >= 0 )
                sdump += " ";
            
            sdump += to_string(ftheta, 4);
            
            if(DST::Math::point::debug)
                sdump += std::string("\033[0m");
            
            return sdump;
        }
        
    }
}

#pragma mark • operator

DST::Math::point operator+(const DST::Math::point& p1, const DST::Math::point& p2)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out+=p2;
    
    return p_out;
}

DST::Math::point operator-(const DST::Math::point& p1, const DST::Math::point& p2)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out-=p2;
    
    return p_out;
}

DST::Math::point operator*(const DST::Math::point& p1, const DST::Math::point& p2)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out*=p2;
    
    return p_out;
}

DST::Math::point operator/(const DST::Math::point& p1, const DST::Math::point& p2)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out/=p2;
    
    return p_out;
}

DST::Math::point operator+(const DST::Math::point& p1, const double s)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out+=s;
    
    return p_out;
}

DST::Math::point operator-(const DST::Math::point& p1, const double s)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out-=s;
    
    return p_out;
}

DST::Math::point operator*(const DST::Math::point& p1, const double s)

{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out*=s;
    
    return p_out;
}

DST::Math::point operator/(const DST::Math::point& p1, const double s)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out/=s;
    
    return p_out;
}

DST::Math::point operator+(const double s, const DST::Math::point& p1)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out+=s;
    
    return p_out;
}

DST::Math::point operator-(const double s, const DST::Math::point& p1)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out-=s;
    p_out*=-1;
    
    return p_out;
}

DST::Math::point operator*(const double s, const DST::Math::point& p1)
{
    DST::Math::point p_out = DST::Math::point(p1);
    p_out*=s;
    
    return p_out;
}

DST::Math::point operator/(const double s, const DST::Math::point& p1)
{
    DST::Math::point p_out = DST::Math::point();
    for(unsigned int n = 0; n < p1.size(); n++)
        p_out.SetPoint(n,1./p1[n]);
    
    p_out*=s;
    
    return p_out;
}

double operator*(const DST::Math::vector2D& v1, const DST::Math::vector2D& v2)
{
    double cosDL = cos(v1.Theta())*cos(v2.Theta())+sin(v1.Theta())*sin(v2.Theta())*cos(v1.Phi()-v2.Phi());
    return v1.Length()*v2.Length()*cosDL;
}

double operator*(const DST::Math::vector2D& v1, const DST::Math::vector3D& v2)
{
    return v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
}

double operator*(const DST::Math::vector3D& v1, const DST::Math::vector2D& v2)
{
    return v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
}

double operator*(const DST::Math::vector3D& v1, const DST::Math::vector3D& v2)
{
    double cosDL = cos(v1.Theta())*cos(v2.Theta())+sin(v1.Theta())*sin(v2.Theta())*cos(v1.Phi()-v2.Phi());
    return v1.Length()*v2.Length()*cosDL;
}

DST::Math::vector2D operator*(const DST::Math::vector2D& v1, const double s)
{
    DST::Math::vector2D v2(v1);
    v2*=s;
    return v2;
}

DST::Math::vector2D operator*(const double s, const DST::Math::vector2D& v1)
{
    DST::Math::vector2D v2(v1);
    v2*=s;
    return v2;
}

DST::Math::vector3D operator*(const DST::Math::vector3D& v1, const double s)
{
    DST::Math::vector3D v2(v1);
    v2*=s;
    return v2;
}

DST::Math::vector3D operator*(const double s, const DST::Math::vector3D& v1)
{
    DST::Math::vector3D v2(v1);
    v2*=s;
    return v2;
}

DST::Math::vector3D operator^(const DST::Math::vector2D& v1, const DST::Math::vector2D& v2)
{
    DST::Math::vector3D vv1(v1);
    DST::Math::vector3D vv2(v2);
    vv1^=vv2;
    
    return vv1;
}

DST::Math::vector3D operator^(const DST::Math::vector2D& v1, const DST::Math::vector3D& v2)
{
    DST::Math::vector3D vv1(v1);
    vv1^=v2;
    
    return vv1;
}

DST::Math::vector3D operator^(const DST::Math::vector3D& v1, const DST::Math::vector2D& v2)
{
    DST::Math::vector3D vv1(v1);
    DST::Math::vector3D vv2(v2);
    vv1^=vv2;
    
    return vv1;
}

DST::Math::vector3D operator^(const DST::Math::vector3D& v1, const DST::Math::vector3D& v2)
{
    DST::Math::vector3D vv1(v1);

    vv1^=v2;
    
    return vv1;
}
