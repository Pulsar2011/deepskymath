//
//  DSTmath.h
//  DeepSkyTools
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#ifndef __DST_MATH_CORE__
#define __DST_MATH_CORE__
#include<cmath>
#include<limits>
#include<vector>
#include <stdexcept>
#include <string>
#if __cplusplus >= 202002L
#include <span>
#endif

namespace DST
{
    std::string to_string(const double&, const int&);
    std::string to_string(const float& , const int&);

    namespace Math
    {
        
#pragma mark - MathCore class definition
        /**
         *  @class MathCore DSTmath.h "DSTmath/DSTmath.h"
         *  @brief Definition of basic math
         *  @ MathCore class is a re-defintion of basic mathematical function required in other class from Math::namspace. Most of the function are re-implementation of ROOT::Math functions.
         */
        class MathCore
        {
            protected:
            
                static const double kBig;
                static const double kBiginv;
                static const double LS2PI;
            
                static const double A[];
                static const double B[];
                static const double C[];
                static const double P[];
                static const double Q[];
                static const double STIR[];
            
                static double igam(double, double);
                static double lgam(double);
            
                static double igamc(double, double);
            
                static double Polynomialeval(double, const double*, unsigned int);
                static double Polynomial1eval(double, const double*, unsigned int);
            
                static double small(double z, double x);
            
                static double stirf( double);
            
            public:
            
                static double Pi();
                static double Gamma(double, double);
                static double Gamma(double);
                static double J1(double);
        };
    
#pragma mark - Func class definition
        
        /*!
            *  @class polynom DSTmath.h "DSTmath/DSTmath.h"
            *  @brief Definition of polynomial function
            *  @details polynom class is a re-defintion of polynomial function required in other class from Math::namspace. Most of the function are re-implementation of Numerical Recipes functions.
        */
        class polynom : public MathCore
        {
            protected:
                static std::vector<double> chebfit(std::vector<double>, double, double);

            public:
                static std::vector<double> chebfit(double (*func)(double), const double&, const double&, const unsigned int&);
                static std::vector<float>  chebfit(float  (*func)(float), const float&, const float&, const unsigned int&);

                static double chebev(const double&, const std::vector<double>&, const double&, const double&, unsigned int n = 0);
                static float  chebev(const float& , const std::vector<float>&  , const float& , const float& , unsigned int n = 0);
                static void   chebder(std::vector<double>, double, double, std::vector<double>&);
                static void   chebder(std::vector<float>, float, float, std::vector<float>&);

                static void chebinv(double, double&, std::vector<double>, double, double, double expsillon = std::numeric_limits<double>::epsilon(), size_t max_iter = 100000);

#if __cplusplus >= 199711L
                static double chebev2(double *, std::vector<double>, std::vector<double>, double *, double *, std::vector<unsigned int> order = {0,0});
                static float  chebev2(float  *, std::vector<float> , std::vector<float> , float  *, float  *, std::vector<unsigned int> order = {0,0});
#else
                static double chebev2(double *, std::vector<double>, std::vector<double>, double *, double *);
                static float  chebev2(float  *, std::vector<float> , std::vector<float> , float  *, float  *);

                static double chebev2(double *, std::vector<double>, std::vector<double>, double *, double *, std::vector<unsigned int> order);
                static float  chebev2(float  *, std::vector<float> , std::vector<float> , float  *, float  *, std::vector<unsigned int> order);
#endif
            
                static double chebev2(double *, std::vector<double>, double *, double *, std::vector<unsigned int> order);
                static float  chebev2(float  *, std::vector<float>,  float  *, float  *, std::vector<unsigned int> order);

                static double polyev(const double&, const std::vector<double>&);
                static double polyev(const double&, const double *, const size_t&);
                static double polyev(const double&, const double *, const unsigned int&);
                static float  polyev(const float& , const std::vector<float>&);
                static float  polyev(const float& , const float *, const size_t&);
                static float  polyev(const float& , const float *, const unsigned int&);
#if __cplusplus >= 202002L
                static double polyev(const double&, std::span<const double>);
                static float  polyev(const float& , std::span<const float>);
#endif
        };
        
#pragma mark - Other mathematical function
        
#if __cplusplus >= 199711L
        /**
         *  @brief Obtain the sign of \f$x\f$
         *  @param x         Numerical value one want to check the sign
         *  @param is_signed Either the type of \f$x\f$ is signed or not
         *  @return 1 if \f$x\f$ null or positive, 0 otherwise
         */
        template <typename T> constexpr
        int sgn(T x, std::false_type is_signed)
        {
            return static_cast<T>(0) < x;
        }
        
        /**
         *  @brief Obtain the sign of \f$x\f$
         *  @param x         Numerical value one want to check the sign
         *  @param is_signed Either the type of \f$x\f$ is signed or not
         *  @return 1 if \f$x\f$ positive, 0 if \f$x = 0\f$ and -1 if \f$x\f$ negative
         */
        template <typename T> constexpr
        int sgn(T x, std::true_type is_signed)
        {
            return (static_cast<T>(0) < x) - (x < static_cast<T>(0));
        }
        
        /**
         *  @brief Obtain the sign of \f$x\f$
         *  @param x         Numerical value one want to check the sign
         *  @return 1 if \f$x\f$ positive, 0 if \f$x = 0\f$ and -1 if \f$x\f$ negativ
         */
        template <typename T>  constexpr
        int sgn(T x)
        {
            return sgn(x, std::is_signed<T>());
        }
#else
        /**
         *  @brief Obtain the sign of \f$x\f$
         *  @param x         Numerical value one want to check the sign 
         *  @param is_signed Either the type of \f$x\f$ is signed or not
         *  @return 1 if \f$x\f$ positiv, 0 if \f$x = 0\f$
         */
        template <typename T>
        int sgn(T x)
        {
            return (static_cast<T>(0) < x) - (x < static_cast<T>(0));
        }
#endif

#ifndef SQT2PI
    #define SQT2PI   std::sqrt(2.*DST::Math::MathCore::Pi()) /* sqrt(2*pi) */
#endif


        namespace function
        {
            template <typename T>
            constexpr bool is_allowed_math_type = std::is_arithmetic_v<T> || 
                std::is_same_v<T,double> ||
                std::is_same_v<T,float>; 

            /**
             * @brief Eval the value of a Gaussian function at \f$G(x) = \frac{A}{\sigma\sqrt{\2\pi}\exp\left(\frac{-(x-m_)^2}{2\sigma^2}\right)\f$.
             * 
             * @param x_ Position at which to evaluate the Gaussian function
             * @param A_ Amplitude of the Gaussian function
             * @param m_ Mean value of the normale distribution
             * @param s_ Standard deviation of the normale distribution
             * @return double value sampled from the normale distribution
             */
            template<typename U>
            U gauss ( const U& _x, const U& A, const U& m_, const U& s_ )
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::normal_distribution] Unsuported datatype or data structure\033[0m");

                const U x = _x - m_;
                const U s2 = s_ * s_;
                const U norm = 1.0 / (s_ * static_cast<U>(SQT2PI) );
                const U exp_part = std::exp(-0.5 * (x * x) / s2);
                return A * norm * exp_part;
            }

            /** 
             * @brief Evaluate the integral of a Gaussian function from \f$x_-\frac{dx}{2}\f$ to \f$x_+\frac{dx}{2}\f$
             * 
             * @param _x Position at which to evaluate the Gaussian function
             * @param A Amplitude of the Gaussian function
             * @param m_ Mean value of the Gaussian function
             * @param s_ Standard deviation of the Gaussian function
             * @param dx_ Width of the integration interval
             * @return value of the integral of the Gaussian function from \f$x_1\f$ to \f$x_2\f$
             */
            template<typename U>
            U integrated_gauss ( const U& _x, const U& A, const U& m_, const U& s_, const U& dx_ )
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::integrated_gauss] Unsuported datatype or data structure\033[0m");

                const U x1 = (_x-dx_/2.) - m_;
                const U x2 = (_x+dx_/2.) - m_;
                const U norm = 0.5 * A ;
                return norm * (std::erf(x2 / (s_ * std::sqrt(2.0))) - std::erf(x1 / (s_ * std::sqrt(2.0))) );
            }

            /** 
             * @brief Evaluate the integral of a Gaussian function from \f$x_-\frac{dx}{2}\f$ to \f$x_+\frac{dx}{2}\f$
             * 
             * @param _x Position at which to evaluate the Gaussian function
             * @param A Amplitude of the Gaussian function
             * @param m_ Mean value of the Gaussian function
             * @param s_ Standard deviation of the Gaussian function
             * @param dx_ Width of the integration interval
             * @return value of the integral of the Gaussian function from \f$x_1\f$ to \f$x_2\f$
             */
            template<typename U>
            U igauss ( const U& _x, const U& A, const U& m_, const U& s_, const U& dx_ )
            {
                return integrated_gauss(_x, A, m_, s_, dx_);
            }

            /**
             * @brief Eval the value of a 2D symmetric Gaussian function \f$ G(x,y) = \frac{A}{2\pi\sigma^2}\exp\left(-\frac{1}{2}\left[\frac{(x-m_x)^2}{\sigma^2} + \frac{(y-m_y)^2}{\sigma^2}\right]\right)\f$
             * 
             * @param x_ Position at which to evaluate the Gaussian function
             * @param y_ Position at which to evaluate the Gaussian function
             * @param A_ Amplitude of the Gaussian function
             * @param m_x Mean value of the normale distribution in x-direction
             * @param m_y Mean value of the normale distribution in y-direction
             * @param s_ Standard deviation of the normale distribution in x-direction
             * @return double value sampled from the normale distribution
             */
            template<typename U>
            U gauss2D ( const U& _x, const U& _y, const U& A, const U& m_x, const U& m_y, const U& s_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::normal_distribution] Unsuported datatype or data structure\033[0m");

                const U r = std::sqrt( (_x - m_x)*(_x - m_x) + (_y - m_y)*(_y - m_y) );
                const U norm = 1.0 / (s_ * static_cast<U>(SQT2PI) );
                return gauss<U>( r, A, 0.0, s_ )*norm;
            }

            /** 
             * @brief Evaluate the integral of a Gaussian function from \f$x_-\frac{dx}{2}\f$ to \f$x_+\frac{dx}{2}\f$ and \f$y_-\frac{dy}{2}\f$ to \f$y_+\frac{dy}{2}\f$
             * 
             * @param _x Position at which to evaluate the Gaussian function
             * @param _y Position at which to evaluate the Gaussian function
             * @param A Amplitude of the Gaussian function
             * @param m_x Mean value of the Gaussian function
             * @param m_y Mean value of the Gaussian function
             * @param s_ Standard deviation of the Gaussian function
             * @param dx_ Width of the integration interval
             * @param dy_ Width of the integration interval
             * @return value of the integral of the Gaussian function from \f$x_1\f$ to \f$x_2\f$
             */
            template<typename U>
            U integrated_gauss2D ( const U& _x, const U& _y, const U& A, const U& m_x, const U& m_y, const U& s_, double dx_, double dy_ )
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::normal_distribution] Unsuported datatype or data structure\033[0m");

                const U iG_x = integrated_gauss<U>(_x, 2., m_x, s_, dx_);
                const U iG_y = integrated_gauss<U>(_y, 2., m_y, s_, dy_);
                return A * 0.25 * iG_x * iG_y;
            }

            /** 
             * @brief Evaluate the integral of a Gaussian function from \f$x_-\frac{dx}{2}\f$ to \f$x_+\frac{dx}{2}\f$ and \f$y_-\frac{dy}{2}\f$ to \f$y_+\frac{dy}{2}\f$
             * 
             * @param _x Position at which to evaluate the Gaussian function
             * @param _y Position at which to evaluate the Gaussian function
             * @param A Amplitude of the Gaussian function
             * @param m_x Mean value of the Gaussian function
             * @param m_y Mean value of the Gaussian function
             * @param s_ Standard deviation of the Gaussian function
             * @param dx_ Width of the integration interval
             * @param dy_ Width of the integration interval
             * @return value of the integral of the Gaussian function from \f$x_1\f$ to \f$x_2\f$
             */
            template<typename U>
            U igauss2D ( const U& _x, const U& _y, const U& A, const U& m_x, const U& m_y, const U& s_, double dx_, double dy_ )
            {
                return integrated_gauss2D(_x, _y, A, m_x, m_y, s_, dx_, dy_);
            }

            /**
             * @brief Eval the value of a 2D asymmetric Gaussian function \f$ G(x,y) = \frac{A}{2\pi\sigma_x\sigma_y\sqrt{1-\rho^2}}\exp\left(-\frac{1}{2}\left[\frac{(x-m_x)^2}{(\cos(\theta)\sigma_x)^2} + \frac{(y-m_y)^2}{\sigma_y^2} - \frac{2\rho(x-m_x)(y-m_y)}{\sigma_x\sigma_y}\right]\right)\f$
             * 
             * @param x_ Position at which to evaluate the Gaussian function
             * @param y_ Position at which to evaluate the Gaussian function
             * @param A_ Amplitude of the Gaussian function
             * @param m_x Mean value of the normale distribution in x-direction
             * @param m_y Mean value of the normale distribution in y-direction
             * @param s_x Standard deviation of the normale distribution in x-direction
             * @param s_y Standard deviation of the normale distribution in y-direction
             * @param _theta Orientation of the Gaussian function in radians. The Gaussian function is rotated by an angle \f$\alpha\f$ with respect to the x-axis.
             * @return double value sampled from the normale distribution
             */
            template<typename U>
            U gauss2D ( const U& _x, const U& _y, const U& A, const U& m_x, const U& m_y, const U& s_x, const U& s_y, U _theta=0)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::normal_distribution] Unsuported datatype or data structure\033[0m");

                const U Ca = 0.5*(std::cos(_theta)*std::cos(_theta)/(s_x*s_x) + std::sin(_theta)*std::sin(_theta)/(s_y*s_y));
                const U Cb = 0.25*(-std::sin(2*_theta)/(s_x*s_x) + std::sin(2*_theta)/(s_y*s_y));
                const U Cc = 0.5*(std::sin(_theta)*std::sin(_theta)/(s_x*s_x) + std::cos(_theta)*std::cos(_theta)/(s_y*s_y));

                const U x = _x - m_x;
                const U y = _y - m_y;
                const U norm = 1.0 / (2.0 * MathCore::Pi() * s_x * s_y);
                const U exp_part = std::exp(-(Ca * x * x - 2.*Cb * x * y + Cc * y * y));
                return A * norm * exp_part;
            }

            /**
             * @brief Evaluate the integral of a 2D asymmetric rotated Gaussian over a rectangular bin of size \f$(dx, dy)\f$ centred on \f$(x, y)\f$.
             *
             * The Gaussian exponent is \f$-(C_a\,\delta x^2 - 2C_b\,\delta x\,\delta y + C_c\,\delta y^2)\f$ where
             * \f$\delta x = x - m_x\f$, \f$\delta y = y - m_y\f$, and
             * \f[C_a = \frac{\cos^2\theta}{2\sigma_x^2}+\frac{\sin^2\theta}{2\sigma_y^2},\quad
             *    C_b = \frac{\sin 2\theta}{4}\!\left(\frac{1}{\sigma_y^2}-\frac{1}{\sigma_x^2}\right),\quad
             *    C_c = \frac{\sin^2\theta}{2\sigma_x^2}+\frac{\cos^2\theta}{2\sigma_y^2}.\f]
             *
             * The integration is semi-analytical: the x-axis integral is computed analytically using
             * the error function, and the resulting 1-D function is integrated numerically along y
             * using the trapezoidal rule with \p nBins_ sample points (minimum 10).
             *
             * @param _x   Centre of the bin along x
             * @param _y   Centre of the bin along y
             * @param A    Amplitude of the Gaussian
             * @param m_x  Mean of the Gaussian along x
             * @param m_y  Mean of the Gaussian along y
             * @param s_x  Standard deviation along x (before rotation)
             * @param s_y  Standard deviation along y (before rotation)
             * @param dx_  Width of the bin along x
             * @param dy_  Width of the bin along y
             * @param _theta  Rotation angle of the Gaussian in radians (default: 0)
             * @param nBins_  Number of sample points for the numerical y-integration (minimum enforced: 10)
             * @return Integral of the Gaussian over the bin \f$[x-dx/2,\,x+dx/2]\times[y-dy/2,\,y+dy/2]\f$
             */
            template<typename U>
            U integrated_gauss2D ( const U& _x, const U& _y, const U& A, const U& m_x, const U& m_y, const U& s_x, const U& s_y, const U& dx_, const U& dy_, U _theta=0, size_t nBins_ = 10)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::normal_distribution] Unsuported datatype or data structure\033[0m");

                size_t nBins = (nBins_ > 10) ? nBins_ : 10;

                const U Ca = 0.5*(std::cos(_theta)*std::cos(_theta)/(s_x*s_x) + std::sin(_theta)*std::sin(_theta)/(s_y*s_y));
                const U Cb = 0.25*(-std::sin(2*_theta)/(s_x*s_x) + std::sin(2*_theta)/(s_y*s_y));
                const U Cc = 0.5*(std::sin(_theta)*std::sin(_theta)/(s_x*s_x) + std::cos(_theta)*std::cos(_theta)/(s_y*s_y));

                const U norm = 1.0 / (2.0 * MathCore::Pi() * s_x * s_y);

                U dbl_n = static_cast<U>(nBins);

                const U x1 = (_x-dx_/2.);
                const U x2 = (_x+dx_/2.);
                const U y1 = (_y-dy_/2.);
                const U y2 = (_y+dy_/2.);

                auto igauss_x = [&x1, &x2, &m_x, &m_y, &Ca, &Cb, &Cc](U y)
                    {
                        return 0.5*std::sqrt(MathCore::Pi()/Ca)*std::exp((y-m_y)*(y-m_y)*(Cb*Cb-Ca*Cc)/Ca)*( std::erf( (Ca*(x2-m_x)-Cb*(y-m_y))/std::sqrt(Ca) ) - std::erf( (Ca*(x1-m_x)-Cb*(y-m_y))/std::sqrt(Ca) ) );
                    };

                U numerical_integral = 0.0;
                for(size_t k=1; k<(nBins-1); ++k)
                {
                    U y0 = y1 + static_cast<U>(k)/(dbl_n-1)*(y2-y1);
                    numerical_integral += igauss_x(y0);
                }

                numerical_integral += igauss_x(y1)/2.;
                numerical_integral += igauss_x(y2)/2.;

                numerical_integral *= (y2-y1) / (dbl_n - 1);

                return A * norm * numerical_integral;
            }

            /**
             * @brief Alias for integrated_gauss2D() — integral of a 2D asymmetric rotated Gaussian over a rectangular bin.
             *
             * @param _x   Centre of the bin along x
             * @param _y   Centre of the bin along y
             * @param A    Amplitude of the Gaussian
             * @param m_x  Mean of the Gaussian along x
             * @param m_y  Mean of the Gaussian along y
             * @param s_x  Standard deviation along x (before rotation)
             * @param s_y  Standard deviation along y (before rotation)
             * @param dx_  Width of the bin along x
             * @param dy_  Width of the bin along y
             * @param _theta  Rotation angle of the Gaussian in radians (default: 0)
             * @param nBins_  Number of sample points for the numerical y-integration (minimum enforced: 10)
             * @return Integral of the Gaussian over the bin \f$[x-dx/2,\,x+dx/2]\times[y-dy/2,\,y+dy/2]\f$
             * @see integrated_gauss2D()
             */
            template<typename U>
            U igauss2D ( const U& _x, const U& _y, const U& A, const U& m_x, const U& m_y, const U& s_x, const U& s_y, const U& dx_, const U& dy_, U _theta=0, size_t nBins_ = 10)
            {
                return integrated_gauss2D(_x, _y, A, m_x, m_y, s_x, s_y, dx_, dy_, _theta, nBins_);
            }

            /**
             * @brief Evaluate the value of a Moffat function at \f$M(x) = \frac{A}\frac{\beta-1}{\pi\alpha^2}\left(1-\frac{(x-x_0)^2}{\alpha^2}\right)^{-\beta}\f$.
             * 
             * @param x_ Position at which to evaluate the Moffat function
             * @param A_ Amplitude of the Moffat function
             * @param x0_ Mean value of the Moffat function
             * @param a_ \f$\alpha\f$ parameter of the Moffat function
             * @param b_ \f$\beta\f$ parameter of the Moffat function
             * @return U 
             */
            template<typename U>
            U moffat ( const U& x_, const U& A_, const U& x0_, const U& a_, const U& b_ )
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::moffat_distribution] Unsuported datatype or data structure\033[0m");

                const U x = x_ - x0_;
                const U norm = (b_-1) / (MathCore::Pi()*a_*a_) ;
                const U moffat_part = std::pow( 1.0 + (x*x)/(a_*a_), -b_);
                return A_ * norm * moffat_part;
            }

            /**
             * @brief Evaluate the value of a 2D Moffat symmetric function at \f$M(x) = \frac{A}\frac{\beta-1}{\pi\alpha^2}\left(1-\frac{(x-x_0)^2+(y-y_0)^2}{\alpha^2}\right)^{-\beta}\f$.
             * 
             * @param x_ Position at which to evaluate the Moffat function
             * @param y_ Position at which to evaluate the Moffat function
             * @param A_ Amplitude of the Moffat function
             * @param x0_ Mean value of the Moffat function
             * @param y0_ Mean value of the Moffat function
             * @param a_ \f$\alpha\f$ parameter of the Moffat function
             * @param b_ \f$\beta\f$ parameter of the Moffat function
             * @return U 
             */
            template<typename U>
            U moffat2D ( const U& x_,const U& y_, const U& A_, const U& x0_, const U& y0_, const U& a_, const U& b_ )
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::normal_distribution] Unsuported datatype or data structure\033[0m");

                const U r = std::sqrt( (x_ - x0_)*(x_ - x0_) + (y_ - y0_)*(y_ - y0_) );
                return moffat<U>( r, A_, 0, a_, b_ );
            }

            /**
             * @brief Evaluate the value of a asymmetric Moffat 2D function \f$M(x) = \frac{A}\left[1+\frac{1}{\alpha^2}\left(x'^2 + \frac{y'^2}{\gamma^2}\right)\right]^{-\beta}\f$ with \f$x'= (x-x_0)\cos(\theta) + (y-y_0)\sin(\theta)\f$ and \f$y'= -(x-x_0)\sin(\theta) + (y-y_0)\cos(\theta)\f$, as defined in Serre et al. Proceedings Volume 7736, Adaptive Optics Systems II; 773649 (2010) https://doi.org/10.1117/12.857277.
             * 
             * @param x_ Position at which to evaluate the Moffat function
             * @param y_ Position at which to evaluate the Moffat function
             * @param A_ Amplitude of the Moffat function
             * @param x0_ Mean value of the Moffat function
             * @param y0_ Mean value of the Moffat function
             * @param a_ \f$\alpha\f$ parameter of the Moffat function
             * @param b_ \f$\beta\f$ parameter of the Moffat function
             * @param g_ \f$\gamma\f$ parameter of the Moffat function
             * @param theta Rotation angle of the Moffat function in radians. The Moffat function is rotated by an angle \f$\theta\f$ with respect to the x-axis.
             * @return U 
             */
            template<typename U>
            U moffat2D ( const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_, const U& a_, const U& b_,const U& g_, U theta=0)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::moffat_distribution] Unsuported datatype or data structure\033[0m");
                if(g_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::moffat_distribution] gamma parameter must be non-zero\033[0m");

                U x1 = (x_-x0_)*std::cos(theta) + (y_-y0_)*std::sin(theta);
                U y1 =-(x_-x0_)*std::sin(theta) + (y_-y0_)*std::cos(theta);

                const U norm = (b_ - 1) / (MathCore::Pi() * a_ * a_ * g_);
                return A_ * norm * std::pow( 1.0 + 1.0/(a_*a_) * (x1*x1 + y1*y1/(g_*g_)) , -b_);
            }

            // ─────────────────────────────────────────────────────────────────
            // Lorentzian (Cauchy) distributions
            // ─────────────────────────────────────────────────────────────────

            /**
             * @brief 1-D Lorentzian (Cauchy) profile.
             *
             * @details Normalised so that \f$\int_{-\infty}^{+\infty} L\,dx = A\f$:
             * \f[
             *   L(x) = A \frac{\gamma/\pi}{(x-x_0)^2 + \gamma^2}
             * \f]
             * Peak value at \f$x_0\f$ is \f$A/(\pi\gamma)\f$.
             *
             * @tparam U   Floating-point type (float or double)
             * @param x_    Evaluation point
             * @param A_    Total flux (integral)
             * @param x0_   Centre
             * @param g_    Half-width at half-maximum \f$\gamma\f$
             * @return      Profile value at \f$x\f$
             */
            template<typename U>
            U lorentzian(const U& x_, const U& A_, const U& x0_, const U& g_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian] Unsupported datatype\033[0m");

                if(g_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian] gamma must be non-zero\033[0m");
                
                    const U dx = x_ - x0_;
                return A_ * (g_ / U(MathCore::Pi())) / (dx*dx + g_*g_);
            }

            /**
             * @brief 2-D symmetric Lorentzian (Cauchy) profile.
             *
             * @details Normalised so that \f$\iint L\,dx\,dy = A\f$:
             * \f[
             *   L(r) = \frac{A\,\gamma}{2\pi\,(r^2+\gamma^2)^{3/2}}
             * \f]
             * where \f$r^2=(x-x_0)^2+(y-y_0)^2\f$.
             *
             * @tparam U   Floating-point type (float or double)
             * @param x_   x coordinate
             * @param y_   y coordinate
             * @param A_   Total flux
             * @param x0_  x centre
             * @param y0_  y centre
             * @param g_   Scale radius \f$\gamma\f$
             * @return     Profile value
             */
            template<typename U>
            U lorentzian2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_, const U& g_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian2D] Unsupported datatype\033[0m");

                if(g_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian2D] gamma must be non-zero\033[0m");

                const U r2 = (x_-x0_)*(x_-x0_) + (y_-y0_)*(y_-y0_);
                return A_ * g_ / (U(2) * U(MathCore::Pi()) * std::pow(r2 + g_*g_, U(1.5)));
            }

            /**
             * @brief 2-D asymmetric (elliptical, rotated) Lorentzian profile.
             *
             * @details Applies a rotation by \f$\theta\f$ and axis-ratio \f$q=b/a\f$:
             * \f[
             *   L = \frac{A\,\gamma}{2\pi\,q\,(r_\mathrm{eff}^2+\gamma^2)^{3/2}},
             *   \quad r_\mathrm{eff}^2 = x_1^2 + y_1^2/q^2
             * \f]
             *
             * @tparam U     Floating-point type
             * @param x_     x coordinate
             * @param y_     y coordinate
             * @param A_     Total flux
             * @param x0_    x centre
             * @param y0_    y centre
             * @param g_     Scale radius \f$\gamma\f$
             * @param q_     Axis ratio \f$b/a \in (0,1]\f$
             * @param theta  Position angle (radians, measured from x-axis)
             * @return       Profile value
             */
            template<typename U>
            U lorentzian2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                           const U& g_, const U& q_, U theta = U(0))
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian2D] Unsupported datatype\033[0m");

                if(g_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian2D] gamma must be non-zero\033[0m");

                if(q_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::lorentzian2D] axis ratio q must be non-zero\033[0m");

                const U x1 =  (x_-x0_)*std::cos(theta) + (y_-y0_)*std::sin(theta);
                const U y1 = -(x_-x0_)*std::sin(theta) + (y_-y0_)*std::cos(theta);
                const U r2 = x1*x1 + y1*y1/(q_*q_);
                return A_ * g_ / (U(2) * U(MathCore::Pi()) * q_ * std::pow(r2 + g_*g_, U(1.5)));
            }

            // ─────────────────────────────────────────────────────────────────
            // Pseudo-Voigt profiles (Thompson-Cox-Hastings approximation)
            // ─────────────────────────────────────────────────────────────────

            /**
             * @brief 1-D pseudo-Voigt profile (Thompson-Cox-Hastings approximation).
             *
             * @details Computes a linear mixture of Gaussian and Lorentzian profiles
             * with a combined FWHM \f$f\f$ and mixing parameter \f$\eta\f$ determined
             * from the Gaussian FWHM \f$f_G\f$ and Lorentzian FWHM \f$f_L\f$:
             * \f[
             *   f^5 = f_G^5 + 2.69269\,f_G^4 f_L + 2.42843\,f_G^3 f_L^2
             *         + 4.47163\,f_G^2 f_L^3 + 0.07842\,f_G f_L^4 + f_L^5
             * \f]
             * \f[\eta = 1.36603\,(f_L/f) - 0.47719\,(f_L/f)^2 + 0.11116\,(f_L/f)^3\f]
             *
             * The total integral equals \f$A\f$.
             *
             * @tparam U   Floating-point type
             * @param x_   Evaluation point
             * @param A_   Total flux (integral)
             * @param x0_  Centre
             * @param fG_  Gaussian FWHM
             * @param fL_  Lorentzian FWHM
             * @return     Profile value at \f$x\f$
             */
            template<typename U>
            U pseudo_voigt(const U& x_, const U& A_, const U& x0_, const U& fG_, const U& fL_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt] Unsupported datatype\033[0m");

                if(fG_ <= U(0) || fL_ <= U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt] FWHM parameters must be positive\033[0m");

                const U f5 = std::pow(fG_, U(5))
                           + U(2.69269) * std::pow(fG_, U(4)) * fL_
                           + U(2.42843) * std::pow(fG_, U(3)) * fL_ * fL_
                           + U(4.47163) * fG_ * fG_ * std::pow(fL_, U(3))
                           + U(0.07842) * fG_ * std::pow(fL_, U(4))
                           + std::pow(fL_, U(5));
                const U f    = std::pow(f5, U(0.2));
                const U rho  = fL_ / f;
                const U eta  = U(1.36603)*rho - U(0.47719)*rho*rho + U(0.11116)*rho*rho*rho;
                const U sig  = f / (U(2) * std::sqrt(U(2) * std::log(U(2))));
                const U gam  = f / U(2);
                const U dx   = x_ - x0_;
                const U G    = std::exp(U(-0.5)*dx*dx/(sig*sig)) / (sig * U(SQT2PI));
                const U L    = (gam / U(MathCore::Pi())) / (dx*dx + gam*gam);
                return A_ * (eta * L + (U(1) - eta) * G);
            }

            /**
             * @brief 2-D symmetric pseudo-Voigt profile.
             *
             * @details Radial extension of the TCH pseudo-Voigt:
             * \f[
             *   V_{2D}(r) = \eta\,L_{2D}(r) + (1-\eta)\,G_{2D}(r)
             * \f]
             * where \f$r^2=(x-x_0)^2+(y-y_0)^2\f$, total integral equals \f$A\f$.
             *
             * @tparam U   Floating-point type
             * @param x_   x coordinate
             * @param y_   y coordinate
             * @param A_   Total flux
             * @param x0_  x centre
             * @param y0_  y centre
             * @param fG_  Gaussian FWHM
             * @param fL_  Lorentzian FWHM
             * @return     Profile value
             */
            template<typename U>
            U pseudo_voigt2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                             const U& fG_, const U& fL_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt2D] Unsupported datatype\033[0m");

                if(fG_ <= U(0) || fL_ <= U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt2D] FWHM parameters must be positive\033[0m");

                const U f5 = std::pow(fG_, U(5))
                           + U(2.69269) * std::pow(fG_, U(4)) * fL_
                           + U(2.42843) * std::pow(fG_, U(3)) * fL_ * fL_
                           + U(4.47163) * fG_ * fG_ * std::pow(fL_, U(3))
                           + U(0.07842) * fG_ * std::pow(fL_, U(4))
                           + std::pow(fL_, U(5));
                const U f    = std::pow(f5, U(0.2));
                const U rho  = fL_ / f;
                const U eta  = U(1.36603)*rho - U(0.47719)*rho*rho + U(0.11116)*rho*rho*rho;
                const U sig  = f / (U(2) * std::sqrt(U(2) * std::log(U(2))));
                const U gam  = f / U(2);
                const U r2   = (x_-x0_)*(x_-x0_) + (y_-y0_)*(y_-y0_);
                const U G2D  = A_ / (U(2) * U(MathCore::Pi()) * sig*sig)
                               * std::exp(U(-0.5) * r2 / (sig*sig));
                const U L2D  = A_ * gam / (U(2) * U(MathCore::Pi()))
                               * std::pow(r2 + gam*gam, U(-1.5));
                return eta * L2D + (U(1) - eta) * G2D;
            }

            /**
             * @brief 2-D asymmetric (elliptical, rotated) pseudo-Voigt profile.
             *
             * @details Same as the symmetric variant but uses the effective radius
             * \f$r_\mathrm{eff}^2 = x_1^2 + y_1^2/q^2\f$ and divides the
             * normalisation by \f$q\f$ to keep the total integral equal to \f$A\f$.
             *
             * @tparam U     Floating-point type
             * @param x_     x coordinate
             * @param y_     y coordinate
             * @param A_     Total flux
             * @param x0_    x centre
             * @param y0_    y centre
             * @param fG_    Gaussian FWHM
             * @param fL_    Lorentzian FWHM
             * @param q_     Axis ratio \f$b/a\f$
             * @param theta  Position angle (radians)
             * @return       Profile value
             */
            template<typename U>
            U pseudo_voigt2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                             const U& fG_, const U& fL_, const U& q_, U theta = U(0))
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt2D] Unsupported datatype\033[0m");

                if(fG_ <= U(0) || fL_ <= U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt2D] FWHM parameters must be positive\033[0m");

                if(q_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::pseudo_voigt2D] axis ratio q must be non-zero\033[0m");

                const U f5 = std::pow(fG_, U(5))
                           + U(2.69269) * std::pow(fG_, U(4)) * fL_
                           + U(2.42843) * std::pow(fG_, U(3)) * fL_ * fL_
                           + U(4.47163) * fG_ * fG_ * std::pow(fL_, U(3))
                           + U(0.07842) * fG_ * std::pow(fL_, U(4))
                           + std::pow(fL_, U(5));
                const U f    = std::pow(f5, U(0.2));
                const U rho  = fL_ / f;
                const U eta  = U(1.36603)*rho - U(0.47719)*rho*rho + U(0.11116)*rho*rho*rho;
                const U sig  = f / (U(2) * std::sqrt(U(2) * std::log(U(2))));
                const U gam  = f / U(2);
                const U x1   =  (x_-x0_)*std::cos(theta) + (y_-y0_)*std::sin(theta);
                const U y1   = -(x_-x0_)*std::sin(theta) + (y_-y0_)*std::cos(theta);
                const U r2   = x1*x1 + y1*y1/(q_*q_);
                const U G2D  = A_ / (U(2) * U(MathCore::Pi()) * sig*sig * q_)
                               * std::exp(U(-0.5) * r2 / (sig*sig));
                const U L2D  = A_ * gam / (U(2) * U(MathCore::Pi()) * q_)
                               * std::pow(r2 + gam*gam, U(-1.5));

                return eta * L2D + (U(1) - eta) * G2D;
            }

            // ─────────────────────────────────────────────────────────────────
            // Sérsic profiles
            // ─────────────────────────────────────────────────────────────────

            /**
             * @brief 1-D (radial) Sérsic profile.
             *
             * @details
             * \f[
             *   I(r) = A \exp\!\left(-b_n\left[\left(\frac{r}{r_e}\right)^{1/n}-1\right]\right)
             * \f]
             * where \f$A = I_e\f$ is the surface brightness at the effective radius \f$r_e\f$
             * and \f$b_n \approx 2n - 1/3 + 4/(405n) + 46/(25515n^2)\f$ (MacArthur 2003).
             *
             * @tparam U   Floating-point type
             * @param r_   Galactocentric radius
             * @param A_   Surface brightness at \f$r_e\f$ (\f$I_e\f$)
             * @param re_  Effective (half-light) radius
             * @param n_   Sérsic index
             * @return     Surface brightness at \f$r\f$
             */
            template<typename U>
            U sersic(const U& r_, const U& A_, const U& re_, const U& n_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::sersic] Unsupported datatype\033[0m");

                if(re_ <= U(0) || n_ <= U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::sersic] re and n must be positive\033[0m");

                const U bn = U(2)*n_ - U(1.0/3.0) + U(4.0)/(U(405)*n_) + U(46.0)/(U(25515)*n_*n_);
                return A_ * std::exp(-bn * (std::pow(r_/re_, U(1)/n_) - U(1)));
            }

            /**
             * @brief 2-D symmetric Sérsic profile.
             *
             * @details Evaluates the Sérsic profile at the projected radius
             * \f$r=\sqrt{(x-x_0)^2+(y-y_0)^2}\f$.  \f$A=I_e\f$ is the
             * surface brightness at \f$r_e\f$.
             *
             * @tparam U   Floating-point type
             * @param x_   x coordinate
             * @param y_   y coordinate
             * @param A_   Surface brightness at \f$r_e\f$
             * @param x0_  x centre
             * @param y0_  y centre
             * @param re_  Effective radius
             * @param n_   Sérsic index
             * @return     Surface brightness
             */
            template<typename U>
            U sersic2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                       const U& re_, const U& n_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::sersic2D] Unsupported datatype\033[0m");

                const U r = std::sqrt((x_-x0_)*(x_-x0_) + (y_-y0_)*(y_-y0_));
                return sersic(r, A_, re_, n_);
            }

            /**
             * @brief 2-D asymmetric (elliptical, rotated) Sérsic profile.
             *
             * @details Uses the effective radius
             * \f$r_\mathrm{eff}=\sqrt{x_1^2 + y_1^2/q^2}\f$.
             * \f$A=I_e\f$ at \f$r_e\f$ along the major axis; total flux scales with \f$q\f$.
             *
             * @tparam U     Floating-point type
             * @param x_     x coordinate
             * @param y_     y coordinate
             * @param A_     Surface brightness at \f$r_e\f$
             * @param x0_    x centre
             * @param y0_    y centre
             * @param re_    Effective radius (along major axis)
             * @param n_     Sérsic index
             * @param q_     Axis ratio \f$b/a\f$
             * @param theta  Position angle (radians)
             * @return       Surface brightness
             */
            template<typename U>
            U sersic2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                       const U& re_, const U& n_, const U& q_, U theta = U(0))
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::sersic2D] Unsupported datatype\033[0m");

                if(q_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::sersic2D] axis ratio q must be non-zero\033[0m");

                const U x1   =  (x_-x0_)*std::cos(theta) + (y_-y0_)*std::sin(theta);
                const U y1   = -(x_-x0_)*std::sin(theta) + (y_-y0_)*std::cos(theta);
                const U reff = std::sqrt(x1*x1 + y1*y1/(q_*q_));
                return sersic(reff, A_, re_, n_);
            }

            // ─────────────────────────────────────────────────────────────────
            // King (1966) profiles
            // ─────────────────────────────────────────────────────────────────

            /**
             * @brief 1-D (radial) King (1966) profile.
             *
             * @details Projected King profile truncated at the tidal radius \f$r_t\f$:
             * \f[
             *   I(r) = \frac{A}{N}\left(\frac{1}{\sqrt{1+(r/r_c)^2}}
             *          - \frac{1}{\sqrt{1+(r_t/r_c)^2}}\right)^2, \quad r \le r_t
             * \f]
             * with normalisation \f$N = \pi r_c^2\left[\ln(1+x_c^2)
             * - 4(1-c_t) + x_c^2 c_t^2\right]\f$, \f$x_c=r_t/r_c\f$,
             * \f$c_t=1/\sqrt{1+x_c^2}\f$, so that \f$\iint I\,dx\,dy = A\f$.
             *
             * @tparam U   Floating-point type
             * @param r_   Projected radius
             * @param A_   Total flux
             * @param rc_  Core radius
             * @param rt_  Tidal radius
             * @return     Surface brightness (0 for \f$r>r_t\f$)
             */
            template<typename U>
            U king(const U& r_, const U& A_, const U& rc_, const U& rt_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::king] Unsupported datatype\033[0m");

                if(rc_ <= U(0) || rt_ <= U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::king] rc and rt must be positive\033[0m");

                if(r_ > rt_) return U(0);
                const U xc   = rt_ / rc_;
                const U ct   = U(1) / std::sqrt(U(1) + xc*xc);
                const U norm = U(MathCore::Pi()) * rc_*rc_
                               * (std::log(U(1) + xc*xc) - U(4)*(U(1)-ct) + xc*xc*ct*ct);
                const U f    = (U(1)/std::sqrt(U(1) + (r_/rc_)*(r_/rc_)) - ct);
                return A_ / norm * f * f;
            }

            /**
             * @brief 2-D symmetric King (1966) profile.
             *
             * @details Evaluates the King profile at
             * \f$r=\sqrt{(x-x_0)^2+(y-y_0)^2}\f$.
             * Total integral equals \f$A\f$.
             *
             * @tparam U   Floating-point type
             * @param x_   x coordinate
             * @param y_   y coordinate
             * @param A_   Total flux
             * @param x0_  x centre
             * @param y0_  y centre
             * @param rc_  Core radius
             * @param rt_  Tidal radius
             * @return     Surface brightness
             */
            template<typename U>
            U king2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                     const U& rc_, const U& rt_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::king2D] Unsupported datatype\033[0m");

                const U r = std::sqrt((x_-x0_)*(x_-x0_) + (y_-y0_)*(y_-y0_));
                return king(r, A_, rc_, rt_);
            }

            /**
             * @brief 2-D asymmetric (elliptical, rotated) King (1966) profile.
             *
             * @details Uses the effective radius
             * \f$r_\mathrm{eff}=\sqrt{x_1^2+y_1^2/q^2}\f$ and adjusts the
             * normalisation by \f$q\f$ so that the total integral remains \f$A\f$.
             *
             * @tparam U     Floating-point type
             * @param x_     x coordinate
             * @param y_     y coordinate
             * @param A_     Total flux
             * @param x0_    x centre
             * @param y0_    y centre
             * @param rc_    Core radius (along major axis)
             * @param rt_    Tidal radius (along major axis)
             * @param q_     Axis ratio \f$b/a\f$
             * @param theta  Position angle (radians)
             * @return       Surface brightness
             */
            template<typename U>
            U king2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                     const U& rc_, const U& rt_, const U& q_, U theta = U(0))
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::king2D] Unsupported datatype\033[0m");

                if(q_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::king2D] axis ratio q must be non-zero\033[0m");

                const U x1   =  (x_-x0_)*std::cos(theta) + (y_-y0_)*std::sin(theta);
                const U y1   = -(x_-x0_)*std::sin(theta) + (y_-y0_)*std::cos(theta);
                const U reff = std::sqrt(x1*x1 + y1*y1/(q_*q_));
                // Compute symmetric King with unit flux then scale by 1/q
                return king(reff, A_ / q_, rc_, rt_);
            }

            // ─────────────────────────────────────────────────────────────────
            // Airy disk profiles
            // ─────────────────────────────────────────────────────────────────

            /**
             * @brief 1-D Airy disk profile.
             *
             * @details
             * \f[
             *   \mathrm{airy}(x) = A\left(\frac{2J_1(a(x-x_0))}{a(x-x_0)}\right)^2
             * \f]
             * where \f$A\f$ is the peak amplitude and \f$a\f$ is the scale parameter
             * (e.g. \f$a=\pi D/\lambda\f$ for diffraction by a circular aperture).
             *
             * @tparam U   Floating-point type
             * @param x_   Evaluation point
             * @param A_   Peak amplitude
             * @param x0_  Centre
             * @param a_   Scale parameter
             * @return     Profile value (equals \f$A\f$ at \f$x=x_0\f$)
             */
            template<typename U>
            U airy(const U& x_, const U& A_, const U& x0_, const U& a_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::airy] Unsupported datatype\033[0m");

                if(a_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::airy] scale parameter a must be non-zero\033[0m");

                const U u = a_ * (x_ - x0_);
                if(u == U(0)) return A_;
                const U j1u = U(MathCore::J1(double(u)));
                return A_ * (U(2) * j1u / u) * (U(2) * j1u / u);
            }

            /**
             * @brief 2-D symmetric Airy disk profile.
             *
             * @details Evaluates the Airy pattern at
             * \f$r=\sqrt{(x-x_0)^2+(y-y_0)^2}\f$.
             * \f$A\f$ is the peak amplitude.
             *
             * @tparam U   Floating-point type
             * @param x_   x coordinate
             * @param y_   y coordinate
             * @param A_   Peak amplitude
             * @param x0_  x centre
             * @param y0_  y centre
             * @param a_   Scale parameter
             * @return     Profile value
             */
            template<typename U>
            U airy2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_, const U& a_)
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::airy2D] Unsupported datatype\033[0m");

                const U r = std::sqrt((x_-x0_)*(x_-x0_) + (y_-y0_)*(y_-y0_));
                return airy(r, A_, U(0), a_);
            }

            /**
             * @brief 2-D asymmetric (elliptical, rotated) Airy disk profile.
             *
             * @details Uses the effective radius
             * \f$r_\mathrm{eff}=\sqrt{x_1^2+y_1^2/q^2}\f$; \f$A\f$ remains
             * the peak amplitude (axis ratio does not change the peak value).
             *
             * @tparam U     Floating-point type
             * @param x_     x coordinate
             * @param y_     y coordinate
             * @param A_     Peak amplitude
             * @param x0_    x centre
             * @param y0_    y centre
             * @param a_     Scale parameter
             * @param q_     Axis ratio \f$b/a\f$
             * @param theta  Position angle (radians)
             * @return       Profile value
             */
            template<typename U>
            U airy2D(const U& x_, const U& y_, const U& A_, const U& x0_, const U& y0_,
                     const U& a_, const U& q_, U theta = U(0))
            {
                if(!is_allowed_math_type<U>)
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::airy2D] Unsupported datatype\033[0m");

                if(q_ == U(0))
                    throw std::invalid_argument("\033[1;35;47m[DST::Math::function::airy2D] axis ratio q must be non-zero\033[0m");
                    
                const U x1   =  (x_-x0_)*std::cos(theta) + (y_-y0_)*std::sin(theta);
                const U y1   = -(x_-x0_)*std::sin(theta) + (y_-y0_)*std::cos(theta);
                const U reff = std::sqrt(x1*x1 + y1*y1/(q_*q_));
                return airy(reff, A_, U(0), a_);
            }
        }

    }
}

#endif
