
//
//  DSTmath.cxx
//  DeepSkyTools
//
//  File created by GILLARD William on 23/10/15.
//  Mathematical method rewitten based on the ROOT::framwork to be ROOT independent compilent
//
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#include<DSTmath/DSTmath.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include<stdexcept>

//---
/* the machine roundoff error */
#define kMACHEP  1.11022302462515654042363166809e-16

/* largest argument for TMath::Exp() */
#define kMAXLOG  709.782712893383973096206318587

/* smallest argument for TMath::Exp() without underflow */
#define kMINLOG  -708.396418532264078748994506896

/* the maximal number that pow(x,x-0.5) has no overflow */
/* we use a (very) conservative portable bound          */
#define kMAXSTIR  108.116855767857671821730036754

#define kMAXLGM 2.556348e305

#define SQTPI    std::sqrt(DST::Math::MathCore::Pi())     /* sqrt(pi) */

#ifndef SQT2PI
#define SQT2PI   std::sqrt(2.*DST::Math::MathCore::Pi()) /* sqrt(2*pi) */
#endif

namespace DST
{
    /*!
        *  @namespace DST::Math
        *  @brief Deep Sky Tool namespace regroubing classes and function from DeepSky project
    */

    /*!
    @brief Convert floating point to string with a given precision
    @param a_value   floating point to convert
    @param n         number of significant digits
    @return string representation of a_value with n significant digits
    @note This function is provided to be compatible with C++ versions prior to C++11
    */
    std::string to_string(const double& a_value, const int& n)
    {
        std::ostringstream sout;
        sout << std::scientific << std::setprecision(n) << a_value;
        return sout.str();
    }

    /*!
    @brief Convert floating point to string with a given precision
    @param a_value   floating point to convert
    @param n         number of significant digits
    @return string representation of a_value with n significant digits
    @note This function is provided to be compatible with C++ versions prior to C++11
    */
    std::string to_string(const float& a_value, const int& n)
    {
        double a = static_cast<double>( a_value );
        return DST::to_string(a, n);
    }

    namespace Math
    {
        /*!
        * @namespace DST::Math
        * @brief Namespace regrouping mathematical function and classes
        */

        const double MathCore::kBig = 4.503599627370496e15;
        const double MathCore::kBiginv =  2.22044604925031308085e-16;
        
        /*!
        @brief Definition of \f$\log\left( \sqrt{2\times\pi} \right)\f$
         */
        const double MathCore::LS2PI  =  0.91893853320467274178;
        
        /*!
        @brief A[]: log gamma function
         */
        const double MathCore::A[] =
        {
            8.11614167470508450300E-4,
            -5.95061904284301438324E-4,
            7.93650340457716943945E-4,
            -2.77777777730099687205E-3,
            8.33333333333331927722E-2
        };

        /*!
        @brief B[]: log gamma function between 2 and 3
         */
        const double MathCore::B[] =
        {
           -1.37825152569120859100E3,
           -3.88016315134637840924E4,
           -3.31612992738871184744E5,
           -1.16237097492762307383E6,
           -1.72173700820839662146E6,
           -8.53555664245765465627E5
        };

        /*!
        @brief C[]: log gamma function between 2 and 3
         */
        const double MathCore::C[] =
        {
           -3.51815701436523470549E2,
           -1.70642106651881159223E4,
           -2.20528590553854454839E5,
           -1.13933444367982507207E6,
           -2.53252307177582951285E6,
           -2.01889141433532773231E6
        };

        /*!
        @brief B[], C[]: log gamma function between 2 and 3
         */
        const double MathCore::P[] =
        {
            1.60119522476751861407E-4,
            1.19135147006586384913E-3,
            1.04213797561761569935E-2,
            4.76367800457137231464E-2,
            2.07448227648435975150E-1,
            4.94214826801497100753E-1,
            9.99999999999999996796E-1
        };

        /*!
        @brief Q[]: log gamma function between 2 and 3
         */
        const double MathCore::Q[] =
        {
            -2.31581873324120129819E-5,
            5.39605580493303397842E-4 ,
            -4.45641913851797240494E-3,
            1.18139785222060435552E-2 ,
            3.58236398605498653373E-2 ,
            -2.34591795718243348568E-1,
            7.14304917030273074085E-2 ,
            1.00000000000000000320E0
        };

        /**
         *  @brief Stirling's formula for the gamma function
         */
        const double MathCore::STIR[] =
        {
            7.87311395793093628397E-4,
            -2.29549961613378126380E-4,
            -2.68132617805781232825E-3,
            3.47222221605458667310E-3,
            8.33333333333482257126E-2,
        };

        /**
         *  @brief compute PI
         *  @return \f$\pi = \acos(-1)\f$
         */
        /**
         * @brief Bessel function of the first kind, order 1.
         *
         * @details Polynomial approximation from Numerical Recipes (Press et al.).
         * Split at |x| = 8: rational polynomial for small |x|, asymptotic
         * trigonometric expansion for large |x|.
         *
         * @param x  Argument
         * @return   J_1(x)
         */
        double MathCore::J1(double x)
        {
            // Coefficients for |x| < 8
            static const double r[] = {
                 72362614232.0, -7895059235.0,  242396853.1,
                -2972611.439,       15704.48260,   -30.16036606
            };
            static const double s[] = {
                144725228442.0, 2300535178.0, 18583304.74,
                     99447.43394,    376.9991397,       1.0
            };
            // Coefficients for |x| >= 8
            static const double p[] = {
                 1.0,             0.183105e-2,   -0.3516396496e-4,
                 0.2457520174e-5, -0.240337019e-6
            };
            static const double q[] = {
                 0.04687499995, -0.2002690873e-3,  0.8449199096e-5,
                -0.88228987e-6,  0.105787412e-6
            };

            const double ax = std::abs(x);
            if (ax < 8.0)
            {
                const double y   = x * x;
                const double num = x * (r[0] + y*(r[1] + y*(r[2] + y*(r[3] + y*(r[4] + y*r[5])))));
                const double den =       s[0] + y*(s[1] + y*(s[2] + y*(s[3] + y*(s[4] + y*s[5]))));
                return num / den;
            }
            else
            {
                const double z   = 8.0 / ax;
                const double y   = z * z;
                const double xx  = ax - 2.356194491;   // ax - 3π/4
                const double pp  = p[0] + y*(p[1] + y*(p[2] + y*(p[3] + y*p[4])));
                const double qq2 = q[0] + y*(q[1] + y*(q[2] + y*(q[3] + y*q[4])));
                const double ans = std::sqrt(2.0 / (Pi() * ax))
                                 * (std::cos(xx)*pp - z*std::sin(xx)*qq2);
                return (x < 0.0) ? -ans : ans;
            }
        }

        double MathCore::Pi()
        {
            return std::acos(-1);
        }


        /**
         *  @details incomplete gamma function (complement integral)
         *  @param a power index
         *  @param x lower limit of the integral
         *  @return \f{eqnarray*}{
            igamc(a,x)  & = &   1 - igam(a,x) \\
                        & = & \frac{1}{\Gamma(a)}\int^{\infty}_{x}e^{-t}t^{a-1}\,dt\f}
         * @note In this implementation both arguments must be positive. The integral is evaluated by either a power series or continued fraction expansion, depending on the relative values of a and x.
         */
        double MathCore::igamc( double a, double x )
        {
            double ans, ax, c, yc, r, t, y, z;
            double pk, pkm1, pkm2, qk, qkm1, qkm2;

            // LM: for negative values returns 0.0
            // This is correct if a is a negative integer since Gamma(-n) = +/- inf
            if (a <= 0)
                return 0.0;

            if (x <= 0)
                return 1.0;

            if( (x < 1.0) || (x < a) )
                return( 1.0 - igam(a,x) );

            ax = a * std::log(x) - x - lgam(a);

            if( ax < -kMAXLOG )
                return( 0.0 );

            ax = std::exp(ax);

            /* continued fraction */

            y = 1.0 - a;
            z = x + y + 1.0;
            c = 0.0;
            pkm2 = 1.0;
            qkm2 = x;
            pkm1 = x + 1.0;
            qkm1 = z * x;
            ans = pkm1/qkm1;
            t   = 1.;
        
            while( t > kMACHEP )
            {
                c += 1.0;
                y += 1.0;
                z += 2.0;
                yc = y * c;
                pk = pkm1 * z  -  pkm2 * yc;
                qk = qkm1 * z  -  qkm2 * yc;

                if(qk)
                {
                    r = pk/qk;
                    t = std::abs( (ans - r)/r );
                    ans = r;
                }
                else
                    t = 1.0;

                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if( std::abs(pk) > kBig )
                {
                    pkm2 *= kBiginv;
                    pkm1 *= kBiginv;
                    qkm2 *= kBiginv;
                    qkm1 *= kBiginv;
                }
            }

            return( ans * ax );
        }

        /**
         *  @details Left tail of incomplete gamma function:
         *
         *  @param a power index
         *  @param x coordinates
         *
         *  @return \f[ x^a\times e^{-x}\sum^{\infty}_{k=0}\frac{x^k}{\Gamma(a+k+1)}\f]
         */
        double MathCore::igam( double a, double x )
        {
            double ans, ax, c, r;

            // LM: for negative values returns 1.0 instead of zero
            // This is correct if a is a negative integer since Gamma(-n) = +/- inf
            if (a <= 0)
                return 1.0;
            if (x <= 0)
                return 0.0;

            if( (x > 1.0) && (x > a ) )
                return( 1.0 - igamc(a,x) );

            /* Compute  x**a * exp(-x) / gamma(a)  */
            ax = a * std::log(x) - x - lgam(a);
            if( ax < -kMAXLOG )
                return( 0.0 );

            ax = std::exp(ax);

            /* power series */
            r = a;
            c = 1.0;
            ans = 1.0;

            while( c/ans > kMACHEP )
            {
                r += 1.0;
                c *= x/r;
                ans += c;
            }

            return( ans * ax/a );
         }
     
        /**
         *  Logarithm of gamma function
         *  @param x coordinate
         * @return \f$\log(\Gamma(x))\f$
         */
        double MathCore::lgam( double x )
        {
            double p, q, u, w, z;
            int i;
        
            int sgngam = 1;
        
            if (x >= std::numeric_limits<double>::infinity())
                return(std::numeric_limits<double>::infinity());
        
            if( x < -34.0 )
            {
                q = -x;
                w = lgam(q);
                p = std::floor(q);
                if( p==q )//_unur_FP_same(p,q)
                    return (std::numeric_limits<double>::infinity());
                i = (int) p;
                if( (i & 1) == 0 )
                    sgngam = -1;
               else
                   sgngam = 1;
                z = q - p;
                if( z > 0.5 )
                {
                    p += 1.0;
                    z = p - q;
                }
                z = q * std::sin( Math::MathCore::Pi() * z );
                if( z == 0 )
                    return (std::numeric_limits<double>::infinity());

                z = std::log(Math::MathCore::Pi()) - std::log( z ) - w;
                return( z );
            }

            if( x < 13.0 )
            {
                z = 1.0;
                p = 0.0;
                u = x;
                while( u >= 3.0 )
                {
                    p -= 1.0;
                    u = x + p;
                    z *= u;
                }

                while( u < 2.0 )
                {
                    if( u == 0 )
                        return (std::numeric_limits<double>::infinity());
                    z /= u;
                    p += 1.0;
                    u = x + p;
                }

                if( z < 0.0 )
                {
                    sgngam = -1;
                    z = -z;
                }
                else
                    sgngam = 1;

                if( u == 2.0 )
                    return( std::log(z) );

                p -= 2.0;
                x = x + p;
                p = x * Polynomialeval(x, &Math::MathCore::B[0], 5 ) / Polynomial1eval( x, Math::MathCore::C, 6);
                return( std::log(z) + p );
            }
        
            if( x > kMAXLGM )
                return( sgngam * std::numeric_limits<double>::infinity() );
        
            q = ( x - 0.5 ) * std::log(x) - x + LS2PI;
            if( x > 1.0e8 )
                return( q );
        
            p = 1.0/(x*x);
            if( x >= 1000.0 )
            {
                q += (
                      (
                         7.9365079365079365079365e-4 * p
                       - 2.7777777777777777777778e-3
                       )
                      * p
                      + 0.0833333333333333333333
                      )
                    / x;
            }
            else
                q += Polynomialeval( p, &Math::MathCore::A[0], 4 ) / x;
        
            return( q );
        }

        /**
         * @details calculates a value of a polynomial of the form: \f$ a[0]x^N+a[1]x^{(N-1)} + ... + a[N] \f$
         * @param x variable
         * @param a array of coefficients
         * @param N order of the polynomial
         * @return polynomial value \f$y\f$
         */
        double MathCore::Polynomialeval(double x, const double* a, unsigned int N)
        {
            double pom = a[0];

            if (N > 0)
            {
                for (unsigned int i=1; i <= N; i++)
                    pom = pom * x + a[i];
            }

            return pom;
        }

        /**
         * @details calculates a value of a polynomial of the form:
          \f$ x^N+a[0]x^{(N-1)}+ ... + a[N-1] \f$
         */
        double MathCore::Polynomial1eval(double x, const double* a, unsigned int N)
        {
            double pom = a[0];

            if (N > 0)
            {
                pom = x + a[0];
                for (unsigned int i=1; i < N; i++)
                    pom = pom * x + a[i];
            }

            return pom;
        }

        /*!
         * @details Stirling formula for the gamma function
         */
        double MathCore::stirf( double x)
        {
            double y, w, v;

            w = 1.0/x;
            w = 1.0 + w * Polynomialeval( w, &Math::MathCore::STIR[0], 4 );
            y = std::exp(x);

            if( x > kMAXSTIR )
            {
                /* Avoid overflow in pow() */
                v = std::pow( x, 0.5 * x - 0.25 );
                y = v * (v / y);
            }
            else
            {
                y = pow( x, x - 0.5 ) / y;
            }
            y = SQT2PI * y * w;
            return( y );
        }

        double MathCore::small(double z, double x)
        {
            if( x == 0 )
                return( std::numeric_limits<double>::infinity() );

            return( z/((1.0 + 0.5772156649015329 * x) * x) );
        }

        /**
         * @details Computation of the normalized lower incomplete gamma function P(a,x) as defined in the Handbook of Mathematical Functions by Abramowitz and Stegun, formula 6.5.1 on page 260 .
         * @note Its normalization is such that DST::Math::Gamma(a,+infinity) = 1.
         *
         *  @param a power index
         *  @param x lower limit of the
         *
         *  @return \f[ P(a, x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt \f]
         */
        double MathCore::Gamma(double a,double x)
        {
            return DST::Math::MathCore::igam(a, x);
        }

        /**
         *  @details Computation of \f$\Gamma(x)\f$ for all \f$x\f$.
         *  @ref C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
         *
         *  @param x coordinate
         *
         *  @return \f$\Gamma(x)\f$
         */
        double MathCore::Gamma( double x )
        {
            double p, q, z;
            int i;
            int sgngam = 1;

            if (x >=std::numeric_limits<double>::infinity())
                return(x);

            q = std::abs(x);

            if( q > 33.0 )
            {
                if( x < 0.0 )
                {
                    p = std::floor(q);
                    if( p == q )
                    {
                        return( sgngam * std::numeric_limits<double>::infinity());
                    }
                    i = (int) p;
                    if( (i & 1) == 0 )
                        sgngam = -1;

                    z = q - p;

                    if( z > 0.5 )
                    {
                        p += 1.0;
                        z = q - p;
                    }
                    z = q * std::sin( std::acos(-1) * z );
                    if( z == 0 )
                        return( sgngam * std::numeric_limits<double>::infinity());

                    z = std::abs(z);
                    z = std::acos(-1)/(z * Math::MathCore::stirf(q) );
                }
                else
                {
                    z = Math::MathCore::stirf(x);
                }
                return( sgngam * z );
            }

            z = 1.0;
            while( x >= 3.0 )
            {
                x -= 1.0;
                z *= x;
            }

            while( x < 0.0 )
            {
                if( x > -1.E-9 )
                    return DST::Math::MathCore::small(z,x);
                z /= x;
                x += 1.0;
            }

            while( x < 2.0 )
            {
                if( x < 1.e-9 )
                     return DST::Math::MathCore::small(z,x);
                z /= x;
                x += 1.0;
            }

            if( x == 2.0 )
                return(z);
            x -= 2.0;

            p = DST::Math::MathCore::Polynomialeval( x, &DST::Math::MathCore::P[0], 6 );
            q = DST::Math::MathCore::Polynomialeval( x, &DST::Math::MathCore::Q[0], 7 );

            return( z * p / q );
        }


            /**
             *  @details Computes the \f$n\f$ Chebychev coefficients \f$c_{0..n-1}\f$ that approximate the function \f$f\f$ by a sery of Chebychev polynom \f$T_k\f$ :
             \f$
             f(x) ≈ \sum_{k=0}^{n-1}\left[ c_k T_k(y)\right] − c_0/2\, ,
             \f$
             where
             \f$
             y = \frac{x-0.5(b+a)}{0.5(b-a)}
             \f$
             This routine is to be used with moderately large \f$ n\f$ (e.g., 30 or 50), the array of coefficient \f$ c\f$ is subsequently truncated at the smaller value \f$m\f$ such that the value for \f$c_m\f$ and subsequent elements are negligible.
             *
             *  @param func Function value at \f$f(x)\f$
             *  @param a Lower limit of the approximation range
             *  @param b upper limit of the approximation range
             *  @return Array of the Chebychev polynom epxansion coeficient
             *  @cite num_rec_C
             */
            std::vector<double> polynom::chebfit(std::vector<double> f, double a, double b)
            {
                std::vector<double> c;

                double fac=2./static_cast<double>( f.size() );

                for (unsigned int j=0;j< f.size(); j++)
                {
                    double sum=0.0;
                    for (unsigned int k = 0;k < f.size(); k++)
                    {

                        sum += f[k]*cos(Pi() * static_cast<double>( j ) * ( static_cast<double>( k ) + 0.5 )/static_cast<double> ( f.size() ));
                    }

                    c.push_back(fac * sum);
                }

                c[0] /= 2.;

                return c;
            }

            /**
             *  @details Computes the \f$n\f$ Chebychev coefficients \f$c_{0..n-1}\f$ that approximate the function \f$f\f$ by a sery of Chebychev polynom \f$T_k\f$ :
             \f$ 
                f(x) ≈ \sum_{k=0}^{n-1}\left[ c_k T_k(y)\right] − c_0/2\, ,
             \f$
             where
             \f$
             y = \frac{x-0.5(b+a)}{0.5(b-a)}
             \f$
             This routine is to be used with moderately large \f$ n\f$ (e.g., 30 or 50), the array of coefficient \f$ c\f$ is subsequently truncated at the smaller value \f$m\f$ such that the value for \f$c_m\f$ and subsequent elements are negligible.
             *
             *  @param func The function to be approximated by Chebychev polynom
             *  @param a Lower limit of the approximation range
             *  @param b upper limit of the approximation range
             *  @param n order of the polynom expansion
             *  @return Array of the Chebychev polynom epxansion coeficient
             *  @cite num_rec_C
             */
            std::vector<double> polynom::chebfit(double (*func)(double), const double& a, const double& b, const unsigned int& n)
            {
            
                double bpa, bma;

                std::vector<double> f = std::vector<double>(n);

                bma=0.5*(b-a);
                bpa=0.5*(b+a);

                for (unsigned int k=0;k<n;k++)
                {
                    double y=cos(Pi()*( static_cast<double>( k ) + 0.5 )/static_cast<double>( n ));
                    f[k]=(*func)(y*bma+bpa);
                }

                std::vector<double> c = chebfit(f, a, b);

                f.clear();

                return c;
            }

            /**
             *  @details Computes the \f$n\f$ Chebychev coefficients \f$c_{0..n-1}\f$ that approximate the function \f$f\f$ by a sery of Chebychev polynom \f$T_k\f$ :
             \f$
             f(x) ≈ \sum_{k=0}^{n-1}\left[ c_k T_k(y)\right] − c_0/2\, ,
             \f$
             where
             \f$
             y = \frac{x-0.5(b+a)}{0.5(b-a)}
             \f$
             This routine is to be used with moderately large \f$ n\f$ (e.g., 30 or 50), the array of coefficient \f$ c\f$ is subsequently truncated at the smaller value \f$m\f$ such that the value for \f$c_m\f$ and subsequent elements are negligible.
             *
             *  @param func The function to be approximated by Chebychev polynom
             *  @param a Lower limit of the approximation range
             *  @param b upper limit of the approximation range
             *  @param n order of the polynom expansion
             *  @return Array of the Chebychev polynom epxansion coeficient
             *  @cite num_rec_C
             */
            std::vector<float>  polynom::chebfit(float (*func)(float), const float& a, const float& b, const unsigned int& n)
            {

                float bpa, bma;

                std::vector<double> f = std::vector<double>(n);

                bma=0.5*(b-a);
                bpa=0.5*(b+a);

                for (unsigned int k=0;k < n;k++)
                {
                    float y=cos(acos(-1)*( static_cast<float>( k ) + 0.5 )/static_cast<float>( n ));
                    f[k]=static_cast<double>( (*func)(y*bma+bpa) );
                }

                std::vector<double> tmp = chebfit(f, a, b);
                std::vector<float> c;

                for(unsigned int k = 0; k < tmp.size(); k++)
                    c.push_back(static_cast<float>( tmp[k] ) );

                f.clear();

                return c;
            }


            /**
             *  Evaluate the Chebychev polynomial expansion at \f$x\f$ include in the range \f$[a, [b\f$. The Chebyshev polynomial expansion \f$\sum_{k=0}^{n-1} c_kT_k(y) − c_0/2\f$ is evaluated at a point \f$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\f$. Note that the bases convertion from x to y is done automatically. Algorithm is based on Clemshaw    algorithm described in numerical reciepies howver, instead of adding \f$+\times C_0\f$ we add \f$C_0\f$ such that \f$T_0(y)=1\f$.
             *
             *  @param x evaluate chebycgev polynom expansion at \f$x\f$
             *  @param a Lower edge of the range
             *  @param b Upper edge of the range
             *  @param c Chebychev polynom coefficients
             *  @param n Truncate Chebychev polynom coefficients. If \f$n < 0\f$ all chebychev coeficient are used.
             *
             *  @return \f$\sum_{k=0}^{n-1} c_kT_k(y)\f$ for \f$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\f$.   
             *  @cite num_rec_C
             */
            double polynom::chebev(const double& x, const std::vector<double>& c, const double& a, const double& b, unsigned int n)
            {
                double d=0.0,dd=0.0,y,y2;
                if ( (x-a)*(x-b) > 0.0 )
                {
                    throw std::invalid_argument("\033[31m[polynom::chebev]\033[0mx @"+std::to_string(x)+" is out-of-range ["+std::to_string(a)+" , "+std::to_string(b)+"].");
                }

                size_t m = (n > 0) ? ( (n <= c.size() )? n:c.size() ) : c.size();


                y=(2.0 * x-a-b)/(b-a);
                y2 = 2.0*(y);
                // Clenshaw recurrence (rewrote from Numerical recipies and use of std::fma for preciser and faster convergence). Iterate j = m-1 down to 1 inclusive.
                for (size_t j = m; j-- > 1; )
                {
                    const double cj = c[j];
                    const double t  = std::fma(y2, d, cj) - dd; // y2*d - dd + cj
                    dd = d;
                    d  = t;
                }

                return std::fma(y, d, c[0] - dd); // y*d + c[0] - dd
            }

            /**
             *  Evaluate the Chebychev polynomial expansion at \f$x\f$ include in the range \f$[a, [b\f$. The Chebyshev polynomial expansion \f$\sum_{k=0}^{n-1} c_kT_k(y) − c_0/2\f$ is evaluated at a point \f$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\f$.
             *
             *  @param x evaluate chebycgev polynom expansion at \f$x\f$
             *  @param a Lower edge of the range
             *  @param b Upper edge of the range
             *  @param c Chebychev polynom coefficients
             *  @param n Truncate Chebychev polynom coefficients. If \f$n < 0\f$ all chebychev coeficient are used.
             *
             *  @return \f$\sum_{k=0}^{n-1} c_kT_k(y) + c_0\f$ for \f$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\f$
             *  @cite num_rec_C
             */
            float polynom::chebev(const float& x, const std::vector<float>& c, const float& a, const float& b, unsigned int n)
            {
                std::vector<double> tmp;

                for(unsigned int j = 0; j < c.size(); j++)
                    tmp.push_back(static_cast<double>( c[j] ));

                return static_cast<float>( polynom::chebev(static_cast<double>( x ), tmp, static_cast<double>( a ), static_cast<double>( b ), n) );
            }

#if __cplusplus < 199711L
        double polynom::chebev2(double *x, std::vector<double> cx, std::vector<double> cy, double *a, double *b)
        {
            std::vector<unsigned int> order;
            order.push_back(0);
            order.push_back(0);
            
            return chebev2(x, cx, cy, a, b, order);
        }
        
        float polynom::chebev2(float *x, std::vector<float> cx, std::vector<float> cy, float *a, float *b)
        {
            std::vector<unsigned int> order;
            order.push_back(0);
            order.push_back(0);
            
            return chebev2(x, cx, cy, a, b, order);
        }
#endif
      
        /**
         * @brief Evaluate derivative of a Chebychev
         * 
         * @details a new Chebyshev object that approximates the derivative of the existing function over the same range [ \c a, \c b ].
         * @param c Cheapychev polynomial coefficients
         * @param a Lower edge of the range
         * @param b Upper edge of the range
         * @param cder Upon return the Chebychev polynomial coefficients of the derivative
         * 
        */
        void polynom::chebder(std::vector<double> c, double a, double b, std::vector<double>& cder)
        {
            int j;
            double con;
            
            if(c.size() < 1)
                return;
            
            if( cder.size() != c.size() )
                cder.resize(c.size(),0);

            for(size_t k=0; k < cder.size(); k++)
                cder[k] = 0.0;
            
            size_t n = c.size();
            
            cder[n-1] = 0.0;
            
            if(n < 2)
                return;
            
            cder[n-2] = 2*(n-1)*c[n-1];
            
            if(n < 3)
                return;
            
            for(j = static_cast<int>(n-3); j>=0; j--)
                cder[static_cast<size_t>(j)] = cder[static_cast<size_t>(j)+2]+2*(j+1)*c[static_cast<size_t>(j)+1];
            
            con = 2.0/(b-a);
            
            for (size_t k=0; k<cder.size() ;k++)
                cder[k] *= con;
            
            return;
        }
        
        void polynom::chebder(std::vector<float> c, float a, float b, std::vector<float>& cder)
        {
            std::vector<double> c_d;
            std::vector<double> der = std::vector<double>(c.size());
            for(size_t k = 0; k < c.size(); k++)
                c_d.push_back(static_cast<double>(c[k]));
            
            chebder(c_d, static_cast<double>(a), static_cast<double>(b), der);
            
            cder.resize(der.size(), 0);
            
            for(size_t k = 0; k < der.size(); k++)
                cder[k] = static_cast<float>(der[k]);
            
            c_d.clear();
            der.clear();
            
            return;
        }
        
        /**
         * @brief Evaluate the inverse of the Chebychev polynomial expansion at \f$y\f$ include in the range \f$[a, [b\f$. The Chebyshev polynomial expansion \f$\sum_{k=0}^{n-1} c_kT_k(y) − c_0/2\f$ is evaluated at a point \f$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\f$. Note that the bases convertion from x to y is done automatically. Algorithm is based on Newton-Raphson method with backtracking from Numerical Recipes in C.
         * 
         * @param y \f$y\f$ coordinate value where to evaluate the inverse of the Chebychev polynomial expansion
         * @param x Upon return the \f$x\f$ coordinate value that satisfy the Chebychev polynomial expansion \f$y\f$. The initial guess must be provided.
         * @param c The list of Chebychev polynomial coefficients
         * @param a The lower edge of the real space domain
         * @param b The upper edge of the real space domain
         * @param expsillon Precision goal of the methods. Once the residual \f$ y - \sum_{k=0}^{n-1} c_kT_k(y) − c_0/2\f$ is below \c expsillon,  the method stops.
         * @param max_iter Max iteration allowed. If the method does not converge within \c max_iter, an exception is thrown.
         */
        void polynom::chebinv(double y, double& x, std::vector<double> c, double a, double b, double expsillon, size_t max_iter)
        {
            std::vector<double> cder(c.size());
            chebder(c,a,b,cder);

            // Clamp initial guess to [a,b]
            auto clamp = [](double v, double lo, double hi){ return std::max(lo, std::min(v, hi)); };
            x = clamp(x, a, b);

            size_t count=0;
            for (size_t it = 0; it < max_iter; ++it)
            {
                count++;
                double fx  = y - chebev(x, c,    a, b);
                if (std::abs(fx) <= expsillon * std::max(1.0, std::abs(y))) return;

                double dfx =      chebev(x, cder, a, b);
                if (!std::isfinite(dfx) || dfx == 0.0)
                {
                    // Derivative too small: damp toward interval center
                    x = 0.5*(x + 0.5*(a+b));
                    continue;
                }

                double step = fx/dfx;
                double xn   = clamp(x + step, a, b);

                // Backtracking if not improving residual
                double fnew = y - chebev(xn, c, a, b);
                int bt = 0;
                while (std::abs(fnew) > std::abs(fx) && bt < 12)
                {
                    step *= 0.5;
                    xn    = clamp(x + step, a, b);
                    fnew  = y - chebev(xn, c, a, b);
                    ++bt;
                }

                x = xn;
            }

            if(count >= max_iter)
            {
                cder.clear();
                throw std::runtime_error("\033[31m[polynom::chebinv]\033[0m Too many iteration. The algorithm did not converge.");
            }
            cder.clear();
        }
        
        
        /**
         *  Evaluate the 2D Chebychev polynomial expansion at \f$(x,y)\f$ include in the range \f$[(a_x,a_y); (b_x,b_y)[\f$. The 2D Chebyshev polynomial expansion \f$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\f$ is evaluated at a point \f$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\f$ and \f$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\f$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \f$(x,y)\f$. x myst be a 2D array.
         *  @param cij 2D Chebychev coefficient \f$x\f$ axis
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \f$x\f$ and \f$y\f$ axis. If \f$n_i < 0\f$ all chebychev coeficient \f$c_i\f$ along the \f$i\f$ axis are used.
         *
         *  @return \f$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\f$
         */
        double polynom::chebev2(double *x, std::vector<double> cx, std::vector<double> cy, double *a, double *b, std::vector<unsigned int> order)
        {
            
            std::vector<unsigned int> m = {1,1};
            if(order.size() > 0)
                m[0] = (order[0] > 0)? ((order[0] < cx.size() )?order[0]:cx.size()) : 1;

            if(order.size() > 1)
                m[1] = (order[1] > 0)? ((order[1] < cy.size() )?order[1]:cy.size()) : 1;

            double valx =  chebev(x[0],cx,a[0],b[0],m[0]);
            double valy =  chebev(x[1],cy,a[1],b[1],m[1]);
           
            m.clear();
            
            return valx*valy;
        }
        
        /**
         *  Evaluate the 2D Chebychev polynomial expansion at \f$(x,y)\f$ include in the range \f$[(a_x,a_y); (b_x,b_y)[\f$. The 2D Chebyshev polynomial expansion \f$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\f$ is evaluated at a point \f$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\f$ and \f$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\f$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \f$(x,y)\f$. x myst be a 2D array.
         *  @param cij 2D Chebychev coefficient along \f$x\f$ axis
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \f$x\f$ and \f$y\f$ axis. If \f$n_i < 0\f$ all chebychev coeficient \f$c_i\f$ along the \f$i\f$ axis are used.
         *
         *  @return \f$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\f$
         */
        float polynom::chebev2(float *x, std::vector<float> cx, std::vector<float> cy, float *a, float *b, std::vector<unsigned int> n)
        {
            std::vector<double> fcx;
            for(unsigned int i = 0; i < cx.size(); i++)
                fcx.push_back(static_cast<double>( cx[i] ));
            
            std::vector<double> fcy;
            for(unsigned int i = 0; i < cy.size(); i++)
                fcy.push_back(static_cast<double>( cy[i] ));
            
            double *fx = new double[2]; fx[0] = static_cast<double>(x[0]); fx[1] = static_cast<double>(x[1]);
            double *fa = new double[2]; fa[0] = static_cast<double>(a[0]); fa[1] = static_cast<double>(a[1]);
            double *fb = new double[2]; fb[0] = static_cast<double>(b[0]); fb[1] = static_cast<double>(b[1]);
            
            double val = chebev2(fx, fcx, fcy, fa, fb, n);
            
            delete [] fx;
            delete [] fa;
            delete [] fb;
            
            fcx.clear();
            fcy.clear();
            
            return static_cast<float>( val );
        }
        
        /**
         *  Evaluate the 2D Chebychev polynomial expansion at \f$(x,y)\f$ include in the range \f$[(a_x,a_y); (b_x,b_y)[\f$. This 2D Chebyshev polynomial expansion, defined as \f$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} a_{ij}T_i(x')T_j(y')\f$ is evaluated at a point \f$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\f$ and \f$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\f$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \f$(x,y)\f$. x myst be a 2D array.
         *  @param aij Coefficient matrix \f$a_{ij}\f$, flattened <b>row major</b>: the coefficient multiplying \f$T_i(x')T_j(y')\f$ sits at index \f$i\,n_y + j\f$, so the matrix holds \f$n_x\f$ rows of \f$n_y\f$ entries each. Entries beyond \f$n_x n_y\f$ are ignored.
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param order Expansion order along \f$x\f$ and \f$y\f$, i.e. \f$\{n_x, n_y\}\f$. A single entry means \f$n_y = n_x\f$. Both must be non zero.
         *
         *  @return \f$\sum_{i=0}^{n_x-1}\sum_{j=0}^{n_y-1}a_{ij}T_i(x')T_j(y')\f$
         *
         *  @throw std::invalid_argument when order is empty, when either order is zero, or when aij is too small to hold \f$n_x n_y\f$ coefficients.
         */
        double polynom::chebev2(double *x, std::vector<double> aij, double *a, double *b, std::vector<unsigned int> order)
        {
            if(order.size() < 1)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2] Errors ***\033[0m Can't compute chebichev polynome without knowing its dimensions. ["+std::to_string(__LINE__)+std::string("]")).c_str());


            const unsigned int nx = order[0];
            const unsigned int ny = (order.size() >= 2)? order[1] : nx;

            // A zero order used to leave the accumulation loop empty and return 0, which is
            // indistinguishable from a genuinely null expansion.
            if(nx == 0 || ny == 0)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2] Errors ***\033[0m The chebychev expansion order is zero along at least one axis. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            // Widened before multiplying: nx*ny in unsigned int can wrap and let an
            // undersized coefficient matrix through the check.
            if(aij.size() < static_cast<size_t>(nx)*static_cast<size_t>(ny))
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2] Errors ***\033[0m The dimensions of the truncated chebychev coeficient 'a' are insuficient. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            // T_i(x') is needed by every term of the j sum and T_j(y') by every term of the
            // i sum, so each axis is evaluated once per index rather than once per (i,j)
            // pair. chebev() applied to the i-th unit vector returns exactly T_i: it returns
            // y*d + c[0] - dd, carrying no Numerical Recipes style halving of the constant
            // term, so T_0 comes back as 1.
            //
            // Building each unit vector fresh is also what removes the reset bug this
            // replaces. That version carried a single cy across the whole double loop and
            // cleared entry j-1 at step j, so the last entry set by one pass of the inner
            // loop was never cleared when the inner loop restarted: from the second i
            // onwards cy held two non-zero entries and had stopped being a unit vector.
            // Every term past the first row was therefore wrong, for square matrices too.
            std::vector<double> Tx(nx, 0.);
            std::vector<double> Ty(ny, 0.);

            {
                std::vector<double> e(nx, 0.);

                for(unsigned int i = 0; i < nx; i++)
                {
                    e[i]  = 1.;
                    Tx[i] = chebev(x[0], e, a[0], b[0], nx);
                    e[i]  = 0.;
                }
            }

            {
                std::vector<double> e(ny, 0.);

                for(unsigned int j = 0; j < ny; j++)
                {
                    e[j]  = 1.;
                    Ty[j] = chebev(x[1], e, a[1], b[1], ny);
                    e[j]  = 0.;
                }
            }

            double val = 0;

            // Row major: the stride of a row is n_y, its length. It used to be n_x, the
            // number of rows, so the indexing was only self consistent on a square matrix:
            // with n_x > n_y it ran past the end and tripped the bounds check, and with
            // n_x < n_y it silently read the wrong coefficients.
            for(size_t i = 0; i < static_cast<size_t>(nx); i++)
                for(size_t j = 0; j < static_cast<size_t>(ny); j++)
                    val += aij[i*static_cast<size_t>(ny)+j]*Tx[i]*Ty[j];

            return val;

        }
        
        /**
         *  Evaluate the 2D Chebychev polynomial expansion at \f$(x,y)\f$ include in the range \f$[(a_x,a_y); (b_x,b_y)[\f$. This 2D Chebyshev polynomial expansion, defined as \f$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} a_{ij}T_i(x')T_j(y')\f$ is evaluated at a point \f$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\f$ and \f$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\f$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \f$(x,y)\f$. x myst be a 2D array.
         *  @param aij 2D Chebychev coefficient of the polynom order ij
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \f$x\f$ and \f$y\f$ axis. If \f$n_i < 0\f$ all chebychev coeficient \f$c_i\f$ along the \f$i\f$ axis are used.
         *
         *  @return \f$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]}a_{ij}T_i(x)T_j(y)\f$
         */
        float  polynom::chebev2(float  *x, std::vector<float> aij,  float  *a, float  *b, std::vector<unsigned int> order)
        {
            double *d_x = new double[2]; d_x[0] = x[0]; d_x[1] = x[1];
            double *d_a = new double[2]; d_a[0] = a[0]; d_a[1] = a[1];
            double *d_b = new double[2]; d_b[0] = b[0]; d_b[1] = b[1];
            
            std::vector<double> d_aij;
            for(size_t i = 0; i < aij.size(); i++)
                d_aij.push_back(aij[i]);
            
            double val = chebev2(d_x, d_aij, d_a, d_b, order);
            
            delete [] d_x;
            delete [] d_a;
            delete [] d_b;
            d_aij.clear();
            
            return static_cast<float>( val );
        }

        // ----------------------------------------------------------------------------
        // Legacy 2D expansion: the pre-2026 convention of A&A 707, A227 (2026).
        // ----------------------------------------------------------------------------

        /**
         *  Evaluate the pre-2026 2D expansion at \f$(x,y)\f$, the basis in which the NISP
         *  spectroscopic calibration coefficients of A&A 707, A227 (2026),
         *  https://doi.org/10.1051/0004-6361/202555859, are expressed.
         *
         *  @details This is NOT a Chebyshev tensor product. The routine that produced those
         *  coefficients did not reset the unit vector driving the second axis between rows,
         *  so the last entry set by one pass of the inner loop survived into the next. From
         *  the second row onwards the second factor was therefore
         *  \f$T_j(y') + T_{n_y-1}(y')\f$ rather than \f$T_j(y')\f$, except at
         *  \f$j = n_y-1\f$ where it was correct. Row 0 was unaffected, and so were all
         *  expansions with \f$n_x = 1\f$ or \f$n_y = 1\f$.
         *
         *  The loop is reproduced verbatim rather than expressed through
         *  chebev2LegacyToStandard(): an independent implementation is what makes the
         *  equivalence between the two representations testable rather than circular. It is
         *  frozen -- it documents an archived convention and must not be "improved".
         *
         *  The legacy functions are an invertible linear transformation of the Chebyshev
         *  tensor basis, hence span the same space. A surface fitted and evaluated wholly
         *  in this convention is exactly the surface a correct implementation would have
         *  found; only the coefficient representation differs. Published results obtained
         *  this way are correct.
         *
         *  @param x evaluate the expansion at cartesian coordinates \f$(x,y)\f$. x must be a 2D array.
         *  @param aij Legacy coefficients, flattened row major, \f$n_x\f$ rows of \f$n_y\f$.
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param order Expansion order along \f$x\f$ and \f$y\f$, \f$\{n_x, n_y\}\f$.
         *
         *  @return The legacy surface at \f$(x,y)\f$.
         *  @throw std::invalid_argument on the same conditions as chebev2().
         *  @see chebev2LegacyToStandard(), chebev2()
         */
        double polynom::chebev2_legacy(double *x, std::vector<double> aij, double *a, double *b, std::vector<unsigned int> order)
        {
            if(order.size() < 1)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2_legacy] Errors ***\033[0m Can't compute chebichev polynome without knowing its dimensions. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            const unsigned int nx = order[0];
            const unsigned int ny = (order.size() >= 2)? order[1] : nx;

            if(nx == 0 || ny == 0)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2_legacy] Errors ***\033[0m The chebychev expansion order is zero along at least one axis. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            if(aij.size() < static_cast<size_t>(nx)*static_cast<size_t>(ny))
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2_legacy] Errors ***\033[0m The dimensions of the truncated chebychev coeficient 'a' are insuficient. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            std::vector<double> cx(nx, 0.);
            std::vector<double> cy(ny, 0.);

            double val = 0;

            // Verbatim from the original: plain assignment, and cy carried across rows.
            // The row stride is ny, which the original had right -- only the reset was
            // wrong. Do not "fix" the reset here; that is the whole point of this function.
            for(unsigned int i = 0; i < nx; i++)
            {
                if(i > 0)
                    cx[i-1] = 0;
                cx[i] = 1;

                for(unsigned int j = 0; j < ny; j++)
                {
                    const size_t k = static_cast<size_t>(i)*static_cast<size_t>(ny) + j;

                    if(j > 0)
                        cy[j-1] = 0;
                    cy[j] = 1;

                    val += aij[k]*chebev(x[0],cx,a[0],b[0],static_cast<unsigned int>(cx.size()))
                                 *chebev(x[1],cy,a[1],b[1],static_cast<unsigned int>(cy.size()));
                }
            }

            return val;
        }

        /// @copydoc polynom::chebev2_legacy(double*, std::vector<double>, double*, double*, std::vector<unsigned int>)
        float polynom::chebev2_legacy(float *x, std::vector<float> aij, float *a, float *b, std::vector<unsigned int> order)
        {
            double d_x[2] = {static_cast<double>(x[0]), static_cast<double>(x[1])};
            double d_a[2] = {static_cast<double>(a[0]), static_cast<double>(a[1])};
            double d_b[2] = {static_cast<double>(b[0]), static_cast<double>(b[1])};

            std::vector<double> d_aij;
            d_aij.reserve(aij.size());
            for(size_t i = 0; i < aij.size(); i++)
                d_aij.push_back(static_cast<double>(aij[i]));

            return static_cast<float>( chebev2_legacy(d_x, d_aij, d_a, d_b, order) );
        }

        /**
         *  Convert legacy coefficients into standard Chebyshev tensor product coefficients.
         *
         *  @details Expanding the legacy basis in the Chebyshev one leaves every coefficient
         *  untouched except the last of each row after the first, which absorbs the
         *  contamination that row contributed to \f$T_{n_y-1}(y')\f$:
         *
         *  \f[ b_{ij} = a_{ij}, \qquad
         *      b_{i,n_y-1} = a_{i,n_y-1} + \sum_{j=0}^{n_y-2} a_{ij} \quad (i \ge 1) \f]
         *
         *  Row 0 is unchanged, and so is everything when \f$n_x = 1\f$ or \f$n_y = 1\f$,
         *  which is consistent with those shapes never having been affected.
         *
         *  The result satisfies, exactly:
         *  chebev2(x, chebev2LegacyToStandard(aij, order), a, b, order)
         *      == chebev2_legacy(x, aij, a, b, order)
         *
         *  Use this to read the coefficients published in A&A 707, A227 (2026), or any
         *  calibration product that carries no convention keyword, with a standard
         *  evaluator.
         *
         *  @param aij Legacy coefficients, flattened row major, \f$n_x\f$ rows of \f$n_y\f$.
         *  @param order Expansion order, \f$\{n_x, n_y\}\f$.
         *  @return Standard Chebyshev coefficients, same layout and size.
         *  @throw std::invalid_argument on the same conditions as chebev2().
         *  @see chebev2StandardToLegacy(), chebev2_legacy()
         */
        std::vector<double> polynom::chebev2LegacyToStandard(const std::vector<double>& aij, std::vector<unsigned int> order)
        {
            if(order.size() < 1)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2LegacyToStandard] Errors ***\033[0m Can't convert without knowing the dimensions. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            const unsigned int nx = order[0];
            const unsigned int ny = (order.size() >= 2)? order[1] : nx;

            if(nx == 0 || ny == 0)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2LegacyToStandard] Errors ***\033[0m The expansion order is zero along at least one axis. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            if(aij.size() < static_cast<size_t>(nx)*static_cast<size_t>(ny))
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2LegacyToStandard] Errors ***\033[0m The coefficient matrix is too small for the requested order. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            std::vector<double> out(aij);

            for(unsigned int i = 1; i < nx; i++)
            {
                double row = 0.;
                for(unsigned int j = 0; j + 1 < ny; j++)
                    row += aij[static_cast<size_t>(i)*static_cast<size_t>(ny) + j];

                out[static_cast<size_t>(i)*static_cast<size_t>(ny) + (ny-1)] += row;
            }

            return out;
        }

        /**
         *  Convert standard Chebyshev tensor product coefficients into the legacy basis.
         *  @details Exact inverse of chebev2LegacyToStandard(). Use it to express a new fit
         *  in the convention of the published products, or to feed a consumer that still
         *  implements the legacy evaluation.
         *  @param aij Standard Chebyshev coefficients, flattened row major.
         *  @param order Expansion order, \f$\{n_x, n_y\}\f$.
         *  @return Legacy coefficients, same layout and size.
         *  @throw std::invalid_argument on the same conditions as chebev2().
         *  @see chebev2LegacyToStandard()
         */
        std::vector<double> polynom::chebev2StandardToLegacy(const std::vector<double>& aij, std::vector<unsigned int> order)
        {
            if(order.size() < 1)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2StandardToLegacy] Errors ***\033[0m Can't convert without knowing the dimensions. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            const unsigned int nx = order[0];
            const unsigned int ny = (order.size() >= 2)? order[1] : nx;

            if(nx == 0 || ny == 0)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2StandardToLegacy] Errors ***\033[0m The expansion order is zero along at least one axis. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            if(aij.size() < static_cast<size_t>(nx)*static_cast<size_t>(ny))
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2StandardToLegacy] Errors ***\033[0m The coefficient matrix is too small for the requested order. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            std::vector<double> out(aij);

            // The forward map only ever reads columns 0..ny-2, which it does not modify, so
            // subtracting the same row sum inverts it exactly.
            for(unsigned int i = 1; i < nx; i++)
            {
                double row = 0.;
                for(unsigned int j = 0; j + 1 < ny; j++)
                    row += aij[static_cast<size_t>(i)*static_cast<size_t>(ny) + j];

                out[static_cast<size_t>(i)*static_cast<size_t>(ny) + (ny-1)] -= row;
            }

            return out;
        }

        /// @copydoc polynom::chebev2LegacyToStandard(const std::vector<double>&, std::vector<unsigned int>)
        std::vector<float> polynom::chebev2LegacyToStandard(const std::vector<float>& aij, std::vector<unsigned int> order)
        {
            std::vector<double> in;
            in.reserve(aij.size());
            for(size_t i = 0; i < aij.size(); i++)
                in.push_back(static_cast<double>(aij[i]));

            const std::vector<double> out = chebev2LegacyToStandard(in, order);

            std::vector<float> f_out;
            f_out.reserve(out.size());
            for(size_t i = 0; i < out.size(); i++)
                f_out.push_back(static_cast<float>(out[i]));

            return f_out;
        }

        /// @copydoc polynom::chebev2StandardToLegacy(const std::vector<double>&, std::vector<unsigned int>)
        std::vector<float> polynom::chebev2StandardToLegacy(const std::vector<float>& aij, std::vector<unsigned int> order)
        {
            std::vector<double> in;
            in.reserve(aij.size());
            for(size_t i = 0; i < aij.size(); i++)
                in.push_back(static_cast<double>(aij[i]));

            const std::vector<double> out = chebev2StandardToLegacy(in, order);

            std::vector<float> f_out;
            f_out.reserve(out.size());
            for(size_t i = 0; i < out.size(); i++)
                f_out.push_back(static_cast<float>(out[i]));

            return f_out;
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         *
         * @param x coordinate where to evaluate the polynome
         * @param a vector of polynomial coefficients
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        double polynom::polyev(const double& x, const std::vector<double>& a)
        {
            double p = 0;
            int    n = 0;
            for(std::vector<double>::const_iterator ik = a.begin(); ik != a.end(); ++ik)
            {
                p += (*ik)*std::pow(x, n);
                n++;
            }

            return p;
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         * 
         * @param x coordinate where to evaluate the polynome
         * @param a array of polynomial coefficients
         * @param n size of the array of polynomial coefficients
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        double polynom::polyev(const double& x, const double *a, const size_t& n)
        {
            std::vector<double> a_vec(a, a + n);

            return polynom::polyev(x, a_vec);
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         * 
         * @param x coordinate where to evaluate the polynome
         * @param a array of polynomial coefficients
         * @param n size of the array of polynomial coefficients
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        double polynom::polyev(const double& x, const double *a, const unsigned int& n)
        {
            return polynom::polyev(x, a, static_cast<size_t>( n ) );
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         * 
         * @param x coordinate where to evaluate the polynome
         * @param a vector of polynomial coefficients
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        float polynom::polyev(const float& x, const std::vector<float>& a)
        {
            std::vector<double> ad;
            for(size_t k = 0; k < a.size(); k++)
                ad.push_back(static_cast<double>( a[k] ));

            return static_cast<float>( polynom::polyev(x, ad) );
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         * 
         * @param x coordinate where to evaluate the polynome
         * @param a array of polynomial coefficients
         * @param n size of the array of polynomial coefficients
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        float polynom::polyev(const float& x, const float *a, const size_t& n)
        {
            std::vector<float> a_vec(a, a + n);
            return polynom::polyev(x, a_vec);
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         * 
         * @param x coordinate where to evaluate the polynome
         * @param a array of polynomial coefficients
         * @param n size of the array of polynomial coefficients
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        float polynom::polyev(const float& x, const float *a, const unsigned int& n)
        {
            return polynom::polyev(x, a, static_cast<size_t>( n ) );
        }

#if __cplusplus >= 202002L
        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         *
         * @param x coordinate where to evaluate the polynome
         * @param a span of polynomial coefficients (carries its own size)
         * @return double value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        double polynom::polyev(const double& x, std::span<const double> a)
        {
            return polynom::polyev(x, a.data(), a.size());
        }

        /**
         * @brief Evaluate polynome \f$ p(x) = \sum_{k=0}^{n-1} a_k x^k\f$ at \f$x\f$.
         *
         * @param x coordinate where to evaluate the polynome
         * @param a span of polynomial coefficients (carries its own size)
         * @return float value of the polynome \f$p(x)\f$ at \f$x\f$.
         */
        float polynom::polyev(const float& x, std::span<const float> a)
        {
            return polynom::polyev(x, a.data(), a.size());
        }
#endif
    }
}
