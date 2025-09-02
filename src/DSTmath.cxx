
//
//  DSTmath.cxx
//  DeepSkyTools
//
//  File created by GILLARD William on 23/10/15.
//  Mathematical method rewitten based on the ROOT::framwork to be ROOT independent compilent
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#include<DSTmath/DSTmath.h>
#include <sstream>
#include <iostream>
#include <iomanip>
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

#define SQTPI   std::sqrt(2.*DST::Math::MathCore::Pi())        /* sqrt(2*pi) */

namespace DST
{
    std::string to_string(const double& a_value, const int& n)
    {
        std::ostringstream sout;
        sout << std::scientific << std::setprecision(n) << a_value;
        return sout.str();
    }

    std::string to_string(const float& a_value, const int& n)
    {
        double a = static_cast<double>( a_value );
        return DST::to_string(a, n);
    }

    namespace Math
    {
    const double MathCore::kBig = 4.503599627370496e15;
    const double MathCore::kBiginv =  2.22044604925031308085e-16;
    
    /** @brief log( sqrt( 2*pi ) )
     */
    const double MathCore::LS2PI  =  0.91893853320467274178;
    
    /** @brief Stirling's formula expansion of log gamma
     */
    const double MathCore::A[] =
    {
        8.11614167470508450300E-4,
        -5.95061904284301438324E-4,
        7.93650340457716943945E-4,
        -2.77777777730099687205E-3,
        8.33333333333331927722E-2
    };
    
    /** @brief B[]: log gamma function between 2 and 3
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
    
    /** @brief C[]: log gamma function between 2 and 3
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
    
    /** @brief B[], C[]: log gamma function between 2 and 3
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
     *  @return \f$\pi = acos(-1)\f$
     */
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
     * @details calculates a value of a polynomial of the form:
       \f$ a[0]x^N+a[1]x^{(N-1)} + ... + a[N] \f$
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
    
    /** @details Stirling formula for the gamma function
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
        y = SQTPI * y * w;
        return( y );
    }
    
    double MathCore::small(double z, double x)
    {
        if( x == 0 )
            return( std::numeric_limits<double>::infinity() );
        
        return( z/((1.0 + 0.5772156649015329 * x) * x) );
    }
    
    /**
     *  @details Computation of the normalized lower incomplete gamma function P(a,x) as defined in the Handbook of Mathematical Functions by Abramowitz and Stegun, formula 6.5.1 on page 260 .
     * @note Its normalization is such that DST::Math::Gamma(a,+infinity) = 1.
     *
     *  @param a power index
     *  @param x lower limit of the
     *
     *  @return \f[ P(a, x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt \f]
     */
    double Func::Gamma(double a,double x)
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
    double Func::Gamma( double x )
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
         *  @details Computes the \$n\$ Chebychev coefficients \$c_{0..n-1}\$ that approximate the function \$f\$ by a sery of Chebychev polynom \$T_k\$ :
         \f
         f(x) ≈ \sum_{k=0}^{n-1}\left[ c_k T_k(y)\right] − c_0/2\, ,
         \f
         where
         \f
         y = \frac{x-0.5(b+a)}{0.5(b-a)}
         \f
         This routine is to be used with moderately large \$ n\$ (e.g., 30 or 50), the array of coefficient \$ c\$ is subsequently truncated at the smaller value \$m\$ such that the value for \$c_m\$ and subsequent elements are negligible.
         *
         *  @param func Function value at \$f(x)\$
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
         *  @details Computes the \$n\$ Chebychev coefficients \$c_{0..n-1}\$ that approximate the function \$f\$ by a sery of Chebychev polynom \$T_k\$ :
         \f 
            f(x) ≈ \sum_{k=0}^{n-1}\left[ c_k T_k(y)\right] − c_0/2\, ,
         \f
         where
         \f
         y = \frac{x-0.5(b+a)}{0.5(b-a)}
         \f 
         This routine is to be used with moderately large \$ n\$ (e.g., 30 or 50), the array of coefficient \$ c\$ is subsequently truncated at the smaller value \$m\$ such that the value for \$c_m\$ and subsequent elements are negligible.
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
            
            std::vector<double> f = std::vector<double>(n-1);
            
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
         *  @details Computes the \$n\$ Chebychev coefficients \$c_{0..n-1}\$ that approximate the function \$f\$ by a sery of Chebychev polynom \$T_k\$ :
         \f
         f(x) ≈ \sum_{k=0}^{n-1}\left[ c_k T_k(y)\right] − c_0/2\, ,
         \f
         where
         \f
         y = \frac{x-0.5(b+a)}{0.5(b-a)}
         \f
         This routine is to be used with moderately large \$ n\$ (e.g., 30 or 50), the array of coefficient \$ c\$ is subsequently truncated at the smaller value \$m\$ such that the value for \$c_m\$ and subsequent elements are negligible.
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
            
            std::vector<double> f = std::vector<double>(n-1);
            
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
         *  Evaluate the Chebychev polynomial expansion at \$x\$ include in the range \$[a, [b\$. The Chebyshev polynomial expansion \$\sum_{k=0}^{n-1} c_kT_k(y) − c_0/2\$ is evaluated at a point \$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\$. Note that the bases convertion from x to y is done automatically. Algorithm is based on Clemshaw algorithm described in numerical reciepies howver, instead of adding \$+\times C_0\$ we add \$C_0\$ such that \$T_0(y)=1\$.
         *
         *  @param x evaluate chebycgev polynom expansion at \$x\$
         *  @param a Lower edge of the range
         *  @param b Upper edge of the range
         *  @param c Chebychev polynom coefficients
         *  @param n Truncate Chebychev polynom coefficients. If \$n < 0\$ all chebychev coeficient are used.
         *
         *  @return \$\sum_{k=0}^{n-1} c_kT_k(y)\$ for \$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\$.   
         *  @cite num_rec_C
         */
        double polynom::chebev(const double& x, const std::vector<double>& c, const double& a, const double& b, unsigned int n)
        {
            double d=0.0,dd=0.0,sv,y,y2;
            if ( (x-a)*(x-b) > 0.0 )
            {
                throw std::invalid_argument("\033[31m[polynom::chebev]\033[0mx @"+std::to_string(x)+" is out-of-range ["+std::to_string(a)+" , "+std::to_string(b)+"].");
            }
            
            unsigned int m = (n > 0) ? ( (n <= c.size() )? n:c.size() ) : c.size();
            
            
            y=(2.0 * x-a-b)/(b-a);
            y2 = 2.0*(y);
            for (unsigned int j = m-1; j>=1 ; j--)
            {
                sv=d;
                d=y2*d-dd+c[j];
                dd=sv;
            }
            
            return y*d-dd+c[0];
        }
        
        /**
         *  Evaluate the Chebychev polynomial expansion at \$x\$ include in the range \$[a, [b\$. The Chebyshev polynomial expansion \$\sum_{k=0}^{n-1} c_kT_k(y) − c_0/2\$ is evaluated at a point \$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\$.
         *
         *  @param x evaluate chebycgev polynom expansion at \$x\$
         *  @param a Lower edge of the range
         *  @param b Upper edge of the range
         *  @param c Chebychev polynom coefficients
         *  @param n Truncate Chebychev polynom coefficients. If \$n < 0\$ all chebychev coeficient are used.
         *
         *  @return \$\sum_{k=0}^{n-1} c_kT_k(y) + c_0\$ for \$ y = \frac{x − 0.5(b + a)}{0.5(b − a)}\$
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
        
        void polynom::chebinv(double y, double& x, std::vector<double> c, double a, double b)
        {
            std::vector<double> cder=std::vector<double>(c.size());
            chebder(c,a,b,cder);
            
            double xn    = x;
            size_t count = 0;
            double fx = std::numeric_limits<double>::max();
            double dfx= std::numeric_limits<double>::min();
            double min_x =0;
            double min_fx=std::numeric_limits<double>::max();
            double min_dfx=std::numeric_limits<double>::max();
            
            while( fabs(fx)/fabs(y) > std::numeric_limits<double>::epsilon())
            {
                x = xn;
                count ++;
                
                try
                {
                    fx  = y - chebev(x,c,a,b);
                    dfx =     chebev(x,cder,a,b);

                    min_fx = (fx < min_fx)?fx:min_fx;
                    min_x  = (fx == min_fx)?x:min_x;
                    min_dfx= (dfx < min_dfx)?dfx:min_dfx;
                    
                    xn = x + fx/dfx;
                }
                catch(...)
                {
                    xn = min_x + min_fx/min_dfx/2;
                }

                if(count > 10000)
                {
                    cder.clear();
                    throw std::runtime_error("\033[31m[polynom::chebinv]\033[0m Too many iteration. The algorithm did not converge.");
                }
            }
            cder.clear();
        }
        
        
        /**
         *  Evaluate the 2D Chebychev polynomial expansion at \$(x,y)\$ include in the range \$[(a_x,a_y); (b_x,b_y)[\$. The 2D Chebyshev polynomial expansion \$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\$ is evaluated at a point \$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\$ and \$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \$(x,y)\$. x myst be a 2D array.
         *  @param cij 2D Chebychev coefficient \$x\$ axis
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \$x\$ and \$y\$ axis. If \$n_i < 0\$ all chebychev coeficient \$c_i\$ along the \$i\$ axis are used.
         *
         *  @return \$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\$
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
         *  Evaluate the 2D Chebychev polynomial expansion at \$(x,y)\$ include in the range \$[(a_x,a_y); (b_x,b_y)[\$. The 2D Chebyshev polynomial expansion \$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\$ is evaluated at a point \$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\$ and \$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \$(x,y)\$. x myst be a 2D array.
         *  @param cij 2D Chebychev coefficient along \$x\$ axis
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \$x\$ and \$y\$ axis. If \$n_i < 0\$ all chebychev coeficient \$c_i\$ along the \$i\$ axis are used.
         *
         *  @return \$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]} (c_iT_i(x') − c_x^0/2)\times(c_jT_j(x') − c_j^0/2)\$
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
         *  Evaluate the 2D Chebychev polynomial expansion at \$(x,y)\$ include in the range \$[(a_x,a_y); (b_x,b_y)[\$. This 2D Chebyshev polynomial expansion, defined as \$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} a_{ij}T_i(x')T_j(y')\$ is evaluated at a point \$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\$ and \$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \$(x,y)\$. x myst be a 2D array.
         *  @param aij 2D Chebychev coefficient of the polynom order ij
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \$x\$ and \$y\$ axis. If \$n_i < 0\$ all chebychev coeficient \$c_i\$ along the \$i\$ axis are used.
         *
         *  @return \$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]}a_{ij}T_i(x)T_j(y)\$
         */
        double polynom::chebev2(double *x, std::vector<double> aij, double *a, double *b, std::vector<unsigned int> order)
        {
            if(order.size() < 1)
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2] Errors ***\033[0m Can't compute chebichev polynome without knowing its dimensions. ["+std::to_string(__LINE__)+std::string("]")).c_str());

            
            unsigned int nx = order[0];
            unsigned int ny = (order.size() >= 2)? order[1] : nx;
            
            if(aij.size() < static_cast<size_t>(nx*ny))
                throw std::invalid_argument(std::string("\033[31m[polynom::chebev2] Errors ***\033[0m The dimensions of the truncated chebychev coeficient 'a' are insuficient. ["+std::to_string(__LINE__)+std::string("]")).c_str());
            
            std::vector<double> cx;
            for(unsigned int i = 0; i < nx; i++)
                cx.push_back(0.);

            std::vector<double> cy;
            for(unsigned int j = 0; j < ny; j++)
                cy.push_back(0.);
            
            double val = 0;
            
            for(size_t i = 0; i < cx.size(); i++)
            {
                if(i>0)
                    cx[i-1] *=0;

                cx[i]+=1.;

                for(size_t j = 0; j < cy.size(); j++)
                {
                    size_t k = i*static_cast<size_t>(nx)+j;
                    if(k >= aij.size())
                        throw std::invalid_argument(std::string("\033[31m[polynom::chebev2] Errors ***\033[0m The dimensions of the truncated chebychev coeficient 'a' are insuficient. ["+std::to_string(__LINE__)+std::string("]")).c_str());

                    if(j>0)
                        cy[j-1] *=0;

                    cy[j]+=1.;

                    val += aij[k]*chebev(x[0],cx,a[0],b[0],cx.size())*chebev(x[1],cy,a[1],b[1],cy.size());

                }
            }

            cx.clear();
            cy.clear();
            
            return val;
            
        }
        
        /**
         *  Evaluate the 2D Chebychev polynomial expansion at \$(x,y)\$ include in the range \$[(a_x,a_y); (b_x,b_y)[\$. This 2D Chebyshev polynomial expansion, defined as \$\sum_{i=0}^{n-1}\sum_{j=0}^{m-1} a_{ij}T_i(x')T_j(y')\$ is evaluated at a point \$ x' = \frac{x − 0.5(b_x + a_x)}{0.5(b_x − a_x)}\$ and \$ y' = \frac{y − 0.5(b_y + a_y)}{0.5(b_y − a_y)}\$.
         *
         *  @param x evaluate chebycgev polynom expansion at cartesian coordinates \$(x,y)\$. x myst be a 2D array.
         *  @param aij 2D Chebychev coefficient of the polynom order ij
         *  @param a Lower edge of the range. a must be a 2D array
         *  @param b Upper edge of the range. b must be a 2D array
         *  @param n Truncate Chebychev polynom coefficients along \$x\$ and \$y\$ axis. If \$n_i < 0\$ all chebychev coeficient \$c_i\$ along the \$i\$ axis are used.
         *
         *  @return \$\sum_{i=0}^{n[0]}\sum_{j=0}^{n[1]}a_{ij}T_i(x)T_j(y)\$
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
        
    
    
}
}
