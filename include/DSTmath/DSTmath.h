//
//  DSTmath.h
//  DeepSkyTools
//
//  File created by GILLARD William on 23/10/15.
//  Mathematical method rewitten based on the ROOT::framwork to be ROOT independent compilent
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#ifndef __DST_MATH_CORE__
#define __DST_MATH_CORE__
#include<cmath>
#include<limits>
#include<vector>

namespace DST
{
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
        };
    
#pragma mark - Func class definition
        /**
         *  @class Func DSTmath.h "DSTmath/DSTmath.h"
         *  @brief Definition of basic math function
         *  @details Func class define common mathematical functions.
         */
        class Func: public MathCore
        {
        public:
            static double Gamma(double, double);
            static double Gamma(double);
    
        };
        
        class polynom : public MathCore
        {
        protected:
            static std::vector<double> chebfit(std::vector<double>, double, double);
            
        public:
            static std::vector<double> chebfit(double (*func)(double), double, double, unsigned int);
            static std::vector<float>  chebfit(float  (*func)(float), float, float, unsigned int);
            
            static double chebev(double, std::vector<double>, double, double, unsigned int n = 0);
            static float  chebev(float, std::vector<float>  , float , float , unsigned int n = 0);
            static void   chebder(std::vector<double>, double, double, std::vector<double>&);
            static void   chebder(std::vector<float>, float, float, std::vector<float>&);
            
            void chebinv(double, double&, std::vector<double>, double, double);

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
        
    }
}

#endif
