//
//  splinefit_utils.h
//  DeepSkyLib
//
//  File created by GILLARD William on 20/07/15.
//  Addaptation and extanded from LxSplineFit function from L. Derome
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#ifndef __DeepSkyLib__splinefit_utils__
#define __DeepSkyLib__splinefit_utils__

#include <stdio.h>
#include <algorithm>
#include <vector>
#include <map>


#ifdef _HAS_ROOT_

#include "TSpline.h"
#include "TGraph.h"
#include "TFitter.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TH1.h"
#include "TDirectory.h"
#include "TMath.h"

#endif

typedef std::vector<double> knots;

namespace DST
{
    namespace Math
    {
        
#pragma mark - spline_functor class definition
        /**
         *  @class spline_functor splinefit_utils.h "DSTmath/splinefit_utils.h"
         *  @brief Spline fuction
         *  @details The DST::math::spline_functor class describe a spline function wich is defined with a finite number of inflection points that are called 'knots'.
         */
        class spline_functor
        {
#pragma mark • Enumerator
            
        public:
            enum spline_mode
            {
                ffull,
                fnorm
            };
            
            enum spline_axes
            {
                LinLin = 0x3,
                LinLog = 0x1,
                LogLin = 0x2,
                LogLog = 0x0
            };
            
            enum spline_extrapolation
            {
                fnormal,
                flin,
                fcontinue,
            };
        
        private:
           
#pragma mark • Private memebers
            spline_functor::spline_axes faxes;
            spline_functor::spline_extrapolation fextrapole;
            spline_functor::spline_mode fmode;
            
            knots fxknots;
            knots fyknots;
            
            double fprecision;
            
            TSpline3 *fspline;
            
            const std::string fopt;
            double fval_start;
            double fval_end;
            
#pragma mark • Private memebers function
            
            static void sort(knots, knots);
        
#pragma mark • Public members
        public:
            
#pragma mark ctor/dtor
            spline_functor();
            spline_functor(const spline_functor&);
            spline_functor(unsigned int nknots,
                           double *xn,
                           spline_axes axes = spline_axes::LinLin,
                           const std::string opt="",
                           double valbeg=0,
                           double valend=0);
            
            spline_functor(unsigned int nknots,
                           double *xn,
                           double *yn,
                           spline_axes axes=spline_axes::LinLin,
                           const std::string opt="",
                           double valbeg=0,
                           double valend=0);
            
            spline_functor(knots xn,
                           spline_axes axes=spline_axes::LinLin,
                           const std::string opt="",
                           double valbeg=0,
                           double valend=0);
            
            spline_functor(knots xn,
                           knots yn,
                           spline_axes axes=spline_axes::LinLin,
                           const std::string opt="",
                           double valbeg=0,
                           double valend=0);
            
            ~spline_functor();
            
#pragma mark Modifiers
            
            void setExtrapolationMethod(spline_functor::spline_extrapolation);
            void setSplineMode         (spline_functor::spline_mode);
            void insertKnots(double xn);
            void insertKnots(double xn, double yn);
    
#pragma mark Operators
            
            double operator() (double *x, double *par);
            
#pragma mark Accessors
            
            inline spline_functor::spline_mode getSplineMode  () const {return fmode;}
            inline spline_functor::spline_axes getAxeMode     () const {return faxes;}
            inline spline_functor::spline_extrapolation getExtrapolationMode() const {return fextrapole;}
            
            inline const knots getXknots() const {return fxknots;}
            inline const knots getYknots() const {return fyknots;}
            
            inline unsigned int getNknots() const {return fxknots.size();}
            
            void print();
        };

        
#ifdef _HAS_ROOT_
        
#pragma mark - spline_fit class definition
        /**
         *  @class spline_fit splinefit_utils.h "DSTmath/splinefit_utils.h"
         *  @brief DST::math::spline_functor wrapper for TRoot
         */
        class spline_fit
        {
        private:
#pragma mark • private member
            
            spline_functor::spline_axes faxes;
            spline_functor::spline_extrapolation fextrapole;
            spline_functor::spline_mode fmode;
            
            knots fxknots;
            knots fyknots;
            
            const std::string fopt;
            double fval_start;
            double fval_end;
        
        public:
            
#pragma mark • ctor/dtor
            
            spline_fit(int nknots,
                        double *xn,
                        spline_functor::spline_axes axes = spline_functor::spline_axes::LinLin,
                        const std::string opt="",
                        double valbeg=0,
                        double valend=0);
            
            spline_fit(knots xn,
                       spline_functor::spline_axes axes = spline_functor::spline_axes::LinLin,
                       const std::string opt="",
                       double valbeg=0,
                       double valend=0);

#pragma mark • Accessors
            TF1 *getFunction(TH1 *h);
            TF1 *getFunction(TH1 &h, double xmin, double  xmax);
            
            TF1 *getFunction(double *y, double xmin, double  xmax);
            TF1 *getFunction(knots y, double xmin, double  xmax);
            TF1 *getFunction(const TGraph &g, double xmin, double xmax);
            TF1 *getFunction(const TGraph &g);
            
            inline unsigned int getNknots() const {return  fxknots.size();}
            inline knots getXknots() const {return  fxknots;}
            inline knots getYknots() const {return  fyknots;}
    
#pragma mark • modifier
            
            void setExtrapolationMethod(spline_functor::spline_extrapolation);
            void setSplineMode(spline_functor::spline_mode);
    
#pragma mark • Static member function
            
            static void getConfidenceIntervals(TObject *obj, Double_t *M=0, Double_t cl = .95);
            static void getConfidenceIntervals(TObject *obj, TMatrixD **MM, Double_t cl = .95);
    
            static TMatrixD* Cov2Cor(const TMatrixD &M);
            
};
        
#endif
        
    }
}





#endif /* defined(__DeepSkyLib__splinefit_utils__) */
