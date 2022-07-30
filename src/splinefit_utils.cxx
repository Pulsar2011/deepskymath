//
//  splinefit_utils.cpp
//  DeepSkyLib
//
//  File created by GILLARD William on 20/07/15.
//  Addaptation and extanded from LxSplineFit function from L. Derome
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#include "splinefit_utils.h"

namespace DST
{
    namespace Math
    {
        
#pragma mark - spline_functor class implementation

#pragma mark • private member function
        
        /**
         *  @brief Sort knots.
         *  @details Sort knots in increasing X coordinates.
         */
        void spline_functor::sort(knots v1, knots v2)
        {
            if(v1.size() < 2)
                return;
            
            while (v2.size() < v1.size())
                v2.push_back(0);
            
            while (v2.size() > v1.size())
                v2.erase(--v2.end());
            
            // sort y according to the order of fxknots
#if __cplusplus >= 199711L
            std::sort(v2.begin(), v2.end(),
                 [&v1](size_t i1, size_t i2) {return v1[i1] < v1[i2];});
            
            // sort x in increasing order
            std::sort(v1.begin(), v1.end());
            
#else
            std::map<double,size_t> index_of_v1;
            unsigned int n = 0;
            for(n = 0; n < v1.size(); n++)
                index_of_v1[v1[n]]=n;
            
            std::vector<double> tmp = std::vector<double>(v1.size());
            n = 0;
            for(std::map<double,size_t>::const_iterator it = index_of_v1.begin(); it != index_of_v1.end(); it++)
            {
                v1[n] = it->first;
                tmp[n]= v2[it->second];
                n++;
            }
            
            for(n = 0; n < tmp.size(); n++)
                v2[n] = tmp[n];
#endif
        }
        
#pragma mark • ctor/dtor
        /**
         *  @brief Default constructor
         */
        spline_functor::spline_functor():fopt("")
        {
            fxknots = knots();
            fyknots = knots();
            
            fprecision = 1e-12;
            fspline = NULL;
            
            fmode = spline_functor::spline_mode::ffull;
            faxes = spline_functor::spline_axes::LinLin;
        }
        
        /**
         *  @brief Copy constructor
         */
        spline_functor::spline_functor(const spline_functor& func): faxes(func.faxes),fmode(func.fmode), fprecision(func.fprecision), fopt(func.fopt),fval_start(func.fval_start), fval_end(func.fval_end)
        {
            fspline = NULL;
            
            if(func.fspline != NULL)
                fspline = static_cast<TSpline3*>(func.fspline->Clone());
            
            for(unsigned int n = 0; n < func.fxknots.size(); n++)
                fxknots.push_back(func.fxknots[n]);
            
            for(unsigned int n = 0; n < func.fyknots.size(); n++)
                    fyknots.push_back(func.fyknots[n]);
        }

        
        /**
         *  @brief Constructor
         *  @details Construct a spline function with finite number of knots positioned at a given x coordinates. The y coordinates of the knots isn't yet defined.
         *
         *  @param nknots       Number of knot
         *  @param xn           X Coordinates of each knots.
         *  @param axes         the nature of the axes scale
         *  @param opt          the option
         *  @param valbeg       Starting x coordinates of the function
         *  @param valend       Ending x coordinates of the function
         *
         *  @note the length of the pointer array xn MUST be equal or larger than nknots.
         */
        spline_functor::spline_functor(unsigned int nknots,
                                       double *xn,
                                       spline_functor::spline_axes axes,
                                       const std::string opt,
                                       double valbeg,
                                       double valend) : faxes(axes), fopt(opt), fval_start(valbeg), fval_end(valend)
        {
            fxknots = knots(nknots);
            
            if( (faxes>>0)&1 )
            {
                for	(unsigned int i = 0; i < fxknots.size() ; i++)
                    fxknots[i] = xn[i];
            }
            else
            {
                for	(unsigned int i = 0; i < fxknots.size() ; i++)
                    fxknots[i] = log(xn[i]);
            }
            
            fyknots = knots(fxknots.size());
            
            fprecision = 1e-12;
            fspline = NULL;
            
            fmode = spline_functor::spline_mode::ffull;
        }
        
        /**
         *  @brief Constructor
         *  @details Construct a spline function with finite number of knots positioned at a given (x,y) coordinates.
         *
         *  @param nknots       Number of knot
         *  @param xn           X Coordinates of each knots.
         *  @param xn           Y Coordinates of each knots.
         *  @param axes         the nature of the axes scale
         *  @param opt          the option
         *  @param valbeg       Starting x coordinates of the function
         *  @param valend       Ending x coordinates of the function
         *
         *  @note the length of the pointer array xn MUST be equal or larger than nknots.
         */
        spline_functor::spline_functor(unsigned int nknots,
                                       double *xn,
                                       double *yn,
                                       spline_functor::spline_axes axes,
                                       const std::string opt,
                                       double valbeg,
                                       double valend) : faxes(axes), fopt(opt), fval_start(valbeg), fval_end(valend)
        {
            fxknots = knots(nknots);
            fyknots = knots(nknots);
            
            if( (faxes>>0)&1 )
            {
                for	(unsigned int i = 0; i < fxknots.size() ; i++)
                    fxknots[i] = xn[i];
            }
            else
            {
                for	(unsigned int i = 0; i < fxknots.size() ; i++)
                    fxknots[i] = log(xn[i]);
            }
            
            if( (faxes>>1)&1 )
            {
                for(unsigned int i = 0; i< fyknots.size() ; i++)
                    fyknots[i] = yn[i] ;
            }
            else
            {
                for(unsigned int i = 0; i< fyknots.size() ; i++)
                    fyknots[i] = log(fabs(yn[i])) ;
            }
            
            fprecision = 1e-12;
            fspline = NULL;
            fmode = spline_functor::spline_mode::fnorm;
        }
        
        /**
         *  @brief Constructor
         *  @details Construct a spline function with finite number of knots positioned at a given x coordinates. The y coordinates of the knots isn't yet defined.
         *  @param nknots       Number of knot
         *  @param xn           X Coordinates of each knots.
         *  @param axes         the nature of the axes scale
         *  @param opt          the option
         *  @param valbeg       Starting x coordinates of the function
         *  @param valend       Ending x coordinates of the function
         *
         *  @note the length of the pointer array xn MUST be equal or larger than nknots.
         */
        spline_functor::spline_functor(knots xn,
                                       spline_axes axes,
                                       const std::string opt,
                                       double valbeg,double valend):faxes(axes), fopt(opt), fval_start(valbeg), fval_end(valbeg)
        {
            if(xn.size() > 0)
            {
                for(unsigned int n = 0; n < xn.size(); n++)
                    fxknots.push_back(xn[n]);
            }
            else
                fxknots = knots();
            
            fyknots = knots(fxknots.size());
            
            fprecision = 1e-12;
            fspline = NULL;
            fmode = spline_functor::spline_mode::fnorm;
        }
        
        /**
         *  @brief Constructor
         *  @details Construct a spline function with finite number of knots positioned at a given x coordinates. The y coordinates of the knots isn't yet defined.
         *  @param xn           X Coordinates of each knots.
         *  @param yn           Y Coordinates of each knots.
         *  @param axes         the nature of the axes scale
         *  @param opt          the option
         *  @param valbeg       Starting x coordinates of the function
         *  @param valend       Ending x coordinates of the function
         *
         *  @note the length of the pointer array xn MUST be equal or larger than nknots.
         */
        spline_functor::spline_functor(knots xn,
                                       knots yn,
                                       spline_axes axes,
                                       const std::string opt,
                                       double valbeg,double valend):faxes(axes),  fopt(opt), fval_start(valbeg), fval_end(valbeg)
        {
            unsigned int nknots = ( xn.size() <= yn.size() )? xn.size() : yn.size();
            
            if(nknots > 0)
            {
                fxknots = knots(nknots);
                fyknots = knots(nknots);
                
                for(unsigned int n = 0; n < nknots; n++)
                {
                    fxknots[n] = xn[n];
                    fyknots[n] = yn[n];
                }
            }
            else
            {
                fxknots = knots();
                fyknots = knots();
            }
            
            fprecision = 1e-12;
            fspline = NULL;
            fmode = spline_functor::spline_mode::fnorm;
        }
        
        /**
         *  @brief Destructor
         */
        spline_functor::~spline_functor()
        {
            if(fspline != NULL)
                delete fspline;
            
            fxknots.clear();
            fyknots.clear();
        }
        
#pragma mark • Modifiers
        
        void spline_functor::setSplineMode(spline_functor::spline_mode _m)
        {
            fmode = _m;
        }
        
        /**
         *  @brief Define the extropolation methods
         *  @details Define either the extrapotation of the spline function, bellow and above the first and last knots, rspectively, is a linear or normal extrapolation.
         *  @param method <#method description#>
         */
        void spline_functor::setExtrapolationMethod(spline_functor::spline_extrapolation method)
        {
            fextrapole = method;
        }
        
        void spline_functor::insertKnots(double xn)
        {
            if( (faxes>>0)&1 )
                fxknots.push_back(xn);
            else
                fxknots.push_back(log(xn));
        }
        
        void spline_functor::insertKnots(double xn, double yn)
        {
            if( (faxes>>0)&1 )
                fxknots.push_back(xn);
            else
                fxknots.push_back(log(xn));
            
            if( (faxes>>1)&1 )
                fyknots.push_back(yn);
            else
                fyknots.push_back(log(yn));
        }
        
#pragma mark • Operators
        
        /**
         *  @brief Eval the spline function at coordinate x[0]
         *
         *  @param x   array coordinates
         *  @param par array that give the y position of each knots
         *  @return <#return value description#>
         */
        double spline_functor::operator() (double *x, double *par)
        {
            double Scale = 1;
            bool Update = false;
            
            sort(fxknots, fyknots);
            
            if(fmode == spline_functor::spline_mode::ffull)
            {
                // Check if spline will have to be updated
                for(unsigned int i = 0; i< fxknots.size() ; i++)
                {
                    if( (faxes>>1)&1 )
                    {
                        if(fabs(fyknots[i] - par[i] ) > fprecision * fabs(par[i]))
                        {
                            Update = true;
                            break;
                        }
                    }
                    else
                    {
                        if(fabs(fyknots[i] - log(fabs(par[i]))) > fprecision * fabs(log(fabs(par[i]))))
                        {
                            Update = true;
                            break;
                        }
                    }
                }
                
                if(fspline == NULL || Update)
                {
                    for(unsigned int i = 0; i< fxknots.size() ; i++)
                        if( (faxes>>1)&1 )
                            fyknots[i] = par[i];
                        else
                            fyknots[i] = log(fabs(par[i])) ;
                }
            }
            else if(fmode == spline_functor::spline_mode::fnorm )
            {
                Scale = par[0];
            }
            
            double val = 0;
            
            // Create or update the spline
            
            if(fspline == NULL || Update)
            {
                if(fspline)
                    delete fspline;
                
                double *_x = &fxknots[0];
                double *_y = &fyknots[0];
                
                fspline = new TSpline3("spline", _x, _y, fxknots.size(), fopt.c_str(),fval_start, fval_end);
            }
            
            
            unsigned int _nknots = fxknots.size();
            
            switch (faxes)
            {
                case spline_functor::spline_axes::LinLin:
                    
                    val = fspline->Eval(x[0]);
                    
                    if(fextrapole == spline_functor::spline_extrapolation::flin)
                    {
                        if(x[0] > fxknots[_nknots-1])
                        {
                            double Dev = fspline->Derivative(fxknots[_nknots-1]);
                            val = (fyknots[_nknots-1] + Dev * ((x[0])-fxknots[_nknots-1]));
                        }
                    
                        if(x[0] < fxknots[0])
                        {
                            double Dev = fspline->Derivative(fxknots[0]);
                            val = (fyknots[0] + Dev * ((x[0])-fxknots[0]));
                        }
                    }
                    else if(fextrapole == spline_functor::spline_extrapolation::fcontinue)
                    {
                        if(x[0] > fxknots[_nknots-1])
                            val = fyknots[_nknots-1];
                    
                        if(x[0] < fxknots[0])
                            val = fyknots[0];
                    }
                    
                    break;
                    
                case spline_functor::spline_axes::LogLin:
                    
                    val = fspline->Eval(log(x[0]));
                    
                    if(fextrapole == spline_functor::spline_extrapolation::flin)
                    {
                        if(log(x[0]) > fxknots[_nknots-1])
                        {
                            double Dev = fspline->Derivative(fxknots[_nknots-1]);
                            val = (fyknots[_nknots-1] + Dev * (log(x[0])-fxknots[_nknots-1]));
                        }
                    
                        if(log(x[0]) < fxknots[0])
                        {
                            double Dev = fspline->Derivative(fxknots[0]);
                            val = ( fyknots[0] + Dev * (log(x[0])-fxknots[0]) );
                        }
                    }
                    else if(fextrapole == spline_functor::spline_extrapolation::fcontinue)
                    {
                        if(log(x[0]) > fxknots[_nknots-1])
                            val = fyknots[_nknots-1];
                    
                        if(log(x[0]) < fxknots[0])
                            val = fyknots[0];
                    }
                
                    break;
            
                case spline_functor::spline_axes::LinLog:
                    
                    val = exp(fspline->Eval(x[0]));
                    
                    if(fextrapole == spline_functor::spline_extrapolation::flin)
                    {
                        if( x[0] > fxknots[_nknots-1] )
                        {
                            double Dev = fspline->Derivative(fxknots[_nknots-1]);
                            val = exp(fyknots[_nknots-1] + Dev*(x[0]-fxknots[_nknots-1]));
                        }
                    
                        if(x[0] < fxknots[0])
                        {
                            double Dev = fspline->Derivative(fxknots[0]);
                            val = exp(fyknots[0] + Dev * (x[0]-fxknots[0]));
                        }
                    }
                    else if(fextrapole == spline_functor::spline_extrapolation::fcontinue)
                    {
                        if(x[0] > fxknots[_nknots-1])
                            val = exp(fyknots[_nknots-1]);
                    
                        if(x[0] < fxknots[0])
                            val = exp(fyknots[0]);
                    }
                    
                    break;
                
                case spline_functor::spline_axes::LogLog:
                    
                    val = exp(fspline->Eval(log(x[0])));
                    
                    if(fextrapole == spline_functor::spline_extrapolation::flin)
                    {
                        if( (log(x[0]) ) > fxknots[_nknots-1])
                        {
                            double Dev = fspline->Derivative(fxknots[_nknots-1]);
                            val = exp( fyknots[_nknots-1] + Dev * ( log(x[0]) - fxknots[_nknots-1] ));
                        }
                    
                        if(log( x[0] ) < fxknots[0])
                        {
                            double Dev = fspline->Derivative(fxknots[0]);
                            val = exp(fyknots[0] + Dev * ( log(x[0]) - fxknots[0] ));
                        }
                    }
                    else if(fextrapole == spline_functor::spline_extrapolation::fcontinue)
                    {
                        if( log(x[0]) > fxknots[_nknots-1])
                            val = exp( fyknots[_nknots-1] );
                    
                        if( log(x[0]) < fxknots[0] )
                            val = exp(fyknots[0]);
                    }
                    
                    break;
            }
            
            //std::cout<<Scale<<" * "<<val<<std::endl;
            val *= Scale;
            
            return val;
        }
        
        void spline_functor::print()
        {
            std::cout<<" • SPLINE FUNCTION"<<std::endl
                     <<"   |- axes    :"<<faxes<<std::endl
                     <<"   |- mode    :"<<fmode<<std::endl
                     <<"   |- ext.    :"<<fextrapole<<std::endl
                     <<"   |- N knots : "<<fxknots.size()<<std::endl;
            
            unsigned int n = 0;
            while( n < fxknots.size()-1)
            {
                if( (faxes>>0)&1 )
                    std::cout<<"      |- "<<fxknots[n]<<"    ";
                else
                    std::cout<<"      |- "<<exp(fxknots[n])<<"    ";
                
                if( (faxes>>1)&1 )
                    std::cout<<fyknots[n]<<std::endl;
                else
                    std::cout<<exp(fyknots[n])<<std::endl;
                
                n++;
            }
            
            if( (faxes>>0)&1 )
                std::cout<<"      `- "<<fxknots[n]<<"    ";
            else
                std::cout<<"      `- "<<exp(fxknots[n])<<"    ";
            
            if( (faxes>>1)&1 )
                std::cout<<fyknots[n]<<std::endl;
            else
                std::cout<<exp(fyknots[n])<<std::endl;
        }


#ifdef _HAS_ROOT_
        
#pragma mark - spline_fit class implementation
        
#pragma mark • ctor/dtor
        
        /**
         *  @brief Constructor
         *  @param nknots number of knots
         *  @param xn     x coordinates of each knots
         *  @param axes   Type of axes scale
         *  @param opt    User's option
         *  @param valbeg Starting point of the x axis
         *  @param valend ending point of the x axes
         */
        spline_fit::spline_fit(int nknots,
                               double *xn,
                               spline_functor::spline_axes axes,
                               const std::string opt,
                               double valbeg,
                               double valend) : faxes(axes),fopt(opt),fval_start(valbeg), fval_end(valend)
        {
            fxknots = knots(nknots);
            
            for	(unsigned int i = 0; i < fxknots.size() ; i++)
            {
                if( (faxes>>0)&1 )
                    fxknots[i] = xn[i];
                else
                    fxknots[i] = log( xn[i] );
            }
            
            fyknots    = knots(fxknots.size());
            fextrapole = spline_functor::spline_extrapolation::fnormal;
        }
        
        /**
         *  @brief Constructor
         *  @param xn     x coordinates of each knots
         *  @param axes   Type of axes scale
         *  @param opt    User's option
         *  @param valbeg Starting point of the x axis
         *  @param valend ending point of the x axes
         */
        spline_fit::spline_fit(knots xn,
                               spline_functor::spline_axes axes,
                               const std::string opt,
                               double valbeg,
                               double valend) : faxes(axes),fopt(opt),fval_start(valbeg), fval_end(valend)
        {
            fxknots = knots( xn.size() );
            
            for	(unsigned int i = 0; i < fxknots.size() ; i++)
            {
                if( (faxes>>0)&1 )
                    fxknots[i] = xn[i];
                else
                    fxknots[i] = log( xn[i] );
            }
            
            fyknots    = knots(fxknots.size());
            fextrapole = spline_functor::spline_extrapolation::fnormal;
        }
        
#pragma mark • Accessors
        
        /**
         *  @brief Generate DST::maht::spline_functor based on an inpute histogram  TRoot::TH1
         *  @param h    Input histogram
         *  @return the generated DST::Math::spline_functor embeded into a TRoot::TF1 objects
         */
        TF1 *spline_fit::getFunction(TH1 *h)
        {
            if(h == NULL)
                return NULL;
            
            return getFunction(*h, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        }
        
        /**
         *  @brief Generate DST::maht::spline_functor based on an inpute histogram  TRoot::TH1
         *  @param h    Input histogram
         *  @param xmin lower edge of the x axis
         *  @param xmax upper edge of the x axis
         *  @return the generated DST::Math::spline_functor embeded into a TRoot::TF1 objects
         */
        TF1 *spline_fit::getFunction(TH1 &h, double xmin,double  xmax)
        {
            for(unsigned int i = 0 ; i < fxknots.size() ; i++)
            {
                double _x = 0;
                if( (faxes>>0)&1 )
                    _x = fxknots[i];
                else
                    _x = exp(fxknots[i]);
                
                fyknots[i] = h.Interpolate(_x);
            }
            
            spline_functor *s= new spline_functor(fxknots, faxes, fopt, fval_start, fval_end);
            
            s->setExtrapolationMethod(fextrapole);
            s->setSplineMode(fmode);
            
            //s->print();
            
            TF1* f = new TF1(Form("f_%s",h.GetName()), s, xmin, xmax, fxknots.size(),"spline_fuction");
            
            f->SetNpx(100000);
            f->SetParameters(&fyknots[0]);
            
            return f;
        }
        
        
        /**
         *  @brief Generate DST::maht::spline_functor based on an inpute array
         *  @param y    Input array
         *  @param xmin lower edge of the x axis
         *  @param xmax upper edge of the x axis
         *
         *  @note The array size MUST be equal to the number of knots (i.e. the size of the fxknots vector)
         *  @return the generated DST::Math::spline_functor embeded into a TRoot::TF1 objects
         */
        TF1 *spline_fit::getFunction(double *y,double xmin,double  xmax)
        {
            for(unsigned int n = 0; n < fxknots.size(); n++)
                fyknots[n] = y[n];
            
            return getFunction(fyknots, xmin, xmax);
        }
        
        /**
         *  @brief Generate DST::maht::spline_functor based on an inpute array
         *  @param y    Input array
         *  @param xmin lower edge of the x axis
         *  @param xmax upper edge of the x axis
         *
         *  @return the generated DST::Math::spline_functor embeded into a TRoot::TF1 objects
         */
        TF1 *spline_fit::getFunction(knots y,double xmin,double  xmax)
        {
            for(unsigned int n = 0; n < y.size(); n++)
                fyknots[n] = y[n];
            
            spline_functor *s= new spline_functor(fxknots, faxes, fopt, fval_start, fval_end);
            
            s->setExtrapolationMethod(fextrapole);
            s->setSplineMode(fmode);
            
            //s->print();
            
            TF1* f = new TF1("f", s, xmin, xmax, fxknots.size(),"spline_function");
            f->SetNpx(10000);
            f->SetParameters(&fyknots[0]);
            
            return f;
        }
        
        /**
         *  @brief Generate DST::maht::spline_functor based on an input TGraph
         *  @param g    Input TGraph
         *  @return the generated DST::Math::spline_functor embeded into a TRoot::TF1 objects
         */
        TF1 *spline_fit::getFunction(const TGraph &g)
        {
            return  getFunction(g, g.GetX()[0], g.GetX()[g.GetN()-1]);
        }

        /**
         *  @brief Generate DST::maht::spline_functor based on an input TGraph
         *  @param g    Input TGraph
         *  @param xmin lower edge of the x axis
         *  @param xmax upper edge of the x axis
         *  @return the generated DST::Math::spline_functor embeded into a TRoot::TF1 objects
         */
        TF1 *spline_fit::getFunction(const TGraph &g, double xmin, double xmax)
        {
            for(unsigned int i = 0 ; i < fxknots.size() ; i++)
            {
                double _x = 0;
                if( (faxes>>0)&1 )
                    _x = fxknots[i];
                else
                    _x = exp(fxknots[i]);
                
                fyknots[i] = g.Eval(_x);
            }
            
            std::string fName =  std::string("f_") + std::string( g.GetName() );
            
            spline_functor *s= new spline_functor(fxknots, faxes, fopt, fval_start, fval_end);
            
            s->setExtrapolationMethod(fextrapole);
            s->setSplineMode(fmode);
            
            //s->print();
            
            TF1* f = new TF1(fName.c_str(), s, xmin, xmax, fxknots.size(),"spline_function");
            
            f->SetNpx(100000);
            f->SetParameters(&fyknots[0]);
            
            return f;
        }
        
#pragma mark • modifier
        
        /**
         *  Define which extrapolation method has to be used bellow and above the first and last knots, respectively.
         *  @param meth The extrapolation methods
         */
        void spline_fit::setExtrapolationMethod(spline_functor::spline_extrapolation meth)
        {
            fextrapole = meth;
        }
        
        void spline_fit::setSplineMode(spline_functor::spline_mode _m)
        {
            fmode = _m;
        }
        
#pragma mark • Static member function
        
        void spline_fit::getConfidenceIntervals(TObject *obj, Double_t *M, Double_t cl)
        {
            // Adapted from TFitter::GetConfidenceIntervals
            TVirtualFitter * fit = (TVirtualFitter::GetFitter());
            TH1 *hfit = (TH1*)obj;
            TF1 *f = (TF1*)fit->GetUserFunc();
            
            Int_t npar = f->GetNpar();
            Double_t *grad1 = new Double_t[npar];
            Double_t *grad2 = new Double_t[npar];
            Double_t *sum_vector = new Double_t[npar];
            Double_t x1[1];
            Double_t x2[1];
            
            Int_t hxfirst = 1;
            Int_t hxlast  = hfit->GetNbinsX();
            Int_t N  = hxlast - hxfirst +1 ;
            
            TAxis *xaxis  = hfit->GetXaxis();
            Double_t t = TMath::StudentQuantile(0.5 + cl/2, f->GetNDF());
            Double_t chidf = TMath::Sqrt(f->GetChisquare()/f->GetNDF());
            std::cout << " hxfirst  = " << hxfirst << " hxlast = " << hxlast << std::endl
            << " chidf  = " << chidf << " t = " << t << std::endl;
            
            Double_t *matr = fit->GetCovarianceMatrix();
            Double_t c=0;
            
            for (Int_t binx1=hxfirst; binx1<=hxlast; binx1++)
            {
                x1[0]=xaxis->GetBinCenterLog(binx1);
                f->GradientPar(x1, grad1);
                
                for (Int_t binx2=binx1; binx2<=hxlast; binx2++)
                {
                    x2[0]=xaxis->GetBinCenterLog(binx2);
                    f->GradientPar(x2, grad2);
                    for (Int_t irow=0; irow<npar; irow++)
                    {
                        sum_vector[irow]=0;
                        for (Int_t icol=0; icol<npar; icol++)
                            sum_vector[irow]+=matr[irow*npar+icol]*grad1[icol];
                    }
                    c = 0;
                    
                    for (Int_t i=0; i<npar; i++) c+=grad2[i]*sum_vector[i];
                    
                    if(binx1==binx2)
                    {
                        hfit->SetBinContent(binx1, f->EvalPar(x1));
                        hfit->SetBinError(binx1 , sqrt(c)*t*chidf);
                    }
                    
                    if(M==0)
                        continue;
                    
                    M[(binx2-1)*N+binx1-1] = c*t*t*chidf*chidf;
                    
                    if(binx1!=binx2)
                        M[(binx1-1)*N+binx2-1] = M[(binx2-1)*N+binx1-1] ;
                }
            }
            
            delete [] grad1;
            delete [] grad2;
            delete [] sum_vector;
        }
        
        void spline_fit::getConfidenceIntervals(TObject *obj, TMatrixD **MM, Double_t cl)
        {
            TH1 *h = (TH1*)obj;
            Double_t *M = new Double_t [h->GetNbinsX()*h->GetNbinsX() ];
            
            getConfidenceIntervals(h,M,cl);
            
            (*MM) = new TMatrixD(h->GetNbinsX(),h->GetNbinsX(),M);
            
            delete [] M;
        }
        
        TMatrixD* spline_fit::Cov2Cor(const TMatrixD &M)
        {
            TMatrixD* C = (TMatrixD*) M.Clone();
            
            for (Int_t binx1=0; binx1 < M.GetNrows(); binx1++)
            {
                for (Int_t binx2=0; binx2 < M.GetNrows(); binx2++)
                {
                    (*C)[binx1][binx2] = M[binx1][binx2]
                                       / sqrt( M[binx1][binx1]* M[binx2][binx2] );
                }
            }
            return C;
        }
#endif
        
    }
}