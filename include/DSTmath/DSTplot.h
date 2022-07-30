//
//  DSTplot.h
//  DeepSkyTools
//
//  File created by GILLARD William on 27/05/16.
//  Plotting Tools
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#ifndef __DST_MATH_PLOT__
#define __DST_MATH_PLOT__
#include <vector>
#include <map>
#include <limits>

#include "DSTgeo.h"

using namespace DST::Math;

namespace DST
{
    namespace Plot
    {
        typedef std::vector <point>  data;
        typedef std::vector <point>::iterator  data_iterator;
        typedef std::vector <point>::const_iterator  data_const_iterator;
        
        bool point_sort(const point& a, const point& b);
        bool point_sort_2D(const point& a, const point& b);
        
        
#pragma mark - scatter plot class definition
        class scatter
        {
            
        private:
#pragma mark • private member
            data fscatter;
            std::string fname;
            std::string ftitle;

        public:
#pragma mark • ctor/dtor
            scatter();
            scatter(const data&);
            scatter(size_t, point, ...);
            scatter(size_t, double *, double *);
            scatter(const scatter&);
            
            virtual ~scatter();
            
            scatter operator=(const scatter&);
            
#pragma mark • modifier

#pragma mark 1) Insertion
            void insert(const scatter&);
            
            void push_back(const data&);
            void push_back(const point&);
            void push_up  (const point&);
        
#pragma mark 2) Deletion
            
            void pop_back ();
            void reset();
            void clear();
            
#pragma mark 3) Sorting
            virtual void sort();
            
            
#pragma mark • acessor
            
            inline point operator[](size_t i){return fscatter[i];}
            
            inline size_t size() const {return fscatter.size();}
            inline double x(size_t i) const {return (i < fscatter.size() )?fscatter[i].X():std::numeric_limits<double>::min();}
            inline double y(size_t i) const {return (i < fscatter.size())?fscatter[i].Y():std::numeric_limits<double>::min();}
            
            double maximum(double max_val = -1.*std::numeric_limits<double>::max()) const;
            double minimum(double min_val = std::numeric_limits<double>::max()) const;
            
            double stepSizeAt(size_t) const;
            double MeanStepSize() const;
            
            std::vector<double> getX() const;
            std::vector<double> getY() const;
            
            
#pragma mark • operator
            virtual void operator*=(const scatter&);
            virtual void operator/=(const scatter&);
            virtual void operator+=(const scatter&);
            virtual void operator-=(const scatter&);
            
            virtual void operator*=(double);
            virtual void operator/=(double);
            virtual void operator+=(double);
            virtual void operator-=(double);

#pragma mark • Mathematic estimator
            virtual double eval(double) const;
            
        };
    }
}

#endif