//
//  DSTplot.cxx
//  DeepSkyTools
//
//  File created by GILLARD William on 27/05/16.
//  Plotting tools
//
//  Centre de Physic des Particules de Marseille
//	Copyright (c) 2015, All rights reserved
//

#include "DSTplot.h"
#include <algorithm> 

namespace DST
{
    namespace Plot
    {
        bool point_sort(const point& a, const point& b)
        {
            return a.X() <= b.X();
        }
        
        bool point_sort_2D(const point& a, const point& b)
        {
            return sqrt(a.X()*a.X()+a.Y()+a.Y()) <= sqrt(b.X()*b.X()+b.Y()+b.Y());
        }
        
#pragma mark - scatter plot class definition
#pragma mark • ctor/dtor
        /**
         *  @brief default constructor
         */
        scatter::scatter()
        {
            fscatter = data();
            fname    = std::string();
            ftitle   = std::string();
        }
        
        /**
         *  @brief Construct from list of scattered point
         *  @param in Input scattered data point
         */
        scatter::scatter(const data& in)
        {
            fname    = std::string();
            ftitle   = std::string();
            
            for(data_const_iterator it = in.begin(); it != in.end(); it++)
                fscatter.push_back(point(*it));
        }
        
        scatter::scatter(size_t N, double *x, double *y)
        {
            fscatter = data();
            fname    = std::string();
            ftitle   = std::string();
            
            for(size_t i = 0; i < N; i++)
                fscatter.push_back(point(x[i],y[i]));
        }

        /**
         *  @brief Construct from list of scattered point
         *  @param in Input scattered data point
         */
        scatter::scatter(const scatter& in):fname(in.fname),ftitle(in.ftitle)
        {
            
            for(data_const_iterator it = (in.fscatter).begin(); it != (in.fscatter).end(); it++)
                fscatter.push_back(point(*it));
        }
        
        /**
         *  @brief destructor
         */
        scatter::~scatter()
        {
            fscatter.clear();
            fname.clear();
            ftitle.clear();
        }
        
        std::vector<double> scatter::getX() const
        {
            std::vector<double> x;
            for(size_t i = 0; i < fscatter.size(); i++)
                x.push_back(fscatter[i].X());
            
            return x;
        }
        
        std::vector<double> scatter::getY() const
        {
            std::vector<double> y;
            for(size_t i = 0; i < fscatter.size(); i++)
                y.push_back(fscatter[i].Y());
                
                return y;
        }
        
        
#pragma mark • modifier
#pragma mark 1) Insertion
        
        /**
         *  @brief Inserts scattered points
         *  @param in point to be added.
         */
        void scatter::insert(const scatter& in)
        {
            push_back(in.fscatter);
        }
        
        /**
         *  @brief Inserts scattered points
         *  @param in point to be added.
         */
        void scatter::push_back(const data& in)
        {
            for(data_const_iterator it = in.begin(); it != in.end(); it++)
                push_back(*it);
        }
        
        /**
         *  @brief Inserte a new point at the end
         *  @param in point to be added.
         */
        void scatter::push_back(const point& in)
        {
            fscatter.push_back(point(in));
        }
        
        /**
         *  @brief Inserte a new point  at the begining
         *  @param in point to be added.
         */
        void scatter::push_up  (const point& in)
        {
            fscatter.insert(fscatter.begin(),point(in));
        }
        
#pragma mark 2) Deletion
        /**
         *  @brief Delete last element
         */
        void scatter::pop_back()
        {
            fscatter.pop_back();
        }
        
        /**
         *  @brief clear scatter point
         */
        void scatter::reset()
        {
            clear();
        }
        
        /**
         *  @brief clear scatter point
         */
        void scatter::clear()
        {
            while (fscatter.size() > 1)
            {
                pop_back();
            }
            
            fscatter.clear();
        }

#pragma mark 3) Sorting
        void scatter::sort()
        {
            std::sort(fscatter.begin(), fscatter.end(),point_sort);
        }
        
#pragma mark • acessor
        
        /**
         *  @details Get the maximum value of the scatter point within range
         *  @param max_val (optional) Threshold bellow which maximum shall be searched
         */
        double scatter::maximum(double max_val) const
        {
            double max = -1.*std::numeric_limits<double>::max();
            for(data_const_iterator it = fscatter.begin(); it != fscatter.end(); it++)
            {
                if((*it).Y() >= max_val)
                    continue;
                
                if((*it).Y() > max)
                    max = (*it).Y();
            }
            
            return max;
        }
        
        /**
         *  @details Get the minimum value of the scatter point within range
         *  @param min_val (optional) Threshold above which minimum shall be searched
         */
        double scatter::minimum(double min_val) const
        {
            double min = std::numeric_limits<double>::max();
            
            for(data_const_iterator it = fscatter.begin(); it != fscatter.end(); it++)
            {
                if((*it).Y() <= min_val)
                    continue;
                
                if((*it).Y() < min)
                    min = (*it).Y();
            }
            
            return min;
        }
        
        double scatter::stepSizeAt(size_t i) const
        {
            if(i > fscatter.size()-1)
                return 0;
            
            return (fscatter[i+1].X()-fscatter[i].X());
        }
        
        double scatter::MeanStepSize() const
        {
            double mean = 0;
            
            for(size_t i = 1; i < fscatter.size(); i++)
            {
                mean += fscatter[i+1].X()-fscatter[i].X();
            }
            
            mean /= static_cast<double>( fscatter.size() -1 );
            return mean;
        }
        
#pragma mark • operator
        void scatter::operator*=(const scatter& in)
        {
            point p = point(fscatter[0].size());
            for(size_t i =0; i< p.size(); i++)
                p.SetPoint(i,1.);
            
            for(data_iterator it = fscatter.begin(); it != fscatter.end(); it++)
            {
                p.SetPoint(1, (in.eval((*it).X()) != 0)? in.eval((*it).X()):1 );
                (*it) *= p;
            }
        }
        
        void scatter::operator/=(const scatter& in)
        {
            point p = point(fscatter[0].size());
            for(size_t i =0; i< p.size(); i++)
                p.SetPoint(i,1.);
            
            for(data_iterator it = fscatter.begin(); it != fscatter.end(); it++)
            {
                p.SetPoint(1, (in.eval((*it).X()) != 0)? 1./in.eval((*it).X()):1 );
                (*it) *= p;
            }
        }
        
        void scatter::operator+=(const scatter& in)
        {
            point p = point(fscatter[0].size());
            for(size_t i =0; i< p.size(); i++)
                p.SetPoint(i,1.);
            
            for(data_iterator it = fscatter.begin(); it != fscatter.end(); it++)
            {
                p.SetPoint(1,  in.eval((*it).X()) );
                (*it)+= p;
            }
        }
        
        void scatter::operator-=(const scatter& in)
        {
            point p = point(fscatter[0].size());
            for(size_t i =0; i< p.size(); i++)
                p.SetPoint(i,1.);
            
            for(data_iterator it = fscatter.begin(); it != fscatter.end(); it++)
            {
                p.SetPoint(1, -1*in.eval((*it).X()) );
                (*it)+= p;
            }
        }
        
        void scatter::operator*=(double a)
        {
            point p = point(fscatter[0].size());
            for(size_t i =0; i< p.size(); i++)
                p.SetPoint(i,1.);
            
            p.SetPoint(1,a);
            
            for(data_iterator it = fscatter.begin(); it != fscatter.end(); it++)
                (*it) *= p;
        }
        
        void scatter::operator/=(double a)
        {
            operator*=(1./a);
        }
        
        void scatter::operator+=(double a)
        {
            point p = point(fscatter[0].size());
            p.SetPoint(1,a);
            
            for(data_iterator it = fscatter.begin(); it != fscatter.end(); it++)
                (*it) += p;
        }
        
        void scatter::operator-=(double a)
        {
            operator+=(-a);
        }
        
#pragma mark • Mathematic estimator
        double scatter::eval(double x) const
        {
            data tmp = std::vector<point>(fscatter);
            std::sort(tmp.begin(), tmp.end(), point_sort);
            
            if(x < tmp[0].X() || x > tmp[tmp.size()-1].X())
                return 0;
            
            if(tmp.size() < 1)
                return 0;
            
            data_iterator it = tmp.end();
            
            double a=0, b=0, fa=0, fb=0, slope=0, cte=0;
            
            if((x - tmp[0].X()) <= (tmp[tmp.size()-1].X() - x))
            {
                it = tmp.begin();
                
                while ((*it).X() < x && it != tmp.end())
                {
                    a = (*it).X();
                    fa= (*it).Y();
                    it++;
                }
                
                if(it == tmp.end())
                    return fa;
                
                b = (*it).X();
                fb= (*it).Y();
            }
            else
            {
                it = --tmp.end();
    
                while ((*it).X() > x && it != tmp.begin())
                {
                    b = (*it).X();
                    fb= (*it).Y();
                    
                    it--;
                }
                
                a = (*it).X();
                fa= (*it).Y();
                
            }
            
            slope = (fb-fa)/(b-a);
            cte   = fb - slope*b;
            
            tmp.clear();
            
            return slope*x + cte;
            
        }
        
    }
}
