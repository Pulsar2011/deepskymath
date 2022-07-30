// -*- mode: c++; -*- 
/* vector_utils.cc
 */

#include <cmath>
#include <stdexcept>

#if __cplusplus < 199711L
#include <cstdarg>
#endif

#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>
#include <algorithm>

#include <DSTmath/vector_utils.h>
#include <stdarg.h>

namespace DST
{
  namespace Math
  {
  bool vector_utils::_debug_ = false;

  bool
  vector_utils::is_debug()
  {
      return DST::Math::vector_utils::_debug_;
  }

  void
  vector_utils::set_debug ( bool debug_ )
  {
    DST::Math::vector_utils::_debug_ = debug_;
  }


  double 
  vector_utils::sum ( const vector_t& v_ )
  {
    double s = 0.;
    for ( vector_t::const_iterator i = v_.begin(); i != v_.end(); i++ )
      s += (*i);
    return s;
  }
	
	double vector_utils::weighted_sum ( const vector_t& v_, const vector_t& w_ )
	{
		if ( v_.size() != w_.size() )
		{
			std::cerr << "vector_utils::weighted_sum: WARNING: size are differents!" << std::endl;
			return 0;
		}
		
		double n = 0;
		for ( unsigned int i = 0; i < v_.size(); i++ )
			n += v_.at ( i ) * w_.at ( i );
		
		return n;
	}

  double
  vector_utils::sum_of_squared ( const vector_t& v_ )
  {
    double s = 0.;
    for ( vector_t::const_iterator i = v_.begin(); i != v_.end(); i++ )
      s += (*i) * (*i);
    return s;
  }


    double
    vector_utils::median ( const vector_t& _v_ )
    {
        vector_t v_ = _v_;
        
        std::sort(v_.begin(), v_.end());
        
        double median = 0;
        if ( v_.size() < 1 )
            throw std::runtime_error ("vector_utils::median: ERROR: Not enough value !");
        else if( v_.size() < 2 )
            median = v_[0];
        else if(!(static_cast<int>( v_.size() )%2))
            median = (v_[v_.size()/2-1] + v_[v_.size()/2])/2.;
        else
            median = v_[v_.size()/2];
        
        v_.clear();
        return median;
    }
    
  double
  vector_utils::mean ( const vector_t& v_ )
  {
    if ( v_.size() < 1 )
      throw std::runtime_error ("vector_utils::mean: ERROR: Not enough value !");
    return sum ( v_ ) / double ( v_.size() );
  }
    
    double vector_utils::trimed_mean ( const vector_t& v_ )
    {
        if ( v_.size() < 1 )
            throw std::runtime_error ("vector_utils::mean: ERROR: Not enough value !");
        
        if(v_.size() < 4)
            return vector_utils::mean(v_);
        
        vector_t cpy = v_;
        std::sort(cpy.begin(), cpy.end());
        
        cpy.erase(cpy.begin());
        cpy.pop_back();
        
        
        double mean = sum ( cpy ) / static_cast<double>( cpy.size() );
        cpy.clear();
        
        return mean;
    }

  double
  vector_utils::quadratic_mean ( const vector_t& v_ )
  {
    if ( v_.size() == 0 )
      throw std::runtime_error ("vector_utils::quadratic_mean: ERROR: Not enough value !");
    return sqrt ( mean_of_squared ( v_ ) );
  }

  double
  vector_utils::mean_of_squared ( const vector_t& v_ )
  {
    if ( v_.size() == 0 )
      throw std::runtime_error ("vector_utils::mean_of_squared: ERROR: Not enough value !");
    return sum_of_squared ( v_ ) / v_.size();
  }

  double
  vector_utils::weighted_mean ( const vector_t& v_, const vector_t& w_ )
  {
    double n = weighted_sum(v_, w_);
    return n / sum ( w_ );
  }

  double
  vector_utils::variance ( const vector_t& v_ )
  {
    double m = mean ( v_ );
    return variance ( v_, m );
  }

  double
  vector_utils::variance ( const vector_t& v_, double m_ )
  {
    double v           = 0.;
    bool   bessel_corr = true;
    for ( vector_t::const_iterator i = v_.begin(); i != v_.end(); i++ )
      v += ( (*i) - m_ ) * ( (*i) - m_ );
    if ( bessel_corr )
      return ( mean_of_squared ( v_ ) - m_ * m_ ) * ( v_.size() / double ( v_.size() - 1. ) );
    else
      return ( mean_of_squared ( v_ ) - m_ * m_ );
  }

  double
  vector_utils::variance2 ( const vector_t& v_ )
  {
    double m = mean ( v_ );
    double v = 0.;

    for ( vector_t::const_iterator i = v_.begin(); i != v_.end(); i++ )
      v += ( (*i) - m ) * ( (*i) - m );
    return v / double ( v_.size() );
  }


  double
  vector_utils::variance_error ( const vector_t& v_ )
  {
    double n = v_.size();
    double v = variance ( v_ );
    double a = 0;

    for ( vector_t::const_iterator i = v_.begin(); i != v_.end(); i++ )
      a += pow ( ( (*i) - mean ( v_ ) ), 4. );
    a /= n;
    double b = ( n - 3 ) / ( n - 1 ) * v * v;
    return sqrt ( ( a - b ) / n );
  }


  double
  vector_utils::covariance ( const vector_t& v_, const vector_t& w_ )
  {
    if ( v_.size() != w_.size() )
      std::cerr << "vector_utils::covariance: WARNING: "
		<< "size are differents!" << std::endl;

    vector_t s;
    for ( size_t i = 0; i < v_.size(); i++ )
      s.push_back ( v_.at(i) * w_.at(i) );
      
      double cov = mean ( s ) - ( mean ( v_ ) * mean ( w_ ) );
      s.clear();
    return cov;
  }


  double
  vector_utils::correlation ( const vector_t& v_, const vector_t& w_ )
  {
    return covariance ( v_, w_ )
      / ( sqrt ( variance2 ( v_ ) ) * sqrt ( variance2 ( w_ ) ) );
  }
    
    
  double
  vector_utils::find_minimum    (vector_t _v_, bool doSort)
  {
      if(doSort)
          return vector_utils::get_minimum(_v_);
      
      double min = _v_[0];
      for(unsigned int i = 0; i < _v_.size(); i++)
          min = (min <= _v_[i])? min : _v_[i];
      
      return min;
  }
  
  double
  vector_utils::find_maximum    (vector_t _v_, bool doSort)
  {
      if(doSort)
          return vector_utils::get_maximum(_v_);
      
      double max = _v_[0];
      for(unsigned int i = 0; i < _v_.size(); i++)
          max = (max >= _v_[i])? max : _v_[i];
      
      return max;
  }
  
  double
  vector_utils::get_minimum     (vector_t _v_)
  {
      std::sort(_v_.begin(), _v_.end());
      return _v_[0];
  }
  
  double
  vector_utils::get_maximum     (vector_t _v_)
  {
      std::sort(_v_.begin(), _v_.end());
      return _v_[_v_.size()-1];
  }
	

  vector_t
  vector_utils::lin_reg ( vector_t x_, vector_t y_, bool display_params_ )
  {
    if ( x_.size() != y_.size() )
      std::cerr << "vector_utils::lin_reg: WARNING: "
		<< "size are differents!" << std::endl;

    vector_t x2;
    vector_t xy;
    vector_t params;

    for ( size_t i = 0; i < x_.size(); i++ )
      {
	x2.push_back ( x_.at(i) * x_.at(i) );
	xy.push_back ( x_.at(i) * y_.at(i) );
      }

    double a = ( sum ( xy ) * double( x_.size() ) - sum ( y_ ) * sum ( x_ ) )
      / ( sum ( x2 ) * double( x_.size() ) - sum ( x_ ) * sum ( x_ ) );
    double b = ( sum ( xy ) * sum ( x_ ) - sum ( y_ ) * sum ( x2 ) )
      / ( sum ( x_ ) * sum ( x_ ) - double( x_.size() ) * sum ( x2 ) );

    double S    = 0.;
    double chi  = 0.;
    double chi2 = 0.;

    for ( size_t i = 0; i < x_.size(); i++ )
      {
	double d = fabs ( y_.at ( i ) - a * x_.at ( i ) - b );
	S += d * d;
	if ( y_.at ( i ) != 0. )
	  {
	    chi = d / sqrt ( fabs ( y_.at ( i ) ) ); // suppose a gaussian error
	    chi2 += chi * chi;
	  }
      }
    double reduced_chi2 = chi2 / ( x_.size() - 2. );

    params.push_back ( a    );
    params.push_back ( b    );
    // params.push_back ( S    );
    params.push_back ( chi2 );
    params.push_back ( reduced_chi2 );

    if ( display_params_ )
      {
	vector_t y_th;
	for ( vector_t::iterator i = x_.begin(); i != x_.end(); i++ )
	  y_th.push_back ( a * (*i) + b );
	
	double cov      = covariance ( y_, y_th );
	double var_y    = variance2  ( y_ );
	double var_y_th = variance2  ( y_th );
	double correlation_factor = cov / ( sqrt ( var_y ) * sqrt ( var_y_th ) );
    
	std::cerr << "vector_utils::lin_reg: NOTICE: y = a * x + b" << std::endl
		  << "|   |-- a     = " << a << std::endl
		  << "|   |-- b     = " << b << std::endl
		  << "|   |-- r2    = " << correlation_factor * correlation_factor << std::endl
		  << "|   `-- sumsq = " << S << std::endl;
	std::cerr << "|-- cov ( y, y_yh ) = " << cov << std::endl;
	std::cerr << "|-- sigma_y         = " << sqrt ( var_y ) << std::endl;
	std::cerr << "|-- sigma_y_th      = " << sqrt ( var_y_th ) << std::endl;
	std::cerr << "`-- rho             = " << correlation_factor << std::endl;
      }

    return params;
  }
  
  vector_t
  vector_utils::lin_reg ( vector_t x_, vector_t y_ )
  {
    return lin_reg ( x_, y_, false );
  }
    
    vector_t
    vector_utils::lin_reg ( vector_t x_, vector_t y_, vector_t s_ )
    {
        return lin_reg ( x_, y_, s_, false );
    }
    
    vector_t
    vector_utils::lin_reg ( vector_t x_, vector_t y_, vector_t s_, bool display_params_ )
    {
        if ( x_.size() != y_.size() || x_.size() != s_.size() )
            std::cerr << "vector_utils::lin_reg: WARNING: "
            << "size are differents!" << std::endl;
        
        vector_t w;
        vector_t x;
        vector_t y;
        vector_t x2;
        vector_t xy;
        vector_t params;
        
        for ( size_t i = 0; i < x_.size(); i++ )
        {
            if(s_.at(i) != 0)
                w .push_back (1./( s_.at(i) * s_.at(i) ));
            else
                w .push_back (1.);
            x .push_back ( x_.at(i) / s_.at(i) / s_.at(i) );
            y .push_back ( y_.at(i) / s_.at(i) / s_.at(i) );
            x2.push_back ( x_.at(i) * x_.at(i) / s_.at(i) / s_.at(i) );
            xy.push_back ( x_.at(i) * y_.at(i) / s_.at(i) / s_.at(i) );
        }
        
        double a = ( sum(x)*sum(y) - sum(w)*sum(xy) )/( sum(x)*sum(x) - sum(x2)*sum(w) );
        double b = ( sum(y) - a*sum(x) )/sum(w);
        double da2= fabs( 1. / (sum(w) * ( mean(x2) - mean(x)*mean(x) )));
        double db2= mean(x2) * da2;
        double cov= - mean(x2) / (sum(w) * (mean(x) - mean(x)*mean(x)));
                        
        
        double S    = 0.;
        double chi  = 0.;
        double chi2 = 0.;
        
        for ( size_t i = 0; i < x_.size(); i++ )
        {
            double d = fabs ( y_.at ( i ) - a * x_.at ( i ) - b );
            S += d * d;
            if ( s_.at ( i ) != 0. )
            {
                chi = d / ( s_.at(i)*s_.at(i) );
                chi2 += chi * chi;
            }
        }
        double reduced_chi2 = chi2 / ( x_.size() - 2. );
        
        params.push_back ( a            );
        params.push_back ( b            );
        params.push_back ( sqrt(da2)    );
        params.push_back ( sqrt(db2)    );
        params.push_back ( cov          );
        params.push_back ( chi2         );
        params.push_back ( reduced_chi2 );
        
        if ( display_params_ )
        {
            vector_t y_th;
            for ( vector_t::iterator i = x_.begin(); i != x_.end(); i++ )
                y_th.push_back ( a * (*i) + b );
            
            double cov      = covariance ( y_, y_th );
            double var_y    = variance2  ( y_ );
            double var_y_th = variance2  ( y_th );
            double correlation_factor = cov / ( sqrt ( var_y ) * sqrt ( var_y_th ) );
            
            std::cerr << "vector_utils::lin_reg: NOTICE: y = a * x + b" << std::endl
            << "|   |-- a       = " << a <<" +/- "<<sqrt(da2)<< std::endl
            << "|   |-- b       = " << b <<" +/- "<<sqrt(db2)<< std::endl
            << "|   |-- cov(a,b)= " << cov << std::endl
            << "|   |-- chi2    = " << chi2 <<" / "<<( x_.size() - 2. )<<" -> "<<reduced_chi2<<std::endl
            << "|   |-- r2      = " << correlation_factor * correlation_factor << std::endl
            << "|   `-- sumsq   = " << S << std::endl;
            std::cerr << "|-- cov ( y, y_yh ) = " << cov << std::endl;
            std::cerr << "|-- sigma_y         = " << sqrt ( var_y ) << std::endl;
            std::cerr << "|-- sigma_y_th      = " << sqrt ( var_y_th ) << std::endl;
            std::cerr << "`-- rho             = " << correlation_factor << std::endl;
        }
        
        return params;
    }

	
  vector_t
  vector_utils::zero_lin_reg ( vector_t x_, vector_t y_, bool display_params_ )
  {
    if ( x_.size() != y_.size() )
      std::cerr << "vector_utils::lin_reg: WARNING: "
		<< "size are differents!" << std::endl;

    vector_t params;

    double a = sum ( y_ ) / sum ( x_ );
    double b = 0.;

    double S = 0.;
    for ( size_t i = 0; i < x_.size(); i++ )
      {
	double d = y_.at ( i ) - a * x_.at ( i ) - b;
	S += d * d;	  
      }

    params.push_back ( a );
    params.push_back ( b );
    params.push_back ( S );

    if ( display_params_ )
      {
	vector_t y_th;
	for ( vector_t::iterator i = x_.begin(); i != x_.end(); i++ )
	  y_th.push_back ( a * (*i) + b );
	
	double cov      = covariance ( y_, y_th );
	double var_y    = variance2  ( y_ );
	double var_y_th = variance2  ( y_th );
	double correlation_factor = cov / ( sqrt ( var_y ) * sqrt ( var_y_th ) );
    
	std::cerr << "vector_utils::lin_reg: NOTICE: y = a * x + b" << std::endl
		  << "|   |-- a     = " << a << std::endl
		  << "|   |-- b     = " << b << std::endl
		  << "|   |-- r2    = " << correlation_factor * correlation_factor << std::endl
		  << "|   `-- sumsq = " << S << std::endl;
	std::cerr << "|-- cov ( y, y_yh ) = " << cov << std::endl;
	std::cerr << "|-- sigma_y         = " << sqrt ( var_y ) << std::endl;
	std::cerr << "|-- sigma_y_th      = " << sqrt ( var_y_th ) << std::endl;
	std::cerr << "`-- rho             = " << correlation_factor << std::endl;
      }

    return params;
  }
    
  
  vector_t
  vector_utils::zero_lin_reg ( vector_t x_, vector_t y_ )
  {
    return zero_lin_reg ( x_, y_, false );
  }


  int
  vector_utils::dbl_counter ( const vector_t& v_, double x_ )
  {
    int a = 0;
    for ( vector_t::const_iterator i = v_.begin(); i != v_.end(); i++ )
      {
	if ( int ( x_ + 0.5 ) == int ( (*i) + 0.5 ) )
	  a++;
      }
    return a;
  }

  std::string
  vector_utils::to_string ( vector_t v_ )
  {
    std::ostringstream oss;
    oss.precision(6);
    oss << "[" << *(v_.begin());
    for ( unsigned int i = 1; i < v_.size(); i++ )
      oss << "," << v_.at ( i );
    oss << "]";
    /*
    std::cerr << "vector_utils::to_char: DEBUG: "
	      << oss.str() << std::endl;
    */
    return oss.str();
  }
	
  vector_t vector_utils::multiply (const void * _v_, ...)
 {
	 const vector_t *arg = static_cast<const vector_t*>( _v_ );
	 
	 vector_t out;
	 
	 va_list argptr;
	 va_start(argptr,_v_);
	 
	 while(arg != NULL)
	 {
		 if(arg->size() > 1)
		 {
			 if(out.size() < 1)
				 for(unsigned int iv = 0 ; iv < arg->size(); iv++)
					 out.push_back(arg->at(iv));
			 else
				 for(unsigned int iv = 0 ; iv < arg->size(); iv++)
					 out[iv] *= arg->at(iv);
		 }
		 arg = va_arg(argptr, const vector_t *);
	 }
	 
	 va_end(argptr);
	
	 return out;

 }
	
	vector_t vector_utils::divide (void * _v_, ...)
 {
	 vector_t *arg = static_cast<vector_t*>( _v_ );
	 
	 vector_t out;
	 
	 va_list argptr;
	 va_start(argptr,_v_);
	 
	 while(arg != NULL)
	 {
		 if(arg->size() > 1)
		 {
			 if(out.size() < 1)
				 for(unsigned int iv = 0 ; iv < arg->size(); iv++)
					 out.push_back(arg->at(iv));
			 else
				 for(unsigned int iv = 0 ; iv < arg->size(); iv++)
					 out[iv] /= arg->at(iv);
		 }
		 arg = va_arg(argptr, vector_t *);
	 }
	 
	 va_end(argptr);
	 
	 return out;
	 
 }
	
	vector_t vector_utils::add (void * _v_, ...)
 {
	 vector_t *arg = static_cast<vector_t*>( _v_ );
	 
	 vector_t out;
	 
	 va_list argptr;
	 va_start(argptr,_v_);
	 
	 while(arg != NULL)
	 {
		 if(arg->size() > 1)
		 {
			 if(out.size() < 1)
				 for(unsigned int iv = 0 ; iv < arg->size(); iv++)
					 out.push_back(arg->at(iv));
			 else
				 for(unsigned int iv = 0 ; iv < arg->size(); iv++)
					 out[iv] += arg->at(iv);
		 }
		 arg = va_arg(argptr, vector_t *);
	 }
	 
	 va_end(argptr);
	 
	 return out;
	 
 }

  void
  vector_utils::display ( vector_t v_, std::ostream & out_ = std::cerr )
  {
    for ( vector_t::iterator i = v_.begin(); i != v_.end(); i++ )
      if ( (*i) == 0 || (*i) == -0 )
	out_ << "0\t";
      else
	out_ << (*i) << "\t";
  }

  void
  vector_utils::display ( vector_t v_ )
  {
    display ( v_, std::cerr );
  }


  // ctor/dtor:
  vector_utils::vector_utils ( bool debug_ )
  {
    _debug_ = debug_;
  }
  
  vector_utils::vector_utils()
  {
    _debug_ = false;
  }
  
  vector_utils::~vector_utils()
  {
  }


} 
}
