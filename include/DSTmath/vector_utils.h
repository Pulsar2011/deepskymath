// -*- mode: c++; -*- 
/* vector_utils.h
 * Author(s):     Arnaud Chapon <chapon@cppm.in2p3.fr>
 * Creation date: 2012-02-09
 * Last modified: 2013-04-04
 * 
 * License: 
 * 
 * Description: 
 * 
 * History: 
 * 
 */

#include <sstream>
#include <iostream>
#include <vector>

#ifndef __utilities__vector_utils_h
#define __utilities__vector_utils_h 1

typedef std::vector<double> vector_t;

namespace DST
{
  
  class vector_utils
  {

  private:

    static bool _debug_;

  public:

    static double sum             ( const vector_t& );
    static double weighted_sum    ( const vector_t&, const vector_t& );
    static double sum_of_squared  ( const vector_t& );
    static double mean            ( const vector_t& );
    static double trimed_mean     ( const vector_t& );
    static double median          ( const vector_t& );
    static double mean_of_squared ( const vector_t& );
    static double quadratic_mean  ( const vector_t& );
    static double weighted_mean   ( const vector_t&, const vector_t& );
    static double variance        ( const vector_t& );
    static double variance        ( const vector_t&, double );
    static double variance2       ( const vector_t& );
    static double variance_error  ( const vector_t& );
    static double covariance      ( const vector_t&, const vector_t& );
    static double correlation     ( const vector_t&, const vector_t& );
    static int    dbl_counter     ( const vector_t&, double   );
      
    static double find_minimum    (vector_t, bool doSort = false);
    static double find_maximum    (vector_t, bool doSort = false);
      
    static double get_minimum     (vector_t);
    static double get_maximum     (vector_t);

    static vector_t lin_reg      ( vector_t, vector_t, bool );
    static vector_t lin_reg      ( vector_t, vector_t );
    static vector_t lin_reg      ( vector_t, vector_t, vector_t, bool );
    static vector_t lin_reg      ( vector_t, vector_t, vector_t );
    static vector_t zero_lin_reg ( vector_t, vector_t, bool );
    static vector_t zero_lin_reg ( vector_t, vector_t );
      
	static vector_t multiply (const void*, ...);
	static vector_t divide (void*, ...);
	static vector_t add (void*, ...);

    static std::string to_string ( vector_t );

    static void display ( vector_t, std::ostream & );
    static void display ( vector_t );

    static bool is_debug();
    static void set_debug ( bool );

    // ctor/dtor:
    vector_utils();
    vector_utils ( bool );
    virtual ~vector_utils();

  };

}; // end of namespace utilities

#endif // __utilities__vector_utils_h

// end of vector_utils.h
