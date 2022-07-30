// -*- mode: c++; -*- 
/* matrix_utils.h
 * Author(s):     Arnaud Chapon <chapon@cppm.in2p3.fr>
 * Creation date: 2012-02-10
 * Last modified: 2012-08-03
 * 
 * License: 
 * 
 * Description: 
 * 
 * History: 
 * 
 */

#include <vector>

#include "vector_utils.h"

#ifndef __DST_MATH_matrix__
#define __DST_MATH_matrix__

typedef std::vector<std::vector<double> > matrix_t;

namespace DST {
  namespace Math{
  
  class matrix_utils
  {

  private:
    
    bool         _debug_;
    vector_utils _vu_;

  public:

    matrix_t init_from ( double *, size_t, size_t );

    void display ( matrix_t, std::ostream & );
    void display ( matrix_t );

    matrix_t create_unitary_matrix     ( unsigned int );
    matrix_t create_THmatrix           ( unsigned int );
    matrix_t create_jacobian_matrix    ( vector_t );
    matrix_t create_jacobian_tr_matrix ( vector_t );
    matrix_t tr_matrix                 ( matrix_t );
    double   get_trace                 ( matrix_t );
    double   get_determinant           ( matrix_t );
    double   elements_sum              ( matrix_t );
    matrix_t product                   ( matrix_t, double );
    matrix_t product                   ( matrix_t, matrix_t );
    matrix_t matrix_part               ( matrix_t, unsigned int, unsigned int );
    matrix_t comatrix                  ( matrix_t );
    matrix_t matinv                    ( matrix_t );
    void     pivot                     ( matrix_t &, matrix_t & );

    bool is_debug();
    void set_debug ( bool );

    // ctor/dtor:
    matrix_utils();
    matrix_utils ( bool );
    virtual ~matrix_utils();

  };

} // end of namespace utilities
}

#endif // __utilities__matrix_utils_h

// end of matrix_utils.h
