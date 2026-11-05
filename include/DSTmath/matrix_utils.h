// -*- mode: c++; -*- 
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under PolyForm Noncommercial 1.0.0
//  You may use, share and adapt this code for noncommercial purposes.
//  Commercial use requires a separate licence from the author.
//  Licence text: https://polyformproject.org/licenses/noncommercial/1.0.0/

#include <vector>
#include <stdexcept>
#include <string>

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
