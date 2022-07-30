// -*- mode: c++; -*- 
/* matrix_utils.cc
 */

#include<iostream>
#include<cmath>
#include<stdexcept>

#include <DSTmath/matrix_utils.h>

namespace DST {
	namespace Math
	{

  bool
  matrix_utils::is_debug()
  {
    return _debug_;
  }

  void
  matrix_utils::set_debug ( bool debug_ )
  {
    _debug_ = debug_;
  }


  matrix_t
  matrix_utils::init_from ( double * tab_, size_t n, size_t m )
  {
    vector_t line;
    matrix_t matrix;

    for ( size_t i = 0; i < n; i++ )
      {
	for ( size_t j = 0; j < m; j++ )
	  {
	    line.push_back ( tab_[m * j + i] );
	  }
	matrix.push_back ( line );
	line.clear();
      }
    return matrix;
  }


  matrix_t
  matrix_utils::create_unitary_matrix ( unsigned int n_ )
  {
    vector_t line;
    matrix_t matrix;
    for ( size_t i = 0; i < n_; i++ )
      {
	for ( size_t k = 0; k < n_; k++ )
	  {
	    if ( i == k )
	      line.push_back ( 1. );
	    else
	      line.push_back ( 0. );
	  }
	matrix.push_back ( line );
	line.clear();
      }
    return matrix;
  }


  matrix_t
  matrix_utils::create_THmatrix ( unsigned int n_ )
  {
    vector_t line;
    matrix_t matrix;
    for ( size_t i = 0; i < n_; i++ )
      {
	for ( size_t k = 0; k < n_; k++ )
	  {
	    if ( i == k )
	      line.push_back ( 2. );
	    else if ( i == k+1 || i == k-1 )
	      line.push_back ( -1. );
	    else
	      line.push_back ( 0. );
	  }
	matrix.push_back ( line );
	line.clear();
      }
    return matrix;
  }

  double
  matrix_utils::get_trace ( matrix_t m_ )
  {
    if ( m_.size() != m_[0].size() )
      throw std::runtime_error ( "matrix_utils::get_trace: ERROR: matrix is not NxN!" );
    unsigned int dim = m_.size();	
    vector_t t;
    for ( size_t i = 0; i < dim; i++ )
      t.push_back ( m_[i][i] );
    return _vu_.sum ( t );
  }


  double
  matrix_utils::get_determinant ( matrix_t m_ )
  {
    if ( m_.size() != m_[0].size() )
      {
	std::cerr << "matrix_utils::get_determinant: ERROR: matrix is "
		  << m_.size() << "x" << m_[0].size() << ":"
		  << std::endl;
	display ( m_ );
	throw std::runtime_error ( "matrix_utils::get_determinant: ERROR: matrix is not NxN!" );
      }

    int    dim = m_.size();
    double det = 0;
  
    matrix_t B;
    if ( dim == 1 )
      return m_[0][0];
    if ( dim == 2 )
      return m_[0][0] * m_[1][1] - m_[0][1] * m_[1][0];
    else
      {
	for ( int j = 0; j < dim; j++ )
	  {
	    B = matrix_part ( m_, 0, j );
	    if ( j % 2 == 0 )
	      det += m_[0][j] * get_determinant ( B );
	    else
	      det -= m_[0][j] * get_determinant ( B );
	  }
      }
    return det; 
  }
  
  double
  matrix_utils::elements_sum ( matrix_t m_ )
  {
    vector_t result;
    unsigned int xdim = m_.size();
    unsigned int ydim = m_[0].size();
    for ( size_t i = 0; i < xdim; i++ )
      for ( size_t j = 0; j < ydim; j++ )
	result.push_back ( m_[i][j] );
    return _vu_.sum ( result );
  }


  matrix_t
  matrix_utils::product ( matrix_t m_, double k_ )
  {
    vector_t line;
    matrix_t nm;
    unsigned int dim = m_.size();
    for ( size_t i = 0; i < dim; i++ )
      {
	for ( size_t j = 0; j < dim; j++ )
	  {
	    line.push_back ( m_[i][j] * k_ );
	  }
	nm.push_back ( line );
	line.clear();
      }
    return nm;
  }

  matrix_t
  matrix_utils::product ( matrix_t m1_, matrix_t m2_ )
  {
    if ( m1_[0].size() != m2_.size() )
      throw std::runtime_error ( "matrix_utils::product: ERROR: bad size." );
    
    vector_t line1, line2;
    matrix_t matrix;

    unsigned int n  = m1_.size();
    unsigned int n1 = m1_[0].size();
    unsigned int n2 = m2_[0].size();

    for ( size_t i = 0; i < n; i++ )
      {
	for ( size_t k = 0; k < n2; k++ )
	  {
	    for ( size_t j = 0; j < n1; j++ )
	      {
		line1.push_back ( m1_[i][j] * m2_[j][k] );
	      }
	    line2.push_back ( _vu_.sum ( line1 ) );
	    line1.clear();
	  }
	matrix.push_back ( line2 );
	line2.clear();
      }
    return matrix;
  }


  matrix_t
  matrix_utils::create_jacobian_matrix ( vector_t t_ )
  {
    matrix_t matrix;
    vector_t line;

    for ( size_t i = 0; i < t_.size(); i++ )
      {
	line.push_back ( t_.at ( i ) );
	line.push_back ( 1. );
	matrix.push_back ( line );
	line.clear();
      }
    return matrix;	
  }
	
  matrix_t
  matrix_utils::create_jacobian_tr_matrix ( vector_t t_ )
  {
    matrix_t matrix;
    vector_t line;

    for ( size_t i = 0; i < t_.size(); i++ )
      line.push_back ( 1. );
    matrix.push_back ( t_ );
    matrix.push_back ( line );
    return matrix;	
  }


  matrix_t
  matrix_utils::tr_matrix ( matrix_t m_ )
  {
    matrix_t mt;
    vector_t line;
    unsigned int xdim = m_.size();
    unsigned int ydim = m_[0].size();
    for ( size_t i = 0; i < ydim; i++ )
      {
	for ( size_t j = 0; j < xdim; j++ )
	  {
	    line.push_back ( m_[j][i] );
	  }
	mt.push_back ( line );
	line.clear();
      }
    return mt;	
  }


  matrix_t
  matrix_utils::matrix_part ( matrix_t m_, unsigned int i_, unsigned int j_ )
  {
    matrix_t mp;
    vector_t line;
    unsigned int xdim = m_.size();
    unsigned int ydim = m_[0].size();
    for ( size_t i = 0; i < xdim; i++ )
      {
	for ( size_t j = 0; j < ydim; j++ )
	  {
	    if ( i != i_ && j != j_ )
	      line.push_back ( m_[i][j] );
	  }
	if ( line.size() > 0 )
	  mp.push_back ( line );
	line.clear();
      }
    return mp;	
  }


  matrix_t
  matrix_utils::comatrix ( matrix_t m_ )
  {
    if ( m_.size() != m_[0].size() )
      throw std::runtime_error ( "matrix_utils::comatrix: ERROR: matrix is not NxN!" );
    
    matrix_t nm, mp;
    vector_t line;
    unsigned int xdim = m_.size();
    unsigned int ydim = m_[0].size();

    for ( unsigned int i = 0; i < xdim; i++ )
      {
	for( unsigned int j = 0; j < ydim; j++ )
	  {
	    mp = matrix_part ( m_, i, j );

	    if ( ( i + j ) % 2 == 0 )
	      line.push_back ( get_determinant ( mp ) );
	    else
	      line.push_back ( -1. * get_determinant ( mp ) );
	  }
	nm.push_back ( line );
	line.clear();
      }
    return nm;
  }
  
  matrix_t
  matrix_utils::matinv ( matrix_t m_ )
  {
    if ( m_.size() != m_[0].size() )
      throw std::runtime_error ( "matrix_utils::matinv: ERROR: matrix is not NxN!" );
    
    double det;
    matrix_t coA, t_coA, inv_A;

    det = get_determinant ( m_ );
    if ( det == 0 )
      throw std::runtime_error ( "matrix_utils::matinv: ERROR: determinant is null." );      

    coA   = comatrix  ( m_ );
    t_coA = tr_matrix ( coA );
    return product ( t_coA, 1./det );
  } 


  void
  matrix_utils::pivot ( matrix_t & m1_, matrix_t & m2_ )
  {
    unsigned int n = m1_.size();
    if ( n != m1_[0].size() || n != m2_.size() )
      throw std::runtime_error ( "matrix_utils::pivot: ERROR: bad size." );      

    for ( size_t j = 0; j < n; j++ )
      {
	double jj = m1_[j][j];

	for ( size_t k = 0; k < n; k++ )
	  {
	    m1_[j][k] = m1_[j][k] * 1.0 / jj;
	    m2_[j][k] = m2_[j][k] / jj;
	  }
	for ( size_t i = j + 1; i < n; i++ )
	  {
	    double ii = m1_[i][j];
	    for ( size_t k = 0; k < n; k++ )
	      {
		m1_[i][k] = m1_[i][k] - ii * m1_[j][k];
		m2_[i][k] = m2_[i][k] - ii * m2_[j][k];
	      }
	  }
      }
    unsigned int count = 0;
    for ( size_t j = n - 1; j > 0; j-- )
      {
	for ( size_t i = j - 1; i != 0; i-- )
	  {
	    double el = m1_[i][n-1-count];

	    for ( size_t k = 0; k < n; k++ )
	      {
		m2_[i][k] = m2_[i][k] - el * m2_[j][k];
		m1_[i][k] = m1_[i][k] - el * m1_[j][k];
	      }
	  }
	count++;
      }
  }

  // GS: verified!
  /*
  double
  matrix_utils::matinv ( int lin_exg_, int col_exg_, matrix_t m1_, matrix_t m2_ )
  {
    double pc[lin_exg_], pl[lin_exg_], cs[lin_exg_];
    double det = 1.; 

    for ( size_t i = 0; i < lin_exg_; i++ )
      {
        pc[i]= 0.;
        pl[i]= 0.;
        cs[i]= 0.;
      }

    int iflagdet = 0;
    int iflagpiv = 0;

    for ( size_t k = 0; k < lin_exg_; k++ )
      {
	// searching greatest pivot 
	if ( iflagdet == 0 ) 
	  {
	    double pv  = m1_[k][k];
	    int    ik  = k;
	    int    jk  = k;
	  
	    for ( size_t i = k; i < lin_exg_; i++ )
	      for ( size_t j = k; j < lin_exg_; j++ ) 
		{     
		  if ( fabs ( m1_[i][j] ) > fabs ( pv ) )
		    { 
		      pv = m1_[i][j];   
		      ik = i;
		      jk = j;
		    }                                                               
		} 
	    // pivot location is ik, jk
	    pc[k] = jk;  
	    pl[k] = ik;

	    // if det == 0, ERROR MESSAGE and STOP                                           
	    if ( ik != k )
	      det = - det;	  
	    if ( jk != k )
	      det = - det;

	    det = det * pv;
	    if ( fabs ( det ) < 1.e-6 )
	      iflagdet = 1;

	    if ( iflagdet != 1 )
	      {   
		// positionning pivot in k, k:		
		if ( ik != k )                                                         
		  for ( size_t i = 0; i < lin_exg_; i++ )
		    {   
		      // exchange lines ik and K of matrix m1_: 
		      double tmp = m1_[ik][i];
		      m1_[ik][i] = m1_[k][i];
		      m1_[k][i]  = tmp;
		    } 
	      
		if ( col_exg_ != 0 )
		  for ( size_t i = 0; i < col_exg_; i++ )
		    {
		      double tmp = m2_[ik][i];
		      m2_[ik][i] = m2_[k][i];
		      m2_[k][i]  = tmp;
		    }

		// pivot is at correct line
		if ( jk != k )
		  for ( size_t i = 0; i < lin_exg_; i++ )
		    {
		      // exchange columns jk and k of matrix m1_:
		      double tmp = m1_[i][jk];
		      m1_[i][jk] = m1_[i][k];
		      m1_[i][k]  = tmp;
		    }
		// pivot is at correct column and is located in k, k.
		for ( size_t i = 0; i < lin_exg_; i++ )
		  {
		    cs[i]     = m1_[i][k];
		    m1_[i][k] = 0.;
		  }
                                                                               
		cs[k]    = 0.;
		m1_[k][k] = 1.;

		// line k of matrix m1_ is modified:
		if ( fabs ( pv ) < 1.e-6 )
		  iflagpiv = 1;
	      
		if ( iflagpiv != 1 )
		  {
		    for ( size_t i = 0; i < lin_exg_; i++ )
		      m1_[k][i] = m1_[k][i] / pv;
		    if ( col_exg_ != 0 )
		      for ( size_t i = 0; i < col_exg_; i++ )
			m2_[k][i] = m2_[k][i] / pv;
		  
		    // other lines of matrix m1_ are modified:		  
		    for ( size_t j = 0; j < lin_exg_; j++ )
		      {
			if ( j == k )
			  j++;                                                  
			for ( size_t i = 0; i < lin_exg_; i++ )
			  // line j of matrix m1_ is modified:
			  m1_[j][i] = m1_[j][i] - cs[j] * m1_[k][i];
			if ( col_exg_ != 0 )
			  for ( size_t i = 0; i < col_exg_; i++ )
			    m2_[j][i] = m2_[j][i] - cs[j] * m2_[k][i];
		      }
		  }
	      }
	  }
      }

    if ( ( iflagdet == 0 ) && ( iflagpiv == 0 ) )
      {
	bool do_exg = true;

	for ( int i = lin_exg_ - 1; i >= 0; i-- )
	  {
	    int ik = (int) pc[i];
	    if ( ik == i )
	      do_exg = false;

	    if ( do_exg )
	      {
		// exchange lines i and pc(i) of matrix m1_
		for ( size_t j = 0; j < lin_exg_; j++ )
		  {
		    double tmp = m1_[i][j];
		    m1_[i][j]  = m1_[ik][j];
		    m1_[ik][j] = tmp;
		  }
		if ( col_exg_ != 0 )
		  for ( size_t j = 0; j < col_exg_; j++ )
		    {
		      double tmp = m2_[i][j];
		      m2_[i][j]  = m2_[ik][j];
		      m2_[ik][j] = tmp;
		    }                                                                 
	      }
	  }  
	
	do_exg = true;

	// exchange columns:  	  
	for ( int j = lin_exg_ - 1; j >= 0; j-- )
	  {     
	    int jk = (int) pl[j];
	    if ( jk == j )
	      do_exg = false;

	    if ( do_exg )
	      {
		// exchange columns j and pl(J) of matrix m1_:
		for ( size_t i = 0; i < lin_exg_; i++ )
		  {
		    double tmp = m1_[i][j];    
		    m1_[i][j]  = m1_[i][jk];  
		    m1_[i][jk] = tmp;   
		  } 
	      }
	  }                                                                     
      }
    return det;
  }
  */

  void
  matrix_utils::display ( matrix_t m_, std::ostream & out_ = std::cerr )
  {
    for ( size_t i = 0; i < m_.size(); i++ )
      {
	_vu_.display ( m_.at(i), out_ );
	out_ << std::endl;
      }
  }

  void
  matrix_utils::display ( matrix_t m_ )
  {
    display ( m_, std::cerr );
  }


  // ctor/dtor:
  matrix_utils::matrix_utils ( bool debug_ )
  {
    _debug_ = debug_;
  }
  
  matrix_utils::matrix_utils()
  {
    _debug_ = false;
  }
  
  matrix_utils::~matrix_utils()
  {
  }

}
} // end of namespace utilities

// end of matrix_utils.cc
