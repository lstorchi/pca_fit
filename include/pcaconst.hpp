#ifndef _INC_PCACONST_H_
#define _INC_PCACONST_H_

#include <stdexcept>
#include <vector>

#define RANGE_CHECK

namespace pca
{
  template<class T> 
  class matrix 
  {
     private:
        std::vector<T> elements; // array of elements
        unsigned int rows;  // number of rows
        unsigned int cols;  // number of columns

        void range_check( unsigned i, unsigned j ) const;

     public:
        matrix(unsigned rows, unsigned cols, const T* elements = 0);
        matrix(const matrix<T>& );

        ~matrix();
    
        matrix<T>& operator=( const matrix<T>& );
    
        bool operator==( const matrix<T>& ) const;

        void clear();

        unsigned int n_rows () const {return rows;};
        unsigned int n_cols () const {return cols;};

        void reset(unsigned rows, unsigned cols);

        T& operator()( unsigned i, unsigned j ) 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements[i*cols+j];
        }

        const T& operator()( unsigned i, unsigned j ) const 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements[i*cols+j];
        }

        const T& element(unsigned i, unsigned j) const 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements[i*cols+j];
        }

        T& element(unsigned i, unsigned j) 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements[i*cols+j];
        }
  };

  template<class T>
  matrix<T>::matrix (unsigned rows, unsigned cols, const T* elements)
  {
    rows = rows;
    cols = cols;
  
    if( rows == 0 || cols == 0 )
      throw std::range_error("attempt to create a degenerate matrix");
  
    if (elements)
      for (unsigned i=0; i<rows*cols; i++)
        this->elements[i] = elements[i];
  };
  
  template<class T>
  matrix<T>::matrix( const matrix<T>& cp )
  {
    rows = cp.rows;
    cols = cp.cols; 
    elements = cp.elements;
  }
  
  template<class T>
  matrix<T>::~matrix()
  {
    clear();
  }
  
  template<class T>
  void matrix<T>::clear() 
  {
    elements.clear();
    rows = 0;
    cols = 0;
  }
  
  template<class T>
  void matrix<T>::reset(unsigned rows, unsigned cols) 
  {
    elements.clear();
    rows = rows;
    cols = cols;
  }
  
  template<class T>
  matrix<T>& pca::matrix<T>::operator=( const pca::matrix<T>& cp )
  {
     if (cp.rows != rows && cp.cols != cols )
        throw std::domain_error("matrix op= not of same order");
     
     for(unsigned i=0; i<rows*cols; i++)
        elements[i] = cp.elements[i];
  
     return *this;
  }
  
  template<class T>
  void matrix<T>::range_check( unsigned i, unsigned j ) const
  {
     if (rows <= i)
       throw std::range_error("matrix access row out of range");
  
     if (cols <= j)
        throw std::range_error("matrix access col out of range");
  }

}

#endif
