#ifndef _INC_PCACONST_H_
#define _INC_PCACONST_H_

#include <stdexcept>
#include <vector>

#define RANGE_CHECK

namespace pca
{
  template<typename T> 
  class matrixpcaconst
  {
     public:
        enum const_type 
        {
          QVEC,
          CMTX,
          AMTX,
          KVEC,
          NONDEF
        };

        enum sector_type 
        {
          BARREL,
          HYBRID,
          ENDCAP,
          UNKN
        };

        enum plane_type 
        {
          RPHI,
          RZ,
          UNDEF
        };

     private:
        std::vector<T> elements_; 
        unsigned int rows_;  
        unsigned int cols_;  

        const_type ctype_;
        sector_type stype_;
        plane_type ptype_;
        int towerid_;
        double ptmin_, ptmax_, etamin_, etamax_;

        void range_check (unsigned i, unsigned j) const;

     public:
        matrixpcaconst(unsigned rows, unsigned cols, const T* elements = 0);
        matrixpcaconst(const matrixpcaconst<T>& );

        ~matrixpcaconst();
    
        matrixpcaconst<T>& operator=( const matrixpcaconst<T>& );
    
        bool operator==( const matrixpcaconst<T>& ) const;

        void set_const_type (const_type in) {ctype_ = in;};
        const_type get_const_type () const {return ctype_;};

        void set_sector_type (sector_type in) {stype_ = in;};
        sector_type get_sector_type () const {return stype_;};

        void set_plane_type (plane_type in) {ptype_ = in;};
        plane_type get_plane_type () const {return ptype_;};

        void set_towerid (int in) {towerid_ = in;};
        int get_towerid () const {return towerid_;};

        void set_ptrange (double min, double max) 
        {
          ptmin_ = min;
          ptmax_ = max;
        };
        void get_ptrange (double & min, double & max) const
        {
          min = ptmin_;
          max = ptmax_;
        };

        void set_etarange (double min, double max) 
        {
          etamin_ = min;
          etamax_ = max;
        };
        void get_etarange (double & min, double & max) const
        {
          min = etamin_;
          min = etamax_;
        };

        void clear();

        unsigned int n_rows () const {return rows_;};
        unsigned int n_cols () const {return cols_;};

        void reset(unsigned rows, unsigned cols);

        T& operator()( unsigned i, unsigned j ) 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements_[i*cols_+j];
        }

        const T& operator()( unsigned i, unsigned j ) const 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements_[i*cols_+j];
        }

        const T& element(unsigned i, unsigned j) const 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements_[i*cols_+j];
        }

        T& element(unsigned i, unsigned j) 
        {
           #ifdef RANGE_CHECK
           range_check(i,j);
           #endif
           return elements_[i*cols_+j];
        }
  };

  template<class T>
  matrixpcaconst<T>::matrixpcaconst (unsigned rows, unsigned cols, 
      const T* elements)
  {
    rows_ = rows;
    cols_ = cols;

    ctype_ = NONDEF;
    stype_ = UNKN;
    ptype_ = UNDEF;

    towerid_ = 0;
    ptmin_ = 0.0;
    ptmax_ = 0.0;
    etamin_ = 0.0;
    etamax_ = 0.0;
  
    if( rows_ == 0 || cols_ == 0 )
      throw std::range_error("attempt to create a degenerate matrix");
  
    if (elements)
      for (unsigned i=0; i<rows_*cols_; i++)
        this->elements_.push_back(elements[i]);
    else
      for (unsigned i=0; i<rows_*cols_; i++)
        this->elements_.push_back((T)0);
  };
  
  template<class T>
  matrixpcaconst<T>::matrixpcaconst( const matrixpcaconst<T>& cp )
  {
    rows_ = cp.rows_;
    cols_ = cp.cols_; 
    elements_ = cp.elements_;

    ctype_ = cp.ctype_;
    stype_ = cp.stype_;
    ptype_ = cp.ptype_;

    ptmin_ = cp.ptmin_;
    ptmax_ = cp.ptmax_;
    etamin_ = cp.etamin_;
    etamax_ = cp.etamax_;
    towerid_ = cp.towerid_;
  }
  
  template<class T>
  matrixpcaconst<T>::~matrixpcaconst()
  {
    clear();
  }
  
  template<class T>
  void matrixpcaconst<T>::clear() 
  {
    elements_.clear();
    rows_ = 0;
    cols_ = 0;

    ctype_ = NONDEF;
    stype_ = UNKN;
    ptype_ = UNDEF;

    towerid_ = 0;
    ptmin_ = 0.0;
    ptmax_ = 0.0;
    etamin_ = 0.0;
    etamax_ = 0.0;
  }
  
  template<class T>
  void matrixpcaconst<T>::reset(unsigned rows, unsigned cols) 
  {
    clear();

    rows_ = rows;
    cols_ = cols;

    for (unsigned i=0; i<rows_*cols_; i++)
      this->elements_.push_back((T)0);
  }
  
  template<class T>
  matrixpcaconst<T>& pca::matrixpcaconst<T>::operator=( 
      const pca::matrixpcaconst<T>& cp )
  {
     if (cp.rows_ != rows_ && cp.cols_ != cols_ )
        throw std::domain_error("matrix op= not of same order");
     
     for(unsigned i=0; i<rows_*cols_; i++)
       elements_.push_back(cp.elements_[i]);

     ctype_ = cp.ctype_;
     stype_ = cp.stype_;
     ptype_ = cp.ptype_;
     
     ptmin_ = cp.ptmin_;
     ptmax_ = cp.ptmax_;
     etamin_ = cp.etamin_;
     etamax_ = cp.etamax_;
     towerid_ = cp.towerid_;
  
     return *this;
  }
  
  template<class T>
  void matrixpcaconst<T>::range_check(unsigned i, unsigned j) const
  {
     if (rows_ <= i)
       throw std::range_error("matrix access row out of range");
  
     if (cols_ <= j)
        throw std::range_error("matrix access col out of range");
  }

}

#endif
