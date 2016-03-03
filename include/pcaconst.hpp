#ifndef _INC_PCACONST_H_
#define _INC_PCACONST_H_

#include <stdexcept>
#include <vector>

#define RANGE_CHECK

/* L. Storchi: quick and dirty PCA const store class */

template< class T >
struct TypeIsInt32_t
{
  static const bool value = false;
};

template< class T >
struct TypeIsDouble
{
  static const bool value = false;
};

template<>
struct TypeIsInt32_t< int32_t >
{
  static const bool value = true;
};

template<>
struct TypeIsDouble< double >
{
  static const bool value = true;
};

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

        enum ttype 
        {
          FLOATPT,
          INTEGPT,
          NOONE
        };

     private:
        std::vector<T> elements_; 
        unsigned int rows_;  
        unsigned int cols_;  

        ttype ttype_;
        const_type ctype_;
        sector_type stype_;
        plane_type ptype_;
        int towerid_;
        double ptmin_, ptmax_, etamin_, etamax_;
        int chargesign_;
        std::string layersids_;

        void range_check (unsigned i, unsigned j) const;

     public:
        matrixpcaconst(unsigned rows, unsigned cols, const T* elements = 0);
        matrixpcaconst(const matrixpcaconst<T>& );

        ~matrixpcaconst();
    
        matrixpcaconst<T>& operator=( const matrixpcaconst<T>& );
    
        bool operator==( const matrixpcaconst<T>& ) const;

        void set_ttype (ttype in) {ttype_ = in;};
        ttype get_ttype () const {return ttype_;};

        void set_const_type (const_type in) {ctype_ = in;};
        const_type get_const_type () const {return ctype_;};

        void set_sector_type (sector_type in) {stype_ = in;};
        sector_type get_sector_type () const {return stype_;};

        void set_plane_type (plane_type in) {ptype_ = in;};
        plane_type get_plane_type () const {return ptype_;};

        void set_towerid (int in) {towerid_ = in;};
        int get_towerid () const {return towerid_;};

        void set_chargesign (int in) {chargesign_ = in;};
        int get_chargesign () const {return chargesign_;};

        void set_layersids (const char * in) {layersids_ = in;};
        std::string get_layersids () const {return layersids_;};

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
          max = etamax_;
        };

        bool compare_nonelements (const matrixpcaconst<T>&) const;

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

        static std::string const_type_to_string (matrixpcaconst<T>::const_type in)
        {
          switch (in)
          {
            case matrixpcaconst<T>::QVEC:
              return "QVEC";
            case matrixpcaconst<T>::CMTX:
              return "CMTX";
            case matrixpcaconst<T>::AMTX:
              return "AMTX";
            case matrixpcaconst<T>::KVEC:
              return "KVEC";
            case matrixpcaconst<T>::NONDEF:
              return "NONDEF";
          }
        
          return "";
        }

        static matrixpcaconst<T>::const_type string_to_const_type (
            const std::string & in)
        {
          if (in == "QVEC")
            return matrixpcaconst<T>::QVEC;
          else if (in == "CMTX")
            return matrixpcaconst<T>::CMTX;
          else if (in == "AMTX")
            return matrixpcaconst<T>::AMTX;
          else if (in == "KVEC")
            return matrixpcaconst<T>::KVEC;

          return matrixpcaconst<T>::NONDEF;
        }


        static std::string sector_type_to_string (matrixpcaconst<T>::sector_type in)
        {
          switch (in)
          {
            case matrixpcaconst<T>::BARREL:
              return "BARREL";
            case matrixpcaconst<T>::HYBRID:
              return "HYBRID";
            case matrixpcaconst<T>::ENDCAP:
              return "ENDCAP";
            case matrixpcaconst<T>::UNKN:
              return "UNKN";
          }
        
          return "";
        }

        static matrixpcaconst<T>::sector_type string_to_sector_type (
            const std::string & in)
        {
          if (in == "BARREL")
            return matrixpcaconst<T>::BARREL;
          else if (in == "HYBRID")
            return matrixpcaconst<T>::HYBRID;
          else if (in == "ENDCAP")
            return matrixpcaconst<T>::ENDCAP;

          return matrixpcaconst<T>::UNKN;
        }


        static std::string plane_type_to_string (matrixpcaconst<T>::plane_type in)
        {
          switch (in)
          {
            case matrixpcaconst<T>::RPHI:
              return "RPHI";
            case matrixpcaconst<T>::RZ:
              return "RZ";
            case matrixpcaconst<T>::UNDEF:
              return "UNDEF";
          }
        
          return "";
        }

        static matrixpcaconst<T>::plane_type string_to_plane_type (
            const std::string & in)
        {
          if (in == "RPHI")
            return matrixpcaconst<T>::RPHI;
          else if (in == "RZ")
            return matrixpcaconst<T>::RZ;

          return matrixpcaconst<T>::UNDEF;
        }


        static std::string ttype_to_string (matrixpcaconst<T>::ttype in)
        {
          switch (in)
          {
            case matrixpcaconst<T>::FLOATPT:
              return "FLOATPT";
            case matrixpcaconst<T>::INTEGPT:
              return "INTEGPT";
            case matrixpcaconst<T>::NOONE:
              return "NOONE";
          }
        
          return "";
        }

        static matrixpcaconst<T>::ttype string_to_ttype (
            const std::string & in)
        {
          if (in == "FLOATPT")
            return matrixpcaconst<T>::FLOATPT;
          else if (in == "INTEGPT")
            return matrixpcaconst<T>::INTEGPT;

          return matrixpcaconst<T>::NOONE;
        }
 
  };

  template<typename T>
  matrixpcaconst<T>::matrixpcaconst (unsigned rows, unsigned cols, 
      const T* elements)
  {
    rows_ = rows;
    cols_ = cols;

    ttype_ = NOONE;
    ctype_ = NONDEF;
    stype_ = UNKN;
    ptype_ = UNDEF;

    towerid_ = 0;
    ptmin_ = 0.0;
    ptmax_ = 0.0;
    etamin_ = 0.0;
    etamax_ = 0.0;
    chargesign_ = 0;
    layersids_ = "unspecified";
  
    /*
    if( rows_ == 0 || cols_ == 0 )
      throw std::range_error("0 0 matx");
    */

    if (elements)
      for (unsigned i=0; i<rows_*cols_; i++)
        this->elements_.push_back(elements[i]);
    else
      for (unsigned i=0; i<rows_*cols_; i++)
        this->elements_.push_back((T)0);
  };
  
  template<typename T>
  matrixpcaconst<T>::matrixpcaconst( const matrixpcaconst<T>& cp )
  {
    rows_ = cp.rows_;
    cols_ = cp.cols_; 
    elements_ = cp.elements_;

    ttype_ = cp.ttype_;
    ctype_ = cp.ctype_;
    stype_ = cp.stype_;
    ptype_ = cp.ptype_;

    ptmin_ = cp.ptmin_;
    ptmax_ = cp.ptmax_;
    etamin_ = cp.etamin_;
    etamax_ = cp.etamax_;
    towerid_ = cp.towerid_;
    chargesign_ = cp.chargesign_;

    layersids_ = cp.layersids_;
  }

  template<typename T>
  matrixpcaconst<T>::~matrixpcaconst()
  {
    clear();
  }

  template<typename T>
  bool matrixpcaconst<T>::compare_nonelements (const matrixpcaconst<T>& cp) const
  {
    return ((rows_ == cp.rows_ ) &&
        (cols_ == cp.cols_ ) && 
        (ttype_ == cp.ttype_ ) &&
        (ctype_ == cp.ctype_ ) &&
        (stype_ == cp.stype_ ) &&
        (ptype_ == cp.ptype_ ) &&
        (ptmin_ == cp.ptmin_ ) &&
        (ptmax_ == cp.ptmax_ ) &&
        (etamin_ == cp.etamin_ ) &&
        (etamax_ == cp.etamax_ ) &&
        (towerid_ == cp.towerid_ ) && 
        (chargesign_ == cp.chargesign_ ) && 
        (cols_ == cp.cols_ ) && 
        (rows_ == cp.rows_ ) &&  
        (layersids_ == cp.layersids_ ));
  }
  
  template<typename T>
  void matrixpcaconst<T>::clear() 
  {
    elements_.clear();
    rows_ = 0;
    cols_ = 0;

    ttype_= NOONE;
    ctype_ = NONDEF;
    stype_ = UNKN;
    ptype_ = UNDEF;

    chargesign_ = 0;
    towerid_ = 0;
    ptmin_ = 0.0;
    ptmax_ = 0.0;
    etamin_ = 0.0;
    etamax_ = 0.0;
    layersids_ = "unspecified";
  }
  
  template<typename T>
  void matrixpcaconst<T>::reset(unsigned rows, unsigned cols) 
  {
    rows_ = rows;
    cols_ = cols;

    this->elements_.clear();

    for (unsigned i=0; i<rows_*cols_; i++)
      this->elements_.push_back((T)0);
  }
  
  template<typename T>
  matrixpcaconst<T>& pca::matrixpcaconst<T>::operator=( 
      const pca::matrixpcaconst<T>& cp )
  {
    /*
    if (cp.rows_ != rows_ && cp.cols_ != cols_ )
       throw std::domain_error("matrix op= not of same order");
    */

    this->reset(cp.rows_, cp.cols_);

    for(unsigned int i=0; i<rows_*cols_; i++)
      this->elements_[i] = cp.elements_[i];

    this->ttype_ = cp.ttype_;
    this->ctype_ = cp.ctype_;
    this->stype_ = cp.stype_;
    this->ptype_ = cp.ptype_;
    
    this->ptmin_ = cp.ptmin_;
    this->ptmax_ = cp.ptmax_;
    this->etamin_ = cp.etamin_;
    this->etamax_ = cp.etamax_;
    this->towerid_ = cp.towerid_;
    this->chargesign_ = cp.chargesign_;

    this->layersids_ = cp.layersids_;
  
    return *this;
  }
  
  template<typename T>
  void matrixpcaconst<T>::range_check(unsigned i, unsigned j) const
  {
    if (rows_ <= i)
      throw std::range_error("matrix access row out of range");
  
    if (cols_ <= j)
      throw std::range_error("matrix access col out of range");
  }

  template<typename T> 
  bool read_pcacosnt_from_file (std::vector<matrixpcaconst<T> > & vct, 
              const char * filename)
  {
    if (std::ifstream(filename))
    {
      std::ifstream infile;
      infile.open(filename);

      int numof;
      infile >> numof;

      for (int i=0; i < numof; ++i)
      {
        matrixpcaconst<T> matt(0, 0);
        std::string instr, instr1;
        int rows, cols;

        infile >> instr;
        matt.set_const_type(matrixpcaconst<T>::string_to_const_type (instr));
        infile >> instr;
        matt.set_towerid(atoi(instr.c_str()));
        infile >> instr;
        matt.set_layersids(instr.c_str());
        infile >> instr;
        matt.set_sector_type(matrixpcaconst<T>::string_to_sector_type(instr));
        infile >> instr;
        matt.set_plane_type(matrixpcaconst<T>::string_to_plane_type(instr));
        infile >> instr;
        matt.set_ttype(matrixpcaconst<T>::string_to_ttype(instr));
        infile >> instr >> instr1;
        matt.set_ptrange(atof(instr.c_str()), atof(instr1.c_str()));
        infile >> instr >> instr1;
        matt.set_etarange(atof(instr.c_str()), atof(instr1.c_str()));
        infile >> instr;
        matt.set_chargesign(atoi(instr.c_str()));
        infile >> rows;
        infile >> cols;
        matt.reset(rows, cols);
        for (unsigned int i = 0; i<matt.n_rows(); ++i)
          for (unsigned int j = 0; j<matt.n_cols(); ++j)
            infile >> matt(i, j);

        vct.push_back(matt);
      }

      infile.close();

      return true;
    }

    return false;
  }

  template<typename T>
  bool write_pcaconst_to_file (const matrixpcaconst<T> & in, 
      const char * filename)
  {
    if (std::ifstream(filename))
    {
      std::vector<matrixpcaconst<T> > vct;
      if (read_pcacosnt_from_file (vct, filename))
      {
        bool copied = false;
        for (unsigned int i = 0; i != vct.size(); ++i)
        {
          if(vct[i].compare_nonelements(in))
          {
            vct[i] = in;
            copied = true;
          }
        }

        if (!copied)
          vct.push_back(in);

        std::ofstream outf;
        outf.open(filename);

        outf << vct.size() << std::endl;
        for (unsigned int i = 0; i != vct.size(); ++i)
        {
          double min, max;
          outf << matrixpcaconst<T>::const_type_to_string(vct[i].get_const_type()) 
            << std::endl;
          outf << vct[i].get_towerid() << std::endl;
          if (vct[i].get_layersids() == "")
            outf << "unspecified" << std::endl;
          else
            outf << vct[i].get_layersids() << std::endl;
          outf << matrixpcaconst<T>::sector_type_to_string(vct[i].get_sector_type()) 
            << std::endl;
          outf << matrixpcaconst<T>::plane_type_to_string(vct[i].get_plane_type()) 
            << std::endl;
          outf << matrixpcaconst<T>::ttype_to_string(vct[i].get_ttype()) 
            << std::endl;
          vct[i].get_ptrange(min, max);
          outf.precision(10);
          outf << std::scientific << min << " " << max << std::endl;
          vct[i].get_etarange(min, max);
          outf << std::scientific << min << " " << max << std::endl;
          outf << vct[i].get_chargesign () << std::endl;
          outf << vct[i].n_rows() << std::endl;
          outf << vct[i].n_cols() << std::endl;

          outf.precision(10);
          for (unsigned int j = 0; j<vct[i].n_rows(); ++j)
            for (unsigned int k = 0; k<vct[i].n_cols(); ++k)
              outf << std::scientific << vct[i](j, k) << std::endl;
        }

        outf.close();
      }
      
      return false;
    }
    else
    {
      double min, max;

      std::ofstream outf;
      outf.open(filename);
      
      outf << "1" << std::endl;
      outf << matrixpcaconst<T>::const_type_to_string(in.get_const_type()) 
        << std::endl;
      outf << in.get_towerid() << std::endl;
      if (in.get_layersids() == "")
        outf << "unspecified" << std::endl;
      else
        outf << in.get_layersids() << std::endl;
      outf << matrixpcaconst<T>::sector_type_to_string(in.get_sector_type()) 
        << std::endl;
      outf << matrixpcaconst<T>::plane_type_to_string(in.get_plane_type()) 
        << std::endl;
      outf << matrixpcaconst<T>::ttype_to_string(in.get_ttype()) 
        << std::endl;
      in.get_ptrange(min, max);
      outf.precision(10);
      outf << std::scientific << min << " " << max << std::endl;
      in.get_etarange(min, max);
      outf << std::scientific << min << " " << max << std::endl;
      outf << in.get_chargesign() << std::endl;
      outf << in.n_rows() << std::endl;
      outf << in.n_cols() << std::endl;

      outf.precision(10);
      for (unsigned int i = 0; i<in.n_rows(); ++i)
        for (unsigned int j = 0; j<in.n_cols(); ++j)
          outf << std::scientific << in(i, j) << std::endl;
      outf.close();
    }

     return true;
  }

}

#endif
