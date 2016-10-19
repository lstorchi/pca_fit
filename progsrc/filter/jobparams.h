#ifndef jobparams_H
#define jobparams_H

#include <string>
#include <cstdio>
#include <tclap/CmdLine.h>
using namespace TCLAP;

class jobparams{

 public:

  /** constructer reading from standard arg */
  jobparams(int argc, char** argv);

  /** default constructer */
  jobparams(){}

  /** copy constructer */
  jobparams(const jobparams& tparams); 

  /** destructer */
  ~jobparams(){}

  /** return value */

  std::string option() const;
  std::string inputfile() const;
  std::string sector() const;
  std::string outfile() const;
  int         nsec() const;
  int         lim() const;

 private:

  std::string  m_opt;  
  std::string  m_inputfile; 
  std::string  m_outfile;
  std::string  m_sector; 
  int          m_nsec;
  int          m_lim;
};

inline std::string jobparams::option() const{
  return m_opt;
}

inline std::string jobparams::inputfile() const{
  return m_inputfile;
}

inline std::string jobparams::outfile() const{
  return m_outfile;
}

inline std::string jobparams::sector() const{
  return m_sector;
}

inline int jobparams::nsec() const{
  return m_nsec;
}

inline int jobparams::lim() const{
  return m_lim;
}

#endif
