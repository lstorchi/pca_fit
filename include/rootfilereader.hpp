#ifndef _PCAROOTFILEREADER_INC_ 
#define _PCAROOTFILEREADER_INC_

#include <armadillo>
#include <vector>
#include <map>

// lstorchi: this should be moved into include as soon as the utility code should be
//           a class... This is not a proper class, just a way to collect some functions 
//           as soon as have time please rewrite it and move main_computation here ... 

namespace pca
{
  class rootfilereader
  {
    public:
      rootfilereader();
      ~rootfilereader();
  
      std::string get_filename();
      void set_filename(const std::string &);
      
      const std::string & get_errmsg () const;
      int get_errnum() const;
  
    private:
  
      int errnum_;
      std::string errmsg_, filename_;

      void set_errmsg (int , const std::string &);
      void reset_error ();
  };
};

#endif
