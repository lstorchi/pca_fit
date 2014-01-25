/******
 ==== The structure of the code ====

	 gf_process
             |
	     ---- gf_fep
	             |
                     ----gf_fep_unpack
                              |
			      ---- gf_iword_decode
                         (dump) dump_gf_ievt
	                  gf_fep_comb
			 (dump) dump_fep_out
	          gf_fit
	             |
                     ----gf_mkaddr
                         gf_fit_proc                         
                         gf_fit_format
                         (dump) dump_fit

			 
                  gf_comparator
                    |  		    
		    -----gf_chi2
                         gf_gfunc
		         gf_formatter
			 (dump)dump_fout
                     

****/


#include <linux/types.h>
#include <string.h>
#include <asm/msr.h>

//#include "main.h"
//#include "net.h"
#include <libgffit/gffit.h>

//#define DEBUG_SVT
//#define DUMP
//#define DEBUG_READFCON

static char svtsim_err_str[10][256];
static int svtsim_n_err;
 
 /*
  * Wrap memory allocation to allow debug, etc.
  */
 
 typedef struct memDebug_s {
   struct memDebug_s *last, *next;
   const char *fileName;
   int lineNum, userBytes;
   int deadBeef;
 } memDebug_t;

 static int n_alloc = 0, t_alloc = 0, m_alloc = 0;
 static memDebug_t *chain_alloc = 0;

#if 0
 static void *badaddr[] = {  /* to hunt down problematic mallocs */
 };
#endif
 
 
/*
 * Wrap memory allocation for SVT simulators
 */
int svtsim_memDebug(int chainDump);
void *svtsim_malloc1(size_t size, const char *filename, int linenum);

void svtsim_free1(void *tmp, const char *filename, int linenum);
void *svtsim_realloc1(void *inptr, size_t size, 
		      const char *filename, int linenum);
#define svtsim_malloc(x) (svtsim_malloc1((x), __FILE__, __LINE__))
#define svtsim_free(x) (svtsim_free1((x), __FILE__, __LINE__))
#define svtsim_mallocForever(x) (svtsim_malloc1((x), __FILE__, -__LINE__))
#define svtsim_realloc(x, y) (svtsim_realloc1((x), (y), __FILE__, __LINE__))
#define svtsim_assert(x) \
  do { \
  if ((x)) continue; \
  svtsim_assert_set( __FILE__, __LINE__ ); \
  } while (0)

#define svtsim_assert(x) \
  do { \
  if ((x)) continue; \
  svtsim_assert_set( __FILE__, __LINE__ ); \
  } while (0)
 

/*
 * Useful macros for parsing input files
 */
#define Get svtsim_assert(fgets(str, sizeof(str), fp))
#define ExpectStr(S) \
  svtsim_assert(fgets(str, sizeof(str), fp) && \
		1==sscanf(str, "%s", s) && !strcmp(s, (S)))
#define ExpectInt(N) \
  svtsim_assert(fgets(str, sizeof(str), fp) && \
		1==sscanf(str, "%d", &n) && n==(N))
#define GetReal(X) do { \
  svtsim_assert(fgets(str, sizeof(str), fp) && \
		1==sscanf(str, "%f", &x)); (X) = x; } while (0)
#define GetInt(X) do { \
  svtsim_assert(fgets(str, sizeof(str), fp) && \
		1==sscanf(str, "%d", &n)); (X) = n; } while (0)
#define GetHex(X) do { \
  svtsim_assert(fgets(str, sizeof(str), fp) && \
		1==sscanf(str, "%x", &n)); (X) = n; } while (0)

#define MAXMKA 8192
