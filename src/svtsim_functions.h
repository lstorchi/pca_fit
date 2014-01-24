#ifndef _SVTSIM_H_
#define _SVTSIM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <cuda_runtime.h>

typedef struct tf_arrays      *tf_arrays_t;

#define SVTSIM_MIN(a, b) ((a)<(b)?(a):(b))
#define SVTSIM_MAX(a, b) ((a)>(b)?(a):(b)) 
#define SVTSIM_NEL(x) (sizeof((x))/sizeof((x)[0]))

#define GFS_ERR_BIT 11

#define OCVR_LSB   10 /* was 11. put to 10 to match TF - OCVR_LSB+1 = LSB of curvature in GF output word */


#define GFS_OFL_HIT 5
#define GFS_OFL_LYR 6
#define GFS_OFL_CMB 7
#define GFS_OFL_CHI 7
#define GFS_INV_DATA 8
#define GFS_OFL_FIT 9
#define GFS_FIFO_ERR 10
#define GFS_ERR_BIT 11


#define GF_TOTCHI2_CUTVAL 0x00000f
#define GF_ERRMASK_SVT 0x0007ff
#define GF_ERRMASK_CDF 0x0007ff
#define GF_ERRMASK_EOE 0x000000


//GF
#define NEVTS 500
//#define NEVTS 1700
//#define NROADS 100
#define NROADS 64
#define EVT 0 //all words from the input file are considered belonging to the same event

/* SVTSIM PARAMETERS */

#define MAXCOE_VSIZE 256 
#define OFF_SUBA_WIDTH 3
#define OFF_SUBA_LSB 19

#define NSHIFTS 3
#define NCHI 3
#define CHI_DWIDTH 10

enum {
  SVTSIM_NWEDGE = 12,/* number of wedges (in phi) */
  SVTSIM_NBAR = 6,/* number of half-barrels (in z) */
  SVTSIM_NLAY = 6,/* number of layers of hits in SVT */
  SVTSIM_NLAYSVX = 5,/* number of silicon layers */
  SVTSIM_NPL = 7,/* number of physical layers (L00+SVX+XFT) */
  SVTSIM_XFTPHIBINS = 2304/* number of XFT phi bins (incl mini-phi) */
};
   
#define PAR_ADDR_WIDTH 5
#define SVTSIM_GF_MKADDR_INVALID 2

  
//from mapset file 512kpatt_20100921000000_20100921104017
//xtfAlgo         1
//xtfP0           8  
//xtfWildVersion0
#define PHIOFFSETBINS 12
//xtfTracksByWedge 1
//xtfTrackDropMode 1
//#patternCenter  x(um)|y(um)|dx(urad)|dy(urad)|beamspot(um)
//patternCenter   -1578|1340|470|-140|1400
  
  
#define MINHITS 5

//#define DEBUG_READFCON 1

enum {
  SVTSIM_GF_OK, 
  SVTSIM_GF_ERR
};
  

#define SVTSIM_NPL 7/* number of physical layers (L00+SVX+XFT) */
#define DIMSPA  6
#define NFITPAR 6
#define WEDGE 6
#define SVTSIM_NBAR  6/* number of half-barrels (in z) */
#define UNKNOWN_ERR 11
#define FITBLOCK 32

#define SVT_WORD_WIDTH 23 /* length of word in input from HB to GF */
#define SVT_SECT_WIDTH 4
#define SVT_SECT_LSB 17
#define SVT_ERR_LSB   9
#define EE_WORD_WIDTH 21 
#define AMROAD_WORD_WIDTH 23


#define SVTNHITS 7

#define NTFWORDS 7

#define MAXROAD 64
#define NSVX_PLANE 5
#define MAX_HIT 7
/* was MAXCOMB 100 - try to run with TrigMon */
#define MAXCOMB 64
#define MAXCOMB5H 5
#define MINLYR 5

#define HIT_PHI 5
#define HIT_CRV 6


/*the following variables are positions in gf_maskdata2[] array  */
#define WIDTH_15B          15 /* 0x7fff */   
#define WIDTH_14B          14 /* 0x3fff */ 
#define WIDTH_15B_NEG_SET     17 /* 0x4000 */
#define WIDTH_6B_NEG_SET     18 /* 0x20 */


#define WIDTH_32B          32 /* 0x3fff */


#define GF_HIT_WIDTH   8  /* position in gf_mask of mask for HIT width  */
#define GF_PHI_WIDTH   12
#define GF_CRV_WIDTH    7
#define GF_ROAD_WIDTH  21

#define SVT_EP_BIT      21
#define SVT_EE_BIT      22
#define SVT_PAR_BIT     8

#define SVT_CRV_WIDTH   6
#define SVT_CRV_LSB     12
#define SVT_D_WIDTH     11
#define SVT_PHI_WIDTH   12
#define SVT_CRV_SHIFT   6
#define SVT_D_SHIFT     0 /* da rimettere a 6 (FC) */
#define SVT_SECT_WIDTH  4
#define SVT_TRKID_LSB   3
#define SVT_TRKID_WIDTH 9
#define SVT_ROAD_WIDTH  21
#define SVT_LYR_LSB     18
#define SVT_LYR_WIDTH   3
#define SVT_Z_LSB       15
#define SVT_Z_WIDTH     3
#define SVT_LCLS_BIT    14
#define SVT_HIT_WIDTH   15

#define MADDR_HITM_LSB 8
#define MADDR_NCLS_LSB 3

#define NFITTER 6 /* number of fit to be done (3 paramters + 3 chisq) */

#define OPHI_WIDTH 13 /* position in gf_mask of mask for PHI width */
#define OCVR_WIDTH 8  /* position in gf_mask of mask for CVR width */
#define OCVR_SIGN  8  /* sign bit of curvature */
#define OIMP_WIDTH 9  /* position in gf_mask of mask for IMP width */
#define OIMP_SIGN  9  /* sign bit of impact parameter */
#define OCHI_WIDTH 14  /* position in gf_mask of mask for CHI width */
#define OAMROAD_WIDTH  17 /* position in gf_mask of mask for AMroad width */
#define OSEC_LSB       17 /* (OSEC_LSB + 1)= LSB of sector in GF output word */ 
#define CHI2SUM_WIDTH 11
#define OCHI2_LSB 10
#define OSTAT_LSB 9 
#define OBP_ERR_BIT 19
#define OBP_ID_BIT 20
#define OX1_LSB 10
#define OX3_LSB 10


#define GF_ZID_WIDTH  6
#define GF_STAT_WIDTH 12
#define GF_SUBZ_WIDTH 3

#define OFLOW_HIT_BIT 1
#define UFLOW_HIT_BIT 2
#define OUTORDER_BIT 3
#define OFLOW_CMB_BIT 4
#define OFLOW_CHI_BIT 9
#define INT_OFLOW_BIT 4


#define PARITY_ERR_BIT 0
#define INV_DATA_BIT 3

/* Layer Definition */
#define XFT_LYR 5
#define XFT_LYRID 7
#define EP_LYR  8
#define EE_LYR 24  

#define FIT_DWIDTH 48
#define FIT_PHI_DWIDTH 48 /* was 22 in TF */
#define FIT_CRV_DWIDTH 48 /* was 20 in TF */


#define FIT_RESULT_OFLOW_BIT 10

#define MAXCHI2A 64
#define CHI_AWIDTH 6 /*  in gf firmware is 1 */

/*
 * Useful macros for parsing input files
 */



#define gf_mask(x) (gf_maskdata[(x)])
static int gf_maskdata[] = {
  0x00000000,
  0x00000001, 0x00000003, 0x00000007, 0x0000000f,
  0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff,
  0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff,
  0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff,
  0x0001ffff, 0x0003ffff, 0x0007ffff, 0x000fffff,
  0x001fffff, 0x003fffff, 0x007fffff, 0x00ffffff,
  0x01ffffff, 0x03ffffff, 0x07ffffff, 0x0fffffff,
  0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff
};


#define gf_mask3(x) (gf_maskdata3[(x)])
static long long int gf_maskdata3[] = {
  0x000000000000ULL,
  0x000000000001ULL, 0x000000000003ULL, 0x000000000007ULL, 0x00000000000fULL,
  0x00000000001fULL, 0x00000000003fULL, 0x00000000007fULL, 0x0000000000ffULL,
  0x0000000001ffULL, 0x0000000003ffULL, 0x0000000007ffULL, 0x000000000fffULL,
  0x000000001fffULL, 0x000000003fffULL, 0x000000007fffULL, 0x00000000ffffULL,
  0x00000001ffffULL, 0x00000003ffffULL, 0x00000007ffffULL, 0x0000000fffffULL,
  0x0000001fffffULL, 0x0000003fffffULL, 0x0000007fffffULL, 0x000000ffffffULL,
  0x000001ffffffULL, 0x000003ffffffULL, 0x000007ffffffULL, 0x00000fffffffULL,
  0x00001fffffffULL, 0x00003fffffffULL, 0x00007fffffffULL, 0x0000ffffffffULL,
  0x0001ffffffffULL, 0x0003ffffffffULL, 0x0007ffffffffULL, 0x000fffffffffULL,
  0x001fffffffffULL, 0x003fffffffffULL, 0x007fffffffffULL, 0x00ffffffffffULL,
  0x01ffffffffffULL, 0x03ffffffffffULL, 0x07ffffffffffULL, 0x0fffffffffffULL,
  0x1fffffffffffULL, 0x3fffffffffffULL, 0x7fffffffffffULL, 0xffffffffffffULL 
};




typedef struct svtsim_cable_s {
  unsigned int *data;          /* pointer to cable data */ 
  int ndata;                    /* number of input words */ 
  int mdata;                    /* number of allocated words */ 
} svtsim_cable_t;


struct fep_arrays {

  int fep_ncmb[NEVTS][MAXROAD];
  int fep_hit[NEVTS][MAXROAD][MAXCOMB][NSVX_PLANE];
  int fep_phi[NEVTS][MAXROAD][MAXCOMB];
  int fep_crv[NEVTS][MAXROAD][MAXCOMB];
  int fep_lcl[NEVTS][MAXROAD][MAXCOMB];
  int fep_lclforcut[NEVTS][MAXROAD][MAXCOMB];
  int fep_hitmap[NEVTS][MAXROAD][MAXCOMB];
  int fep_zid[NEVTS][MAXROAD];
  int fep_road[NEVTS][MAXROAD];
  int fep_sect[NEVTS][MAXROAD];
  int fep_cable_sect[NEVTS][MAXROAD];
  int fep_err[NEVTS][MAXROAD][MAXCOMB][MAXCOMB5H];
  int fep_crv_sign[NEVTS][MAXROAD][MAXCOMB];
  int fep_ncomb5h[NEVTS][MAXROAD][MAXCOMB];
  int fep_hitZ[NEVTS][MAXROAD][MAXCOMB][NSVX_PLANE];
  int fep_nroads[NEVTS];
  int fep_ee_word[NEVTS];
  int fep_err_sum[NEVTS];

};

struct extra_data {

  int totEvts;
  int wedge[NEVTS];
  int whichFit[SVTSIM_NBAR][FITBLOCK];/* handles TF mkaddr degeneracy */
  int evt_road[NEVTS][NROADS];
  long long int lfitparfcon[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK]; /* full-precision coeffs, P0s as read from fcon files SA*/

};

struct fit_arrays {

  long long int fit_fit[NEVTS][6][MAXROAD][MAXCOMB][MAXCOMB5H];
  int fit_err[NEVTS][MAXROAD][MAXCOMB][MAXCOMB5H];
  int fit_err_sum[NEVTS];

};

struct fout_arrays {

  int fout_parity[NEVTS];
  int fout_ntrks[NEVTS];
//  int fout_iroad[NEVTS][MAXROAD*MAXCOMB];
//  int fout_icmb[NEVTS][MAXROAD*MAXCOMB];
  unsigned int fout_gfword[NEVTS][MAXROAD*MAXCOMB][NTFWORDS];
  int fout_cdferr[NEVTS];
  int fout_svterr[NEVTS];
  int fout_ee_word[NEVTS];
  int fout_err_sum[NEVTS];

};

struct tf_arrays {

  int totEvts;

  int evt_hit[NEVTS][NROADS][NSVX_PLANE][MAX_HIT];
  int evt_hitZ[NEVTS][NROADS][NSVX_PLANE][MAX_HIT];
  int evt_lcl[NEVTS][NROADS][NSVX_PLANE][MAX_HIT];
  int evt_lclforcut[NEVTS][NROADS][NSVX_PLANE][MAX_HIT];
  int evt_layerZ[NEVTS][NROADS][NSVX_PLANE];
  int evt_zid[NEVTS][NROADS];
  int evt_nhits[NEVTS][NROADS][NSVX_PLANE+1];
  int evt_err[NEVTS][NROADS];
  int evt_crv[NEVTS][NROADS][MAX_HIT];
  int evt_crv_sign[NEVTS][NROADS][MAX_HIT];
  int evt_phi[NEVTS][NROADS][MAX_HIT];
  int evt_err_sum[NEVTS];
  int evt_cable_sect[NEVTS][NROADS];
  int evt_sect[NEVTS][NROADS];
  int evt_nroads[NEVTS];
  int evt_road[NEVTS][NROADS];
  int evt_ee_word[NEVTS];
  
  int fep_ncmb[NEVTS][MAXROAD];
  int fep_hit[NEVTS][MAXROAD][MAXCOMB][NSVX_PLANE];
  int fep_phi[NEVTS][MAXROAD][MAXCOMB];
  int fep_crv[NEVTS][MAXROAD][MAXCOMB];
  int fep_lcl[NEVTS][MAXROAD][MAXCOMB];
  int fep_lclforcut[NEVTS][MAXROAD][MAXCOMB];
  int fep_hitmap[NEVTS][MAXROAD][MAXCOMB];
  int fep_zid[NEVTS][MAXROAD];
  int fep_road[NEVTS][MAXROAD];
  int fep_sect[NEVTS][MAXROAD];
  int fep_cable_sect[NEVTS][MAXROAD];
  int fep_err[NEVTS][MAXROAD][MAXCOMB][MAXCOMB5H];
  int fep_crv_sign[NEVTS][MAXROAD][MAXCOMB];
  int fep_ncomb5h[NEVTS][MAXROAD][MAXCOMB];
  int fep_hitZ[NEVTS][MAXROAD][MAXCOMB][NSVX_PLANE];
  int fep_nroads[NEVTS];
  int fep_ee_word[NEVTS];
  int fep_err_sum[NEVTS];
  
  long long int fit_fit[NEVTS][6][MAXROAD][MAXCOMB][MAXCOMB5H];
  int fit_err[NEVTS][MAXROAD][MAXCOMB][MAXCOMB5H];
  int fit_err_sum[NEVTS];
 
  // formatter_out 
  int fout_parity[NEVTS];
  int fout_ntrks[NEVTS]; 
  int fout_iroad[NEVTS][MAXROAD*MAXCOMB];
  int fout_icmb[NEVTS][MAXROAD*MAXCOMB];
  unsigned int fout_gfword[NEVTS][MAXROAD*MAXCOMB][NTFWORDS];
  int fout_cdferr[NEVTS];
  int fout_svterr[NEVTS];
  int fout_ee_word[NEVTS];
  int fout_err_sum[NEVTS];
  int fout_found[NEVTS];

 
 //output "cable"
  svtsim_cable_t *out;
  
  //gf_memory
  int wedge[NEVTS];
  unsigned int *mem_mkaddr[NEVTS];
  short mem_coeff[NFITTER][SVTNHITS][MAXCOE_VSIZE]; 
  int mem_nintcp;
  short (*intcp)[NFITTER];
  int mem_fitsft[NFITTER][NSHIFTS]; 
  /* int chi2[NCHI][MAXCHI2A]; */
  int minhits; /* The minimum number of hits that we require 
		  (including XFT hit)*/
  int chi2cut;
  int svt_emsk;
  int gf_emsk;
  int cdf_emsk;
  int eoe_emsk; /* MASK for the errors */

  int whichFit[SVTSIM_NBAR][FITBLOCK];/* handles TF mkaddr degeneracy */
  int ifitpar[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK];       /* TF coefficients, P0s */
  long long int lfitpar[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK];/* full-precision coeffs, P0s */
  long long int lfitparfcon[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK]; /* full-precision coeffs, P0s as read from fcon files SA*/
  float gcon[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK];     /* [fdc012][I0123PC][z][whichfit] */
  int xftphi_shiftbits[NFITPAR];/* TF bit shifts */
  int xftcrv_shiftbits[NFITPAR];
  int result_shiftbits[NFITPAR];
  int lMap[SVTSIM_NBAR][SVTSIM_NPL];/* map physical layer => cable layer */
  int oX[SVTSIM_NBAR], oY[SVTSIM_NBAR];/* origin used for fcon */
  float phiUnit, dvxUnit, crvUnit;/* units for fit parameters */
  float k0Unit, k1Unit, k2Unit;/* units for constraints */
  int dphiNumer, dphiDenom;/* dphi(wedge) = 2pi*N/D */
  int mkaddrBogusValue;

  
  
};


//int svtsim_fconread(tf_arrays_t tf);

int svtsim_fconread_GPU(tf_arrays_t tf);

//int gf_fit(tf_arrays_t tf);

//int gf_fep_CPU(tf_arrays_t tf, int n_words_in, unsigned int* data_in); 

//int gf_fep(tf_arrays_t tf, int n_words_in, void* data); 

int gf_fep_unpack(tf_arrays_t tf, int n_words_in, void* data);

int gf_fep_unpack_CPU(tf_arrays_t tf, int n_words_in, unsigned int* data_in);

int  gf_fep_comb(tf_arrays_t tf);

//int gf_init(tf_arrays_t* ptr_tf);


//void svtsim_cable_copywords(svtsim_cable_t *cable, unsigned int *word, int nword);

//void svtsim_cable_addwords(svtsim_cable_t *cable, unsigned int *word, int nword);

//void svtsim_cable_addword(svtsim_cable_t *cable, unsigned int word);

int  gf_mkaddr(tf_arrays_t tf,
	       int hitmap, int lclmap, int zmap, 
	       int *coe_addr, int *int_addr, int *addr, int *err);

int  svtsim_get_gfMkAddr(tf_arrays_t tf, int *d, int nd, int d0);


int gf_fit_proc(int hit[], int sign_crv, long long int coeff[], long long int intcp, long long int *result, int *err);

//int gf_fit_format(tf_arrays_t tf);

int gf_chi2(tf_arrays_t tf, long long int chi[], int * trk_err, long long int *chi2);

int gf_gfunc(int ncomb5h, int icomb5h, int hitmap, int lcmap, int chi2);

int gf_getq(int layer_config);

int gf_formatter(tf_arrays_t tf,int ie, int ir, int ic, int ich,int chi2); 

void dump_fout(tf_arrays_t tf); 

static int gf_formatter_err(int err_sum, int cdf_emsk, int svt_emsk, int eoe_emsk , int * eoe_err, int * cdf_err, int * svt_err);


int gf_stword(int id, int err, int mask);

//int gf_comparator(tf_arrays_t tf);

int whichFit[SVTSIM_NBAR][FITBLOCK];/* handles TF mkaddr degeneracy */
int ifitpar[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK];       /* TF coefficients, P0s */
long long int lfitpar[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK];/* full-precision coeffs, P0s */
long long int lfitparfcon[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK]; /* full-precision coeffs, P0s as read from fcon files SA*/
float gcon[NFITPAR][DIMSPA+1][SVTSIM_NBAR][FITBLOCK];     /* [fdc012][I0123PC][z][whichfit] */
int xftphi_shiftbits[NFITPAR];/* TF bit shifts */
int xftcrv_shiftbits[NFITPAR];
int result_shiftbits[NFITPAR];
int lMap[SVTSIM_NBAR][SVTSIM_NPL];/* map physical layer => cable layer */
int oX[SVTSIM_NBAR], oY[SVTSIM_NBAR];/* origin used for fcon */
float phiUnit, dvxUnit, crvUnit;/* units for fit parameters */
float k0Unit, k1Unit, k2Unit;/* units for constraints */
int dphiNumer, dphiDenom;/* dphi(wedge) = 2pi*N/D */
int mkaddrBogusValue;

#endif
