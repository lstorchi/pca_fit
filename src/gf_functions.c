#include "private_includes.h"

int
gf_chi2 (tf_arrays_t tf, long long int chi[], int *trk_err,
         long long int *chi2)
{
  int i;
  int oflow = 0;
  long long int temp = 0;
  long long int chi2memdata = 0;
  *chi2 = 0;

  /* --------- Executable starts here ------------ */
  for (i = 0; i < NCHI; i++)

  {

#ifdef DEBUG_SVT
    printf ("chi[%d]: %.6x \n", i, chi[i]);

#endif /*   */
    temp = abs (chi[i]);
    if (chi[i] < 0)
      temp++;

    /* 
       if (chi[i]>=MAXCHI2A)  
       oflow = 1; 

       temp = (temp & gf_mask(CHI_AWIDTH));       
     */

#ifdef DEBUG_SVT
    printf ("temp: %.6x \n", temp);

#endif /*   */
    chi2memdata = temp * temp;
    *chi2 += chi2memdata;

#ifdef DEBUG_SVT
    printf ("chi2memdata: %.6llx, *chi2: %.6llx \n", chi2memdata, *chi2);

#endif /*   */
  }
  *chi2 = (*chi2 >> 2);
  if ((*chi2 >> 2) > gf_mask (CHI_DWIDTH))

  {

    /*      *chi2 = (*chi2 & gf_mask(CHI_DWIDTH)) + (1<<CHI_DWIDTH); */
    *chi2 = 0x7ff;
    *trk_err |= (1 << OFLOW_CHI_BIT);
  }

#ifdef DEBUG_SVT
  printf ("*chi2: %.6x \n", *chi2);

#endif /*   */
  return SVTSIM_GF_OK;
}

int
gf_getq (int lyr_config)
{
  int q = 0;
  switch (lyr_config)
  {
  case 0x01e:                  /* lcmap = 00000, hitmap = 11110 */
    q = 3;
    break;
  case 0x01d:                  /* lcmap = 00000, hitmap = 11101 */
    q = 2;
    break;
  case 0x01b:                  /* lcmap = 00000, hitmap = 11011 */
    q = 1;
    break;
  case 0x017:                  /* lcmap = 00000, hitmap = 10111 */
    q = 2;
    break;
  case 0x00f:                  /* lcmap = 00000, hitmap = 01111 */
    q = 2;
    break;
  case 0x03e:                  /* lcmap = 00001, hitmap = 11110 */
    q = 2;
    break;
  case 0x03d:                  /* lcmap = 00001, hitmap = 11101 */
    q = 1;
    break;
  case 0x03b:                  /* lcmap = 00001, hitmap = 11011 */
    q = 1;
    break;
  case 0x037:                  /* lcmap = 00001, hitmap = 10111 */
    q = 1;
    break;
  case 0x02f:                  /* lcmap = 00001, hitmap = 01111 */
    q = 1;
    break;
  case 0x05e:                  /* lcmap = 00010, hitmap = 11110 */
    q = 7;
    break;
  case 0x05d:                  /* lcmap = 00010, hitmap = 11101 */
    q = 1;
    break;
  case 0x05b:                  /* lcmap = 00010, hitmap = 11011 */
    q = 2;
    break;
  case 0x057:                  /* lcmap = 00010, hitmap = 10111 */
    q = 2;
    break;
  case 0x04f:                  /* lcmap = 00010, hitmap = 01111 */
    q = 2;
    break;
  case 0x09e:                  /* lcmap = 00100, hitmap = 11110 */
    q = 7;
    break;
  case 0x09d:                  /* lcmap = 00100, hitmap = 11101 */
    q = 2;
    break;
  case 0x09b:                  /* lcmap = 00100, hitmap = 11011 */
    q = 1;
    break;
  case 0x097:                  /* lcmap = 00100, hitmap = 10111 */
    q = 2;
    break;
  case 0x08f:                  /* lcmap = 00100, hitmap = 01111 */
    q = 3;
    break;
  case 0x11e:                  /* lcmap = 01000, hitmap = 11110 */
    q = 7;
    break;
  case 0x11d:                  /* lcmap = 01000, hitmap = 11101 */
    q = 2;
    break;
  case 0x11b:                  /* lcmap = 01000, hitmap = 11011 */
    q = 2;
    break;
  case 0x117:                  /* lcmap = 01000, hitmap = 10111 */
    q = 1;
    break;
  case 0x10f:                  /* lcmap = 01000, hitmap = 01111 */
    q = 3;
    break;
  case 0x21e:                  /* lcmap = 10000, hitmap = 11110 */
    q = 7;
    break;
  case 0x21d:                  /* lcmap = 10000, hitmap = 11101 */
    q = 2;
    break;
  case 0x21b:                  /* lcmap = 10000, hitmap = 11011 */
    q = 2;
    break;
  case 0x217:                  /* lcmap = 10000, hitmap = 10111 */
    q = 2;
    break;
  case 0x20f:                  /* lcmap = 10000, hitmap = 01111 */
    q = 1;
    break;
  case 0x0de:                  /* lcmap = 00110, hitmap = 11110 */
    q = 7;
    break;
  case 0x0dd:                  /* lcmap = 00110, hitmap = 11101 */
    q = 1;
    break;
  case 0x0db:                  /* lcmap = 00110, hitmap = 11011 */
    q = 2;
    break;
  case 0x0d7:                  /* lcmap = 00110, hitmap = 10111 */
    q = 3;
    break;
  case 0x0cf:                  /* lcmap = 00110, hitmap = 01111 */
    q = 4;
    break;
  case 0x19e:                  /* lcmap = 01100, hitmap = 11110 */
    q = 7;
    break;
  case 0x19d:                  /* lcmap = 01100, hitmap = 11101 */
    q = 2;
    break;
  case 0x19b:                  /* lcmap = 01100, hitmap = 11011 */
    q = 1;
    break;
  case 0x197:                  /* lcmap = 01100, hitmap = 10111 */
    q = 1;
    break;
  case 0x18f:                  /* lcmap = 01100, hitmap = 01111 */
    q = 3;
    break;
  case 0x31e:                  /* lcmap = 11000, hitmap = 11110 */
    q = 7;
    break;
  case 0x31d:                  /* lcmap = 11000, hitmap = 11101 */
    q = 3;
    break;
  case 0x31b:                  /* lcmap = 11000, hitmap = 11011 */
    q = 3;
    break;
  case 0x317:                  /* lcmap = 11000, hitmap = 10111 */
    q = 1;
    break;
  case 0x30f:                  /* lcmap = 11000, hitmap = 01111 */
    q = 2;
    break;
  case 0x15e:                  /* lcmap = 01010, hitmap = 11110 */
    q = 7;
    break;
  case 0x15d:                  /* lcmap = 01010, hitmap = 11101 */
    q = 1;
    break;
  case 0x15b:                  /* lcmap = 01010, hitmap = 11011 */
    q = 3;
    break;
  case 0x157:                  /* lcmap = 01010, hitmap = 10111 */
    q = 2;
    break;
  case 0x14f:                  /* lcmap = 01010, hitmap = 01111 */
    q = 4;
    break;
  case 0x25e:                  /* lcmap = 10010, hitmap = 11110 */
    q = 7;
    break;
  case 0x25d:                  /* lcmap = 10010, hitmap = 11101 */
    q = 1;
    break;
  case 0x25b:                  /* lcmap = 10010, hitmap = 11011 */
    q = 2;
    break;
  case 0x257:                  /* lcmap = 10010, hitmap = 10111 */
    q = 2;
    break;
  case 0x24f:                  /* lcmap = 10010, hitmap = 01111 */
    q = 1;
    break;
  case 0x29e:                  /* lcmap = 10100, hitmap = 11110 */
    q = 7;
    break;
  case 0x29d:                  /* lcmap = 10100, hitmap = 11101 */
    q = 2;
    break;
  case 0x29b:                  /* lcmap = 10100, hitmap = 11011 */
    q = 1;
    break;
  case 0x297:                  /* lcmap = 10100, hitmap = 10111 */
    q = 2;
    break;
  case 0x28f:                  /* lcmap = 10100, hitmap = 01111 */
    q = 1;
    break;
  default:
    q = 7;
    break;
  }
  return q;
}

int
gf_gfunc (int ncomb5h, int icomb5h, int hitmap, int lcmap, int chi2)
{
  int lyr_config;
  int gvalue;
  int newhitmap = 0;
  int newlcmap = 0;
  int q = 0;

#ifdef DEBUG_SVT
  printf ("in gf_gfunc: ncomb5h = %d, icomb5h = %d\n", ncomb5h, icomb5h);

#endif /*   */
  if (ncomb5h == 1)
  {
    newhitmap = hitmap;
    newlcmap = lcmap;
  }

  else if (ncomb5h == 5)
  {
    switch (icomb5h)
    {
    case 0:                    /*  11110 */
      newhitmap = 0x1e;
      newlcmap = (lcmap & 0x1e);
      break;
    case 1:                    /*  11101 */
      newhitmap = 0x1d;
      newlcmap = lcmap & 0x1d;
      break;
    case 2:                    /*  11011 */
      newhitmap = 0x1b;
      newlcmap = lcmap & 0x1b;
      break;
    case 3:                    /*  10111 */
      newhitmap = 0x17;
      newlcmap = lcmap & 0x17;
      break;
    case 4:                    /*  01111 */
      newhitmap = 0x0f;
      newlcmap = lcmap & 0x0f;
      break;
    }
  }

  else
  {
    printf ("ERROR in gf_gfunc: wrong number of combinations\n");
  }

  /* end  if(ncomb5h == xx) */
#ifdef DEBUG_SVT
  printf ("in gf_gfunc: newhitmap = %.6x, newlcmap = %.6x, chi2 = %.6x\n",
          newhitmap, newlcmap, chi2);

#endif /*   */
  lyr_config = newhitmap + (newlcmap << 5);

#ifdef DEBUG_SVT
  printf ("layer configuration: %.6x\n", lyr_config);

#endif /*   */
  q = gf_getq (lyr_config);

#ifdef DEBUG_SVT
  printf ("quality value: %d\n", q);
  printf ("q << 4  : %.6x\n", q << 4);
  printf ("chi2 >> 7: %.6x\n", chi2 >> 7);

#endif /*   */
  gvalue = (q << 4) + ((chi2 & 0x3ff) >> 6);
  return gvalue;
}

int
gf_formatter_err (int err, int cdfmsk, int svtmsk, int eoemsk, int *eoe,
                  int *cdf, int *svt)
     /* 
        Simulate the board error conditions (CDF-ERR, SVT-ERR and EOE-ERR) 
        INPUT: err; error summary. 
        cdfmsk; Mask for the CDF-ERR. 
        svtmsk; Mask for the SVT-ERR. 
        eoemsk; Mask for the EOE-ERR. 
        OUTPUT: *eoe; EOE error 
        *cdf; CDF error 
        *svt; SVT error 
      */
{
  int i;

  /* --------- Executable starts here ------------ */
  *cdf = 0;                     /* never turned ON except for the FIFO overflow */
  *svt = 0;
  *eoe = 0;
  for (i = 0; i <= FIT_RESULT_OFLOW_BIT; i++)
  {
    if ((err >> i) & gf_mask (1))
    {
      if (((svtmsk >> i) & gf_mask (1)) == 0)
        *svt = 1;
      if (i == 0)
      {
        if (((eoemsk >> PARITY_ERR_BIT) & gf_mask (1)) == 0)
        {
          *eoe |= (1 << PARITY_ERR_BIT);
        }
      }
      else if ((i == 2) || (i == 3))
      {
        if (((eoemsk >> INV_DATA_BIT) & gf_mask (1)) == 0)
        {
          *eoe |= (1 << INV_DATA_BIT);
        }
      }
      else
      {
        if (((eoemsk >> INT_OFLOW_BIT) & gf_mask (1)) == 0)
        {
          *eoe |= (1 << INT_OFLOW_BIT);
        }
      }
    }                           /* if ((err>>i)&gf_mask(1)) { */
  }                             /* for (i=0; i<= FIT_RESULT_OFLOW_BIT; i++) { */
  return SVTSIM_GF_OK;
}

int
gf_stword (int id, int err, int mask)
     /* 
        Compose the GF status word in the 7th word from the GF  
        INPUT : err; error summary 
        mask; MASK 
        OUTPUT : return the gf_stword 

        NOTE: Currently this code does not support the parity error and 
        FIFO error. 
      */
{
  int word;
  word = id;
  if ((err >> OFLOW_HIT_BIT) & gf_mask (1))
    word |= (1 << GFS_OFL_HIT);
  if ((err >> OFLOW_CHI_BIT) & gf_mask (1))
    word |= (1 << GFS_OFL_CHI);
  if (((err >> UFLOW_HIT_BIT) & gf_mask (1)) ||
      ((err >> OUTORDER_BIT) & gf_mask (1)))
    word |= (1 << GFS_INV_DATA);

  /*  
     if ((err>>FIT_RESULT_OFLOW_BIT)&gf_mask(1)) 
     word |= (1<<GFS_OFL_FIT); 
   */
  return word;
}

int
gf_formatter (tf_arrays_t tf, int ie, int ir, int ic, int ich, int chi2)
{

  /* 
     Pack the GF data output from the (raw) fit values. 
     INPUT: hits;  
     fit; FIT output 
     ir; the road ID  (=-1 and ic=-1 means the EE word.) 
     ic; the combination ID 
     chi2; the Chi^2 value. 
     OUTPUT: *form_out;  
   */
  int i, it, err, obp_err, eoe_err, tmp_phi;
  int hit_form[NSVX_PLANE];
  int icomb5h = 0;              /* number of 4/5 tracks out of a 5/5 track */
  int z = 0;                    /* z should be 6 bits large */
  int sector = 0;
  int AMroad = 0x6789a;
  int gf_stat = 0;
  int trackID = 1;

  /* --------- Executable starts here ------------ */

#ifdef DEBUG_SVT
  printf ("in gf_formatter: ie: %d, ir: %d, ic: %d\n", ie, ir, ic);
  printf ("in gf_formatter: gf->fep->ee_word:%.6x\n", tf->fep_ee_word[ie]);

#endif /*   */
  if ((ir < 0) && (ic < 0))
  {
    tf->fout_err_sum[ie] = (tf->fep_err_sum[ie] | tf->fit_err_sum[ie]);
    gf_formatter_err (tf->fout_err_sum[ie], GF_ERRMASK_CDF,
                      GF_ERRMASK_SVT, GF_ERRMASK_EOE, &eoe_err,
                      &tf->fout_cdferr[ie], &tf->fout_svterr[ie]);
    tf->fout_ee_word[ie] =
      (tf->fep_ee_word[ie] &
       (gf_mask (SVT_WORD_WIDTH) & ~(1 << SVT_PAR_BIT)));
    tf->fout_ee_word[ie] |= (eoe_err << SVT_ERR_LSB);
    tf->fout_ee_word[ie] |= (tf->fout_parity[ie] << SVT_PAR_BIT);

#ifdef DEBUG_SVT
    printf
      ("in gf_formatter(%d): fout->parity = %.6x, eoe_err = %.6x, ee_word = %.6x\n",
       ie, tf->fout_parity[ie], eoe_err, tf->fout_ee_word[ie]);

#endif /*   */

    /*    fout->ee_word = (gf->fep->ee_word &  
       (gf_mask(SVT_WORD_WIDTH))); */
    return SVTSIM_GF_OK;
  }
  it = tf->fout_ntrks[ie];

#ifdef DEBUG_SVT
  printf ("in gf_formatter: formatting word %d\n", it);

#endif /*   */
  (tf->fout_ntrks[ie])++;
  tf->fout_iroad[ie][it] = tf->fep_road[ie][ir];
  tf->fout_icmb[ie][it] = ic + 1;
  err = (tf->fep_err[ie][ir][ic][ich] | tf->fit_err[ie][ir][ic][ich]);

#ifdef DEBUG_SVT
  printf ("Formatting SVX hits for output words\n");
  printf ("EE word: fout->ee_word:%.6x\n", tf->fout_ee_word[ie]);

#endif /*   */
  for (i = 0; i < NSVX_PLANE; i++)
  {

    /* Hit coordinate */
#ifdef DEBUG_SVT
    printf ("i = %d, ich = %d\n", i, ich);

#endif /*   */
    if (tf->fep_ncomb5h[ie][ir][ic] == 5)
    {
      if (i != ich)
      {
        hit_form[i] = tf->fep_hit[ie][ir][ic][i] & gf_mask (GF_HIT_WIDTH);

        /* Long Cluster bit */
        hit_form[i] +=
          (((tf->fep_hit[ie][ir][ic][i] & 0x4000) ? 1 : 0) << GF_HIT_WIDTH);

        /* Hit existence bit */
        hit_form[i] +=
          (((tf->fep_hitmap[ie][ir][ic] >> i) & gf_mask (1)) <<
           (GF_HIT_WIDTH + 1));
        hit_form[i] = (hit_form[i] & gf_mask (GF_HIT_WIDTH + 2));
      }

      else
      {
        hit_form[i] = 0;
      }
    }

    else
    {
      hit_form[i] = tf->fep_hit[ie][ir][ic][i] & gf_mask (GF_HIT_WIDTH);

      /* Long Cluster bit */
      hit_form[i] +=
        (((tf->fep_hit[ie][ir][ic][i] & 0x4000) ? 1 : 0) << GF_HIT_WIDTH);

      /* Hit existence bit */
      hit_form[i] +=
        (((tf->fep_hitmap[ie][ir][ic] >> i) & gf_mask (1)) <<
         (GF_HIT_WIDTH + 1));
      hit_form[i] = (hit_form[i] & gf_mask (GF_HIT_WIDTH + 2));
    }

#ifdef DEBUG_SVT
    printf ("hit[%d]: %.6x\n", i, tf->fep_hit[ie][ir][ic][i]);
    printf ("hit_form[%d]: %.6x\n", i, hit_form[i]);

#endif /*   */
  }
  if (1)
  {

    /* 
     * Compute hit status 
     */
    int presentmask = 0;
    int newhitmap = 0;
    if (tf->fep_ncomb5h[ie][ir][ic] == 1)
    {
      presentmask = tf->fep_hitmap[ie][ir][ic];
    }

    else if (tf->fep_ncomb5h[ie][ir][ic] == 5)
    {
      switch (ich)
      {
      case 0:                  /*  11110 */
        newhitmap = 0x1e;
        break;
      case 1:                  /*  11101 */
        newhitmap = 0x1d;
        break;
      case 2:                  /*  11011 */
        newhitmap = 0x1b;
        break;
      case 3:                  /*  10111 */
        newhitmap = 0x17;
        break;
      case 4:                  /*  01111 */
        newhitmap = 0x0f;
        break;
      }
      presentmask = newhitmap;
    }

    {
      int longmask = presentmask & tf->fep_lcl[ie][ir][ic];
      int goodmask = presentmask & ~longmask;
      int badmask = 0x1f & ~goodmask;
      static int badmap[] = { 0x0,      /* 00000: all layers good */
        0x5,                    /* 10000: layer 0 bad */
        0x4,                    /* 01000: layer 1 bad */
        0xe,                    /* 11000: layers 0,1 bad  (changed from f to e) */
        0x3,                    /* 00100: layer 2 bad */
        0xe,                    /* 10100: layers 0,2 bad */
        0xb,                    /* 01100: layers 1,2 bad */
        0xf,                    /* 11100: >2 layers bad */
        0x2,                    /* 00010: layer 3 bad */
        0xd,                    /* 10010: layers 0,3 bad */
        0xa,                    /* 01010: layers 1,3 bad */
        0xf,                    /* 11010: >2 layers bad */
        0x8,                    /* 00110: layers 2,3 bad */
        0xf,                    /* 10110: >2 layers bad */
        0xf,                    /* 01110: >2 layers bad */
        0xf,                    /* 11110: >2 layers bad */
        0x1,                    /* 00001: layer 4 bad */
        0xc,                    /* 10001: layers 0,4 bad */
        0x8,                    /* 01001: layers 1,4 bad  (oops: doc says 0x9 not 0x8) */
        0xf,                    /* 11001: >2 layers bad */
        0x7,                    /* 00101: layers 2,4 bad */
        0xf,                    /* 10101: >2 layers bad */
        0xf,                    /* 01101: >2 layers bad */
        0xf,                    /* 11101: >2 layers bad */
        0x6,                    /* 00011: layers 3,4 bad */
        0xf,                    /* 10011: >2 layers bad */
        0xf,                    /* 01011: >2 layers bad */
        0xf,                    /* 11011: >2 layers bad */
        0xf,                    /* 00111: >2 layers bad */
        0xf,                    /* 10111: >2 layers bad */
        0xf,                    /* 01111: >2 layers bad */
        0xf                     /* 11111: all layers bad! */
      };
      gf_stat = badmap[badmask];
  }}

#ifdef DEBUG_SVT
  printf ("gfstat before gf_stword: %.6x\n", gf_stat);

#endif /*   */
  gf_stat = gf_stword (gf_stat, err, tf->gf_emsk);

  /* obp_err is the MSB of the status word, Error Summary = OR of the  
     errors contained in the status word 
   */
  obp_err = (gf_stat >> GFS_ERR_BIT) & gf_mask (1);

  /* output word (25 bits) (from CDFnote 5026) 
     4-3-2-1-0-9-8-7-6-5-4-3-2-1-0-9-8-7-6-5-4-3-2-1-0                 
   */
#ifdef DEBUG_SVT
  printf ("\n OUTPUT WORDS \n");

#endif /*   */
  /* 1st word  
     24-23-22-21- 20- 19- 18-17-16-15-14-13- 12-11-10-9-8-7-6-5-4-3-2-1-0  
     --------     1   -  z                  phi      
   */

  /* phi is already formatted by the fitter (13 bits) */
  if (tf->fep_ncomb5h[ie][ir][ic] == 1)
  {
    z = tf->fep_zid[ie][ir];
  }

  else if (tf->fep_ncomb5h[ie][ir][ic] == 5)
  {

#ifdef DEBUG_SVT
    printf ("gf->fep->hitZ[%d][%d][%d][%d] = %d\n", ie, ir, ic, ich,
            tf->fep_hitZ[ie][ir][ic][ich]);

#endif /*   */
    if (ich == 0)
    {
      z =
        ((tf->fep_hitZ[ie][ir][ic][4] & gf_mask (GF_SUBZ_WIDTH)) <<
         GF_SUBZ_WIDTH) +
        (tf->fep_hitZ[ie][ir][ic][1] & gf_mask (GF_SUBZ_WIDTH));
    }

    else if (ich == 4)
    {
      z =
        ((tf->fep_hitZ[ie][ir][ic][3] & gf_mask (GF_SUBZ_WIDTH)) <<
         GF_SUBZ_WIDTH) +
        (tf->fep_hitZ[ie][ir][ic][0] & gf_mask (GF_SUBZ_WIDTH));
    }

    else
    {
      z =
        ((tf->fep_hitZ[ie][ir][ic][4] & gf_mask (GF_SUBZ_WIDTH)) <<
         GF_SUBZ_WIDTH) +
        (tf->fep_hitZ[ie][ir][ic][0] & gf_mask (GF_SUBZ_WIDTH));
    }
  }
  tf->fout_gfword[ie][it][0] =
    (tf->fit_fit[ie][0][ir][ic][ich] & gf_mask (OPHI_WIDTH)) + ((z &
                                                                 gf_mask
                                                                 (GF_ZID_WIDTH))
                                                                << OPHI_WIDTH)
    /*   + (obp_err << OBP_ERR_BIT)  */
    + (0 << OBP_ERR_BIT)
    /* we follow the word structure in  http://www-cdf.fnal.gov/internal/upgrades/daq_trig/trigger/svt/BoardDocs/data_words/tracks_bits.html */
    + (1 << (OBP_ID_BIT));

#ifdef DEBUG_SVT
  printf ("phi: %.6x\n", tf->fit_fit[ie][0][ir][ic][ich]);
  printf ("zin: %d, zout:%d, z=%.6x\n", z & gf_mask (GF_SUBZ_WIDTH),
          (z >> GF_SUBZ_WIDTH) & gf_mask (GF_SUBZ_WIDTH), z);
  printf ("obp_err: %d\n", obp_err);
  printf ("1st word out of formatter: %.7x\n", tf->fout_gfword[ie][it][0]);
  printf ("\n\n");

#endif /*   */

  /* 2nd word  
     4-3-2-1-0-9-8   -7-6-5-4-3-2-1-0 -9   -8-7-6-5-4-3-2-1-0  
     24-23-22-21- 20- 19-  18-  17-16-15-14-13- 12-11-  10-9-8-7-6-5-4-3-2-1-0  
     ------------  rID      sign c                       d 
     17mo bit di roadID -> 19 
     18mo               -> 20 
   */
  tf->fout_gfword[ie][it][1] = tf->fit_fit[ie][1][ir][ic][ich]
    + (tf->fit_fit[ie][2][ir][ic][ich] << OCVR_LSB)
    + ((tf->evt_road[ie][ir] & 0x60000) << 2);

#ifdef DEBUG_SVT
  printf ("imp w sign: %.6x\n", tf->fit_fit[ie][1][ir][ic][ich]);
  printf ("cvr w sign: %.6x\n", tf->fit_fit[ie][2][ir][ic][ich]);
  printf ("cvr w sign << OCVR_LSB: %.6x\n",
          tf->fit_fit[ie][2][ic][ic][ich] << OCVR_LSB);
  printf ("AMroad : %.6x\n", AMroad);
  printf ("(AMroad & 0x30000) : %.6x\n", AMroad & 0x30000);
  printf ("2nd word out of formatter: %.12x\n", tf->fout_gfword[ie][it][1]);
  printf ("\n\n");

#endif /*   */

  /* 3rd word  
     4-3-2-1-0-9-8-7 -6-5-4-3-2-1-0-9-8-7-6-5-4-3-2-1-0  
     --------sector   AM road id (17 LSB) 
   */
  tf->fout_gfword[ie][it][2] =
    (tf->evt_road[ie][ir] & gf_mask (OAMROAD_WIDTH)) +
    ((tf->fep_cable_sect[ie][ir] & gf_mask (SVT_SECT_WIDTH)) << OSEC_LSB);

#ifdef DEBUG_SVT
  printf ("AMroad & gf_mask(OAMROAD_WIDTH): %.12x\n",
          tf->evt_road[ie][ir] & gf_mask (OAMROAD_WIDTH));
  printf ("sector & gf_mask(SVT_SECT_WIDTH):%.12x\n",
          tf->fep_cable_sect[ie][ir] & gf_mask (SVT_SECT_WIDTH));
  printf ("sector & gf_mask(SVT_SECT_WIDTH))<< OSEC_LSB):%.12x\n",
          (tf->fep_cable_sect[ie][ir] & gf_mask (SVT_SECT_WIDTH)) <<
          OSEC_LSB);
  printf ("3rd word out of formatter: %.12x\n", tf->fout_gfword[ie][it][2]);
  printf ("\n\n");

#endif /*   */

  /* 4th word  
     4-3-2-1-0-9-8-7-6-5-4-3-2-1-0 -9-8-7-6-5-4-3-2-1-0  
     ----------x1                   x0 
     bit 21 = bit 19 del roadID 
     hit = 8 bassi e 2 alti      
   */
  tf->fout_gfword[ie][it][3] = hit_form[0] + (hit_form[1] << OX1_LSB)
    + ((tf->evt_road[ie][ir] & 0x80000) << 1);

#ifdef DEBUG_SVT
  printf ("x0: %.12x\n", hit_form[0]);
  printf ("x1: %.12x\n", hit_form[1]);
  printf ("x1 << OX1_LSB : %.12x\n", hit_form[1] << OX1_LSB);
  printf ("AMroad & 0x40000: %.12x\n", tf->evt_road[ie][ir] & 0x40000);
  printf ("AMroad : %.6x\n", AMroad);
  printf ("4th word out of formatter: %.12x\n", tf->fout_gfword[ie][it][3]);
  printf ("\n\n");

#endif /*   */

  /* 5th word  
     4-3-2-1-0-9-8-7-6-5-4-3-2-1-0 -9-8-7-6-5-4-3-2-1-0  
     ----------x3                   x2 
     bit 21 = road ID 20 
   */
  tf->fout_gfword[ie][it][4] = hit_form[2] + (hit_form[3] << OX3_LSB)
    + ((tf->evt_road[ie][ir] & 0x100000));

#ifdef DEBUG_SVT
  printf ("x2: %.12x\n", hit_form[2]);
  printf ("x3: %.12x\n", hit_form[3]);
  printf ("x3 << OX3_LSB : %.12x\n", hit_form[2] << OX3_LSB);
  printf ("5th word out of formatter: %.12x\n", tf->fout_gfword[ie][it][4]);
  printf ("\n\n");

#endif /*   */

  /* 6th word  
     4-3-2-1-0-9-8-7-6-5-4-3-2-1-0 -9-8-7-6-5-4-3-2-1-0  
     ----------chisq                x4 
   */
  tf->fout_gfword[ie][it][5] = hit_form[4]
    + ((chi2 & gf_mask (CHI2SUM_WIDTH)) << OCHI2_LSB);

#ifdef DEBUG_SVT
  printf ("x4: %.12x\n", tf->fep_hit[ie][ir][ic][4]);
  printf ("chiSq: %.12x\n", chi2);
  printf ("chiSq & gf_mask(CHI2SUM_WIDTH): %.12x\n",
          chi2 & gf_mask (CHI2SUM_WIDTH));
  printf ("chiSq & gf_mask(CHI2SUM_WIDTH)) << OCHI2_LSB): %.12x\n",
          ((chi2 & gf_mask (CHI2SUM_WIDTH)) << OCHI2_LSB));
  printf ("6th word out of formatter: %.12x\n", tf->fout_gfword[ie][it][5]);
  printf ("\n\n");

#endif /*   */

  /* 7th word  
     4-3-2-1 -0-9-8-7-6-5-4-3-2-1-0-9 -8-7-6-5-4-3-2-1-0  
     ------0  TrackFitter status       Track Number                 
     Track Num = identificativo della traccia XFT 
     phi - 3 bit meno significativi del phi della traccia XFT 
   */
  tf->fout_gfword[ie][it][6] = ((tf->fep_phi[ie][ir][ic] >> SVT_TRKID_LSB)
                                & gf_mask (SVT_TRKID_WIDTH))
    + ((gf_stat & gf_mask (GF_STAT_WIDTH)) << OSTAT_LSB) + (1 << SVT_EP_BIT);

#ifdef DEBUG_SVT
  printf ("track Number: %.12x\n", trackID);
  printf ("GigaFitter status: %.12x\n", gf_stat);
  printf ("gf_stat & gf_mask(GF_STAT_WIDTH): %.12x\n",
          gf_stat & gf_mask (GF_STAT_WIDTH));
  printf ("gf_stat & gf_mask(GF_STAT_WIDTH)) << OSTAT_LSB: %.12x\n",
          ((gf_stat & gf_mask (GF_STAT_WIDTH)) << OSTAT_LSB));
  printf ("7th word out of formatter: %.12x\n", tf->fout_gfword[ie][it][6]);
  printf ("\n\n");

#endif /*   */
  for (i = 0; i < NTFWORDS; i++)
    tf->fout_parity[ie] ^= cal_parity (tf->fout_gfword[ie][it][i]);

#ifdef DEBUG_SVT
  printf ("in gf_formatter (end): fout->parity = %.6x\n",
          tf->fout_parity[ie]);

#endif /*   */
  return SVTSIM_GF_OK;
}


int
gf_comparator (tf_arrays_t tf)
{

  /* 

     Pack the GF data output from the (raw) fit values. 
     INPUT: hits;  
     fit; FIT output 
     ir; the road ID  (=-1 and ic=-1 means the EE word.) 
     ic; the combination ID 
     chi2; the Chi^2 value. 
     OUTPUT: *form_out;  
   */
  int ie;
  int i, it, err, obp_err, eoe_err, tmp_phi;
  long long int chi[3], chi2;
  int chipass = 0;

  /* int ChiSqCut = 0xfffffff; *//*No ChiSq cut */
  /* int ChiSqCut = 0x0; *//*No tracks in output, only EE */
  /* int ChiSqCut = 0xf40; */
  int ChiSqCut;
  int ir, ic, ich = 0;
  int gvalue;
  int gvalue_best;
  int ind_best = 0;
  int chi2_best = 0;
  int gvalue_cut = 0x70;
  int bestTrackFound = 0;

  /* --------- Executable starts here ------------ */
  svtsim_cable_copywords (tf->out, 0, 0);
  for (ie = 0; ie < tf->totEvts; ie++)
  {

    /*  ChiSqCut = mem->par.chi2cut;  */
    ChiSqCut = 0x40;
    tf->fout_ntrks[ie] = 0;
    tf->fout_parity[ie] = 0;

#ifdef DEBUG_SVT
    printf ("\n\n *********** COMPARATOR ************* \n\n");
    printf ("chi2 cut = %.6x\n", ChiSqCut);
    printf ("tf->fep_nroads[%d] = %d\n", ie, tf->fep_nroads[ie]);

#endif /*   */
    for (ir = 0; ir < tf->fep_nroads[ie]; ir++)
    {                           /* Loop for all roads */
      for (ic = 0; ic < tf->fep_ncmb[ie][ir]; ic++)
      {                         /* Loop for all combinations */
        gvalue_best = 0x70;

#ifdef DEBUG_SVT
        printf
          ("-------------> Processing road %d (ID = %.6x), combination %d\n",
           ir, tf->evt_road[ie][ir], ic);

#endif /*   */
        if (tf->fep_ncomb5h[ie][ir][ic] == 1)
        {

#ifdef DEBUG_SVT
          printf ("In gf_comparator, 4/5 track\n");
          printf ("Calculating ChiSq\n");

#endif /*   */
          for (i = 0; i < NCHI; i++)
          {
            chi[i] = tf->fit_fit[ie][i + 3][ir][ic][0];

#ifdef DEBUG_SVT
            printf ("ir = %d, ic = %d, i = %d, chi[i] = %.6x \n",
                    ir, ic, i, chi[i]);

#endif /*   */
          }

          /*  calculate chisq */
          if (gf_chi2 (tf, chi, &tf->fit_err[ie][ir][ic][0], &chi2)
              == SVTSIM_GF_ERR)
          {

            //    svtsim_gf_err(gf, "gf_chi2 has an error in tf_opp"); 
            return SVTSIM_GF_ERR;
          }

          /* check chiSq and select the track having the best chiSq */
#ifdef DEBUG_SVT
          printf ("chi2 in comparator: %.6x, %d \n", chi2, chi2);
          printf ("ChiSqCut in comparator: %.6x \n", ChiSqCut);

#endif /*   */
          if (chi2 <= ChiSqCut)
          {
            chi2 = chi2 >> 2;   /* FC - hack .. see matching shift in gf_chi2 */

#ifdef DEBUG_SVT
            printf ("ChiSq: %.6x, ir: %d, ic: %d \n", chi2, ir, ic);
            printf ("fep->ncomb5h[%d][%d][%d] = %d\n", ie, ir, ic,
                    tf->fep_ncomb5h[ie][ir][ic]);
            printf ("Good chi2 for this 4/5 track...calling gfunc \n");
            printf ("hitmap = %.6x, lcmap = %.6x\n",
                    tf->fep_hitmap[ie][ir][ic], tf->fep_lcl[ie][ir][ic]);

#endif /*   */

            gvalue = 
              gf_gfunc (tf->fep_ncomb5h[ie][ir][ic], ich,
                        tf->fep_hitmap[ie][ir][ic],
                        tf->fep_lcl[ie][ir][ic],
                        (chi2 & gf_mask (CHI2SUM_WIDTH)));

#ifdef DEBUG_SVT
            printf ("g_func value = %.6x\n", gvalue);

#endif /*   */
            if (gvalue < gvalue_cut)
            {

#ifdef DEBUG_SVT
              printf ("Good gfunc for this 4/5 track...calling formatter \n");

#endif /*   */
              if (gf_formatter (tf, ie, ir, ic, 0, chi2) == SVTSIM_GF_ERR)
              {

                //            svtsim_gf_err(gf, "gf_comparator has an error in gf_formatter"); 
                return SVTSIM_GF_ERR;
              }
              svtsim_cable_addwords (tf->out,
                                     tf->fout_gfword[ie]
                                     [tf->fout_ntrks[ie] - 1], NTFWORDS);
            }                   /* end  if(gvalue <  gvalue_cut) { */
          }
        }

        else if (tf->fep_ncomb5h[ie][ir][ic] == 5)
        {

#ifdef DEBUG_SVT
          printf
            ("In gf_comparator, 5/5 track. Start loop on 4/5 sub-tracks\n");

#endif /*   */
          bestTrackFound = 0;
          gvalue_best = 999;
          ind_best = 999;
          chi2_best = 999;
          for (ich = 0; ich < tf->fep_ncomb5h[ie][ir][ic]; ich++)
          {

#ifdef DEBUG_SVT
            printf ("5/5 track: sub track # %d\n", ich);
            printf ("Calculating ChiSq\n");

#endif /*   */
            for (i = 0; i < NCHI; i++)
            {
              chi[i] = tf->fit_fit[ie][i + 3][ir][ic][ich];

#ifdef DEBUG_SVT
              printf
                ("ie = %d, ir = %d, ic = %d, i = %d, chi[i] = %.6x \n",
                 ie, ir, ic, i, chi[i]);

#endif /*   */
            }

            /*  calculate chisq */
            if (gf_chi2
                (tf, chi, &tf->fit_err[ie][ir][ic][ich],
                 &chi2) == SVTSIM_GF_ERR)
            {

              //            svtsim_gf_err(gf, "gf_chi2 has an error in tf_opp"); 
              return SVTSIM_GF_ERR;
            }

            /* check chiSq  */
#ifdef DEBUG_SVT
            printf ("chi2 in comparator: %.6x, %d \n", chi2, chi2);
            printf ("ChiSqCut in comparator: %.6x \n", ChiSqCut);

#endif /*   */
            if (chi2 <= ChiSqCut)
            {
              chi2 = chi2 >> 2; /* FC - hack .. see matching shift in gf_chi2 */

#ifdef DEBUG_SVT
              printf ("ChiSq: %.6x, ir: %d, ic: %d \n", chi2, ir, ic);
              printf ("fep->ncomb5h[%d][%d] = %d\n", ir, ic,
                      tf->fep_ncomb5h[ie][ir][ic]);
              printf
                ("Good chi2 for track # %d of 5/5 track...calling gfunc\n",
                 ich);
              printf ("hitmap = %.6x, lcmap = %.6x\n",
                      tf->fep_hitmap[ie][ir][ic], tf->fep_lcl[ie][ir][ic]);

#endif /*   */
              gvalue =
                gf_gfunc (tf->fep_ncomb5h[ie][ir][ic], ich,
                          tf->fep_hitmap[ie][ir][ic],
                          tf->fep_lcl[ie][ir][ic],
                          (chi2 & gf_mask (CHI2SUM_WIDTH)));

#ifdef DEBUG_SVT
              printf ("g_func value for track %d = %.6x\n", ich, gvalue);

#endif /*   */
              if ((gvalue < gvalue_cut) && (gvalue < gvalue_best))
              {
                gvalue_best = gvalue;
                ind_best = ich;
                chi2_best = chi2;
                bestTrackFound = 1;
              }
            }                   /*  end if(chi2 <= ChiSqCut) { */
          }                     /*  end for(ich = 0; ich < gf->fep->ncomb5h[ir][ic]; ich++){ */
          if (bestTrackFound)
          {

#ifdef DEBUG_SVT
            printf
              ("gf_gfunc calculation completed, now calling gf_formatter for the best track %d\n",
               ind_best);
            printf ("\n ****************** Formatter *************** \n");

#endif /*   */
            if (gf_formatter (tf, ie, ir, ic, ind_best, chi2_best)
                == SVTSIM_GF_ERR)
            {

              //            svtsim_gf_err(gf, "gf_comparator has an error in gf_formatter"); 
              return SVTSIM_GF_ERR;
            }

            // #ifdef DEBUG_SVT 
            //printf("calling svtsim_cable_addwords for track %d",tf->fout_ntrks[ie]-1); 
            //printf("\n ****************** Formatter *************** \n"); 
            // #endif            
            svtsim_cable_addwords (tf->out,
                                   tf->fout_gfword[ie][tf->fout_ntrks
                                                       [ie] - 1], NTFWORDS);
          }                     /* end if(bestTrackFound) { */
        }                       /* end  if(gf->fep->ncomb5h[ir][ic] == 1) { */
      }                         /* end of for (ic=0; ic<gf->fep->ncmb[ir]; ic++) - loop on combinations */
    }                           /* end of for (ir=0; ir<gf->fep->nroads; ir++) - loop on roads */

    /* End of Events */
    ir = -1;
    ic = -1;

    //    printf("in gf_comparator: event # %d\n", ie); 
    gf_formatter (tf, ie, ir, ic, ich, chi2);
    svtsim_cable_addword (tf->out, tf->fout_ee_word[ie]);
  }                             //end loop on events 
#ifdef DUMP
  dump_fout (tf);

#endif /*   */
  return 0;
}

/* 
 * Simulate the performance of the FEP processor on the board. 
 * 
 */
int
gf_fep (tf_arrays_t tf, int n_words_in, void *data)
{
  int ir, ic;                   /* Loop variables for the roads and the combinations */
  int status;
  int gfwd[4];

  //  struct fep_out *fep = tf->fep; 

  /* --------- Executable starts here ------------ */

  /* Unpack the input data */
  gf_fep_unpack (tf, n_words_in, data);

  /* dump the unpack result if requested */
#ifdef DUMP
  dump_gf_ievt (tf);

#endif /*   */
  /* Try all possible combinations of track candidates */
  gf_fep_comb (tf);

  /* dump the combination result if requested */
#ifdef DUMP
  dump_fep_out (tf);

#endif /*   */

  return 0;
}

  /* 
   * Simulate the performance of the FEP processor on the board. 
   * NOTE: see the CDF note xxxx for the details 
   */
int
gf_fep_CPU (tf_arrays_t tf, int n_words_in, unsigned int *data_in)
{
  int ir, ic;                   /* Loop variables for the roads and the combinations */
  int status;
  int gfwd[4];

  //  struct fep_out *fep = tf->fep; 

  /* --------- Executable starts here ------------ */
  printf ("in gf_fep_CPU\n");
  int i;

  //   for(i = 0; i < n_words_in; i++) { 
  // printf("data[%d] = %.6x\n", i, data_in[i]);  
  //} 
  /* Unpack the input data */
  gf_fep_unpack_CPU (tf, n_words_in, data_in);

  /* dump the unpack result if requested */
#ifdef DUMP
  dump_gf_ievt (tf);

#endif /*   */
  /* Try all possible combinations of track candidates */
  gf_fep_comb (tf);

  /* dump the combination result if requested */
#ifdef DUMP
  dump_fep_out (tf);

#endif /*   */

  return 0;
}

  /*  
   * Unpack the input words into the FEP  
   */
int
gf_fep_unpack (tf_arrays_t tf, int n_words_in, void *data_in)
{
  int id, gf_xft, id_last = -1;
  int ee_word, parity_in, sector, amroad, crv, crv_sign, phi, zid, lcl, hit,
    err_sum;
  int evt, rd, ie;
  evt = EVT;
  gf_xft = 0;
  tf->totEvts = 0;

  /* initializing arrays */
  for (ie = 0; ie < NEVTS; ie++)
  {
    tf->evt_nroads[ie] = 0;
    tf->evt_err_sum[ie] = 0;
    for (id = 0; id <= NSVX_PLANE; id++)
    {
      tf->evt_layerZ[ie][tf->evt_nroads[ie]][id] = 0;
    }
    for (id = 0; id <= XFT_LYR; id++)
    {
      tf->evt_nhits[ie][tf->evt_nroads[ie]][id] = 0;
    }
    tf->evt_err[ie][tf->evt_nroads[ie]] = 0;
    tf->evt_zid[ie][tf->evt_nroads[ie]] = -1;

    //    printf("tf->evt_nroads[%d] = %d, tf->evt_zid[%d][tf->evt_nroads[%d]] = %d\n", ie, tf->evt_nroads[ie], ie, ie, tf->evt_zid[ie][tf->evt_nroads[ie]]); 
  }

  //GF Now unpack the data 

  //   unsigned int* data = (unsigned int*) (n->rx_bufs[n->my_rank].in); 
  unsigned int *data = (unsigned int *) data_in;
  int i;
  for (i = 0; i < n_words_in; i++)
  {
    id =
      word_decode (data[i], gf_xft, &ee_word, &parity_in, &sector, &amroad,
                   &crv, &crv_sign, &phi, &zid, &lcl, &hit, &err_sum);

#ifdef DEBUG_SVT
    printf ("in fep_unpack, after word_decode: word %d, id = %d\n", i, id);

#endif /*   */
    /* SVX Data */
    if (id < XFT_LYR)
    {
      tf->evt_hit[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                [tf->evt_nroads[evt]][id]] =
        hit;
      tf->evt_hitZ[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                 [tf->evt_nroads[evt]][id]] =
        zid;
      tf->evt_lcl[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                [tf->evt_nroads[evt]][id]] =
        lcl;
      tf->evt_lclforcut[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                      [tf->evt_nroads[evt]]
                                                      [id]] = lcl;
      tf->evt_layerZ[evt][tf->evt_nroads[evt]][id] = zid;
      if (tf->evt_zid[evt][tf->evt_nroads[evt]] == -1)
      {
        tf->evt_zid[evt][tf->evt_nroads[evt]] = zid & gf_mask (GF_SUBZ_WIDTH);
      }
      else
      {
        tf->evt_zid[evt][tf->evt_nroads[evt]] =
          (((zid & gf_mask (GF_SUBZ_WIDTH)) << GF_SUBZ_WIDTH)
           +
           (tf->evt_zid[evt][tf->evt_nroads[evt]] & gf_mask (GF_SUBZ_WIDTH)));
      }
      tf->evt_nhits[evt][tf->evt_nroads[evt]][id]++;

      /* Error Checking */
      if (tf->evt_nhits[evt][tf->evt_nroads[evt]][id] == MAX_HIT)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OFLOW_HIT_BIT);
      if (id < id_last)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OUTORDER_BIT);
    }
    else if ((id == XFT_LYR) && (gf_xft == 0))
    {
      gf_xft = 1;
    }
    else if ((id == XFT_LYR) && (gf_xft == 1))
    {
      gf_xft = 0;
      tf->evt_crv[evt][tf->evt_nroads[evt]][tf->evt_nhits[evt][tf->evt_nroads
                                                               [evt]][id]] =
        crv;
      tf->evt_crv_sign[evt][tf->evt_nroads[evt]][tf->evt_nhits[evt]
                                                 [tf->evt_nroads[evt]][id]] =
        crv_sign;
      tf->evt_phi[evt][tf->evt_nroads[evt]][tf->evt_nhits[evt][tf->evt_nroads
                                                               [evt]][id]] =
        phi;
      tf->evt_nhits[evt][tf->evt_nroads[evt]][id]++;

      /* Error Checking */
      if (tf->evt_nhits[evt][tf->evt_nroads[evt]][id] == MAX_HIT)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OFLOW_HIT_BIT);
      if (id < id_last)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OUTORDER_BIT);
    }
    else if (id == EP_LYR)
    {
      tf->evt_cable_sect[evt][tf->evt_nroads[evt]] = sector;
      tf->evt_sect[evt][tf->evt_nroads[evt]] = sector;
      tf->evt_road[evt][tf->evt_nroads[evt]] = amroad;
      tf->evt_err_sum[evt] |= tf->evt_err[evt][tf->evt_nroads[evt]];
      tf->evt_nroads[evt]++;
      if (tf->evt_nroads[evt] > MAXROAD)
      {
        printf
          ("The limit on the number of roads fitted by the TF is %d\n",
           MAXROAD);
        printf ("You reached that limit evt->nroads = %d\n",
                tf->evt_nroads[evt]);
      }
      for (id = 0; id <= XFT_LYR; id++)
        tf->evt_nhits[evt][tf->evt_nroads[evt]][id] = 0;
      tf->evt_err[evt][tf->evt_nroads[evt]] = 0;
      tf->evt_zid[evt][tf->evt_nroads[evt]] = -1;
      id = -1;
      id_last = -1;
    }
    else if (id == EE_LYR)
    {

#ifdef DEBUG_SVT
      printf ("id==EE_LYR --> ee_word = %.6x\n", ee_word);

#endif /*   */
      tf->evt_ee_word[evt] = ee_word;
      tf->totEvts++;

#ifdef DEBUG_SVT
      printf ("tf->evt_ee_word[%d] = %.6x\n", evt, tf->evt_ee_word[evt]);

#endif /*   */
      evt++;
      id = -1;
      id_last = -1;
    }
    else
    {
      tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << INV_DATA_BIT);
    }
    id_last = id;
  }                             //end loop on input words 

#ifdef DEBUG_SVT
  printf ("in gf_fep_unpack: found %d events in this file\n", tf->totEvts);

#endif /*   */

  return 0;
}


  /*  
   * Unpack the input words into the FEP  
   */
int
gf_fep_unpack_CPU (tf_arrays_t tf, int n_words_in, unsigned int *data_in)
{

  //   printf("in gf_fep_unpack_CPU\n"); 
  int id, gf_xft, id_last = -1;
  int ee_word, parity_in, sector, amroad, crv, crv_sign, phi, zid, lcl, hit,
    err_sum;
  int evt, rd, ie;
  evt = EVT;
  gf_xft = 0;
  tf->totEvts = 0;

  /* initializing arrays */
  for (ie = 0; ie < NEVTS; ie++)
  {
    tf->evt_nroads[ie] = 0;
    tf->evt_err_sum[ie] = 0;
    for (id = 0; id <= NSVX_PLANE; id++)
    {
      tf->evt_layerZ[ie][tf->evt_nroads[ie]][id] = 0;
    }
    for (id = 0; id <= XFT_LYR; id++)
    {
      tf->evt_nhits[ie][tf->evt_nroads[ie]][id] = 0;
    }
    tf->evt_err[ie][tf->evt_nroads[ie]] = 0;
    tf->evt_zid[ie][tf->evt_nroads[ie]] = -1;

    //    printf("tf->evt_nroads[%d] = %d, tf->evt_zid[%d][tf->evt_nroads[%d]] = %d\n", ie, tf->evt_nroads[ie], ie, ie, tf->evt_zid[ie][tf->evt_nroads[ie]]); 
  }

  //GF Now unpack the data 

  //#if DEBUG_SVT 
  printf ("in gf_fep_unpack_CPU: unpack data\n");

  //#endif 
  int i;
  for (i = 0; i < n_words_in; i++)
  {
    printf ("in gf_fep_unpack_CPU: data_in[%d] = %.6x\n", i, data_in[i]);
    id =
      word_decode (data_in[i], gf_xft, &ee_word, &parity_in, &sector,
                   &amroad, &crv, &crv_sign, &phi, &zid, &lcl, &hit,
                   &err_sum);

    /* SVX Data */
    if (id < XFT_LYR)
    {
      tf->evt_hit[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                [tf->evt_nroads[evt]][id]] =
        hit;
      tf->evt_hitZ[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                 [tf->evt_nroads[evt]][id]] =
        zid;
      tf->evt_lcl[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                [tf->evt_nroads[evt]][id]] =
        lcl;
      tf->evt_lclforcut[evt][tf->evt_nroads[evt]][id][tf->evt_nhits[evt]
                                                      [tf->evt_nroads[evt]]
                                                      [id]] = lcl;
      tf->evt_layerZ[evt][tf->evt_nroads[evt]][id] = zid;
      if (tf->evt_zid[evt][tf->evt_nroads[evt]] == -1)
      {
        tf->evt_zid[evt][tf->evt_nroads[evt]] = zid & gf_mask (GF_SUBZ_WIDTH);
      }
      else
      {
        tf->evt_zid[evt][tf->evt_nroads[evt]] =
          (((zid & gf_mask (GF_SUBZ_WIDTH)) << GF_SUBZ_WIDTH)
           +
           (tf->evt_zid[evt][tf->evt_nroads[evt]] & gf_mask (GF_SUBZ_WIDTH)));
      }
      tf->evt_nhits[evt][tf->evt_nroads[evt]][id]++;

      /* Error Checking */
      if (tf->evt_nhits[evt][tf->evt_nroads[evt]][id] == MAX_HIT)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OFLOW_HIT_BIT);
      if (id < id_last)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OUTORDER_BIT);
    }
    else if ((id == XFT_LYR) && (gf_xft == 0))
    {
      gf_xft = 1;
    }
    else if ((id == XFT_LYR) && (gf_xft == 1))
    {
      gf_xft = 0;
      tf->evt_crv[evt][tf->evt_nroads[evt]][tf->evt_nhits[evt][tf->evt_nroads
                                                               [evt]][id]] =
        crv;
      tf->evt_crv_sign[evt][tf->evt_nroads[evt]][tf->evt_nhits[evt]
                                                 [tf->evt_nroads[evt]][id]] =
        crv_sign;
      tf->evt_phi[evt][tf->evt_nroads[evt]][tf->evt_nhits[evt][tf->evt_nroads
                                                               [evt]][id]] =
        phi;
      tf->evt_nhits[evt][tf->evt_nroads[evt]][id]++;

      /* Error Checking */
      if (tf->evt_nhits[evt][tf->evt_nroads[evt]][id] == MAX_HIT)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OFLOW_HIT_BIT);
      if (id < id_last)
        tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << OUTORDER_BIT);
    }
    else if (id == EP_LYR)
    {
      tf->evt_cable_sect[evt][tf->evt_nroads[evt]] = sector;
      tf->evt_sect[evt][tf->evt_nroads[evt]] = sector;
      tf->evt_road[evt][tf->evt_nroads[evt]] = amroad;
      tf->evt_err_sum[evt] |= tf->evt_err[evt][tf->evt_nroads[evt]];
      tf->evt_nroads[evt]++;
      if (tf->evt_nroads[evt] > MAXROAD)
      {
        printf
          ("The limit on the number of roads fitted by the TF is %d\n",
           MAXROAD);
        printf ("You reached that limit evt->nroads = %d\n",
                tf->evt_nroads[evt]);
      }
      for (id = 0; id <= XFT_LYR; id++)
        tf->evt_nhits[evt][tf->evt_nroads[evt]][id] = 0;
      tf->evt_err[evt][tf->evt_nroads[evt]] = 0;
      tf->evt_zid[evt][tf->evt_nroads[evt]] = -1;
      id = -1;
      id_last = -1;
    }
    else if (id == EE_LYR)
    {
      tf->totEvts++;
      tf->evt_ee_word[evt] = ee_word;
      evt++;
      id = -1;
      id_last = -1;
    }
    else
    {
      tf->evt_err[evt][tf->evt_nroads[evt]] |= (1 << INV_DATA_BIT);
    }
    id_last = id;
  }                             //end loop on input words 

#ifdef DEBUG_SVT
  printf ("in gf_fep_unpack_CPU: found %d events in this file\n",
          tf->totEvts);

#endif /*   */

  return 0;
}

int
gf_fep_comb (tf_arrays_t tf)
{

  /* 
     combiner 
     hitmap (6 bits) = each bit corresponds to a single layer (c and phi 
     from XFT are considered together) 
     lcl (lclmap, 5 bits) = long cluster map; each bit corresponds to a  
     single SVXII layer. 
     zid (3 bits) = z coordinate of the innermost hit  
   */
  int ir, id, ic;               /* Several loop variables 
                                   ir : roads 
                                   id : layer 
                                   ic : combinations */
  int nlyr;                     /* The number of layers with a hit */
  int nmlyr;                    /* The number of layers with multiple hits */
  int icomb, ncomb;             /* The number of combinations */
  int ic2 = 0;                  /* ic2 == counter of combs that we will fit */
  int icomb5h = 0;              /* number of 4/5 tracks out of a 5/5 track */
  int comb2skip;
  int minlyr = MINHITS;
  int ie = 0;

  //  struct gf_ievt *evt = tf->evt; 
  //  struct fep_out *fep = tf->fep; 

  /* --------- Executable starts here ------------ */
#ifdef DEBUG_SVT
  printf ("minlyr: %d \n", minlyr);

#endif /*   */
  for (ie = 0; ie < tf->totEvts; ie++)
  {
    for (ir = 0; ir < tf->evt_nroads[ie]; ir++)
    {                           /* Loop for the number of roads */
      ncomb = 1;
      nlyr = 0;
      nmlyr = 0;
      ic2 = 0;

      /* At first, we calculate how many combinations are there */
      for (id = 0; id < (XFT_LYR + 1); id++)    /* Loop for the number of layers */
        if (tf->evt_nhits[ie][ir][id] != 0)
        {
          ncomb *= tf->evt_nhits[ie][ir][id];
          nlyr++;
          if (tf->evt_nhits[ie][ir][id] > 1)
            nmlyr++;
        }

      /* do we check this in the gf? */
      if (nlyr < minlyr)
      {
        tf->evt_err[ie][ir] |= (1 << UFLOW_HIT_BIT);
      }
      for (ic = 0; ic < ncomb; ic++)
      {                         /* Loop for all combinations */
        comb2skip = 0;
        icomb = ic;
        tf->fep_lcl[ie][ir][ic2] = 0;
        tf->fep_hitmap[ie][ir][ic2] = 0;
        tf->fep_phi[ie][ir][ic2] = 0;
        tf->fep_crv[ie][ir][ic2] = 0;
        tf->fep_lclforcut[ie][ir][ic2] = 0;
        for (id = 0; id < XFT_LYR; id++)
        {                       /* SVX Layers */
          tf->fep_hit[ie][ir][ic2][id] = 0;
          tf->fep_hitZ[ie][ir][ic2][id] = 0;
          if (tf->evt_nhits[ie][ir][id] != 0)
          {

#ifdef DEBUG_SVT
            printf
              ("in gf_fep_comb: tf->evt_hitZ[%d][%d][%d][%d] = %d\n",
               ie, ir, ic, icomb % tf->evt_nhits[ie][ir][id],
               tf->evt_hitZ[ie][ir][id][icomb % tf->evt_nhits[ie][ir][id]]);

#endif /*   */
            tf->fep_hit[ie][ir][ic2][id] =
              tf->evt_hit[ie][ir][id][icomb % tf->evt_nhits[ie][ir][id]];
            tf->fep_hitZ[ie][ir][ic2][id] =
              tf->evt_hitZ[ie][ir][id][icomb % tf->evt_nhits[ie][ir][id]];
            tf->fep_lcl[ie][ir][ic2] |=
              ((tf->evt_lcl[ie][ir][id][icomb %
                                        tf->evt_nhits[ie][ir][id]] &
                gf_mask (1)) << id);
            tf->fep_lclforcut[ie][ir][ic2] |=
              ((tf->evt_lclforcut[ie][ir][id][icomb %
                                              tf->evt_nhits[ie][ir]
                                              [id]] & gf_mask (1)) << id);
            icomb /= tf->evt_nhits[ie][ir][id];
            tf->fep_hitmap[ie][ir][ic2] |= (1 << id);
          }                     /* if (tf->h_NHITS[ie][ir][id] |= 0) { */

#ifdef DEBUG_SVT
          printf ("in gf_fep_comb: tf->fep_hitZ[%d][%d][%d] = %d\n",
                  ir, ic2, id, tf->fep_hitZ[ie][ir][ic2][id]);

#endif /*   */
        }                       /* for (id=0; id<XFT_LYR; id++) { */

        /* check if this is a 5/5 track */
        if (tf->fep_hitmap[ie][ir][ic2] != 0x1f)
        {
          tf->fep_ncomb5h[ie][ir][ic2] = 1;
        }

        else
        {
          tf->fep_ncomb5h[ie][ir][ic2] = 5;
        }

#ifdef DEBUG_SVT
        printf
          ("tf->fep_hitmap[%d][%d] = %.6X, tf->fep_ncomb5h = %d\n", ir,
           ic2, tf->fep_hitmap[ie][ir][ic2], tf->fep_ncomb5h[ie][ir][ic2]);

#endif /*   */
        id = XFT_LYR;           /* XFT layers */
        if (tf->evt_nhits[ie][ir][id] != 0)
        {
          tf->fep_phi[ie][ir][ic2] =
            (tf->evt_phi[ie][ir][icomb %
                                 tf->evt_nhits[ie][ir][id]] &
             gf_mask (SVT_PHI_WIDTH));
          tf->fep_crv[ie][ir][ic2] =
            (tf->evt_crv[ie][ir][icomb %
                                 tf->evt_nhits[ie][ir][id]] &
             gf_mask (SVT_CRV_WIDTH));
          tf->fep_crv_sign[ie][ir][ic2] =
            (tf->evt_crv_sign[ie][ir][icomb % tf->evt_nhits[ie][ir][id]]);
          icomb /= tf->evt_nhits[ie][ir][id];
        }
        tf->fep_zid[ie][ir] = (tf->evt_zid[ie][ir] & gf_mask (GF_ZID_WIDTH));
        tf->fep_road[ie][ir] =
          (tf->evt_road[ie][ir] & gf_mask (SVT_ROAD_WIDTH));
        tf->fep_sect[ie][ir] =
          (tf->evt_sect[ie][ir] & gf_mask (SVT_SECT_WIDTH));
        tf->fep_cable_sect[ie][ir] =
          (tf->evt_cable_sect[ie][ir] & gf_mask (SVT_SECT_WIDTH));
        for (icomb5h = 0; icomb5h < tf->fep_ncomb5h[ie][ir][ic2]; icomb5h++)
        {
          tf->fep_err[ie][ir][ic2][icomb5h] = tf->evt_err[ie][ir];
        }
        ic2++;

        /* we will cut on the LC in the g_func 
           if (comb2skip == 0) 
           ic2++; 
         */
      }                         /* for (ic=0; ic<ncomb; ic++) { */
      tf->fep_ncmb[ie][ir] = ic2;
      tf->evt_err_sum[ie] |= tf->evt_err[ie][ir];

#ifdef DEBUG_SVT
      printf ("in gf_fep_comb: tf->evt_err[ie][%d] = %.6x\n", ir,
              tf->evt_err[ie][ir]);

#endif /*   */
    }                           /* for (ir=0; ir<evt->nroads; ir++) { */
    tf->fep_nroads[ie] = tf->evt_nroads[ie];
    tf->fep_ee_word[ie] = tf->evt_ee_word[ie];
    tf->fep_err_sum[ie] = tf->evt_err_sum[ie];

#ifdef DEBUG_SVT
    printf ("in gf_fep_comb: tf->fep_err_sum = %.6x\n", tf->fep_err_sum[ie]);

#endif /*   */
  }                             //end loop on events 
  return SVTSIM_GF_OK;
}

int
gf_fit (tf_arrays_t tf)
{

  /* this function should  
     1) fill the following arrays: 
     - hits[SVTNHITS] from the combiner 
     - coeff[NFITTER][SVTNHITS] from the combiner 
     2)call gf_fit_proc which does the scalar product  
     it can handle both 4/5 and 5/5 tracks. 5/5 tracks are splitted 
     into 5 4/5 tracks.  
   */
  int ir, ic, ip, ih, ihit, il, i;
  int rc = 1;
  int hit[SVTNHITS];
  long long int coeff[NFITTER][SVTNHITS];
  int coe_addr, int_addr;       /* Address for coefficients and intercept */
  int mka_addr;                 /* Address for MKADDR memory */
  long long int theintcp = 0;
  int sign_crv = 0;
  int which = 0, lwhich = 0;
  int iz = 0;
  int newhitmap = 0;
  int g = 0;
  int p0[6], ix[7];
  int ie;

  //  struct fep_out *fep; 
  //struct fit_out *trk; 
  int map[7][7] = { {0, 1, 2, 3, -1, 4, 5},     /* 01235 */
  {0, 1, 2, -1, 3, 4, 5},       /* 01245 */
  {0, 1, -1, 2, 3, 4, 5},       /* 01345 */
  {0, -1, 1, 2, 3, 4, 5},       /* 02345 */
  {-1, 0, 1, 2, 3, 4, 5},       /* 12345 */
  {0, 1, 2, 3, -1, 4, 5},       /* (??) */
  {0, 1, 2, 3, -1, 4, 5}        /* (??) */
  };

  /* --------- Executable starts here ------------ */
  //  fep = tf->fep; 
  //trk = tf->trk; 
  for (ie = 0; ie < tf->totEvts; ie++)
  {

#ifdef DEBUG_SVT
    printf ("\n ******** Starting fit procedure ******** \n");

#endif /*   */
    tf->fit_err_sum[ie] = tf->fep_err_sum[ie];
    for (ir = 0; ir < tf->fep_nroads[ie]; ir++)
    {                           /* Loop for all roads */
      for (ic = 0; ic < tf->fep_ncmb[ie][ir]; ic++)
      {                         /* Loop for all combinations */

#ifdef DEBUG_SVT
        printf
          ("\n-------------> Processing road %d (ID = %.6x), combination %d\n",
           ir, tf->evt_road[ie][ir], ic);
        printf ("tf->fep_hitmap[%d][%d]: %.6x\n", ir, ic,
                tf->fep_hitmap[ie][ir][ic]);

#endif /*   */
        if (tf->fep_hitmap[ie][ir][ic] != 0x1f)
        {

#ifdef DEBUG_SVT
          printf ("4/5 track\n");

#endif /*   */
          for (ip = 0; ip < NFITTER; ip++)
          {
            for (ih = 0; ih < SVTNHITS; ih++)
            {
              if (gf_mkaddr
                  (tf, tf->fep_hitmap[ie][ir][ic],
                   tf->fep_lcl[ie][ir][ic],
                   tf->fep_zid[ie][ir], &coe_addr, &int_addr,
                   &mka_addr, tf->fit_err_sum))
              {

                //            svtsim_gf_err(gf, "gf_mkaddr has an error in tf_fit"); 
                return SVTSIM_GF_ERR;
              }

#ifdef DEBUG_SVT
              printf
                ("coe_addr = %d, int_addr = %d, mka_addr = %.6x \n",
                 coe_addr, int_addr, mka_addr);
              printf ("hit[%d] = %.6x\n", ih, tf->fep_hit[ie][ir][ic][ih]);

#endif /*   */
              int_addr = (int_addr << OFF_SUBA_LSB) + tf->fep_road[ie][ir];
              iz = tf->fep_zid[ie][ir] & 7;

#ifdef DEBUG_SVT
              printf ("in gf_fit: iz = %d\n", iz);

#endif /*   */
              which = coe_addr / 6;
              lwhich = which;
              svtsim_assert (which >= 0
                             && which < SVTSIM_NEL (tf->whichFit[iz]));
              which = tf->whichFit[iz][which];

#ifdef DEBUG_SVT
              printf ("in gf_fit: lwhich = %d, which = %d\n", lwhich, which);

#endif /*   */
              svtsim_assert (which >= 0
                             && which <
                             SVTSIM_NEL (tf->lfitpar[ip]
                                         [map[lwhich][ih]][iz]));

#ifdef DEBUG_SVT
              printf
                ("in gf_fit: ip = %d, map[%d][%d] = %d, iz = %d, which = %d\n",
                 ip, ip, ih, map[lwhich][ih], iz, which);

#endif /*   */
              coeff[ip][ih] = map[lwhich][ih] < 0 ? 0 :
                (tf->lfitparfcon[ip][map[lwhich][ih]][iz][which]);

#ifdef DEBUG_SVT
              printf ("in gf_fit: coeff[%d][%d] = %.6x\n", ip,
                      ih, coeff[ip][ih]);

#endif /*   */
              if (ih < NSVX_PLANE)
              {
                hit[ih] =
                  ((tf->fep_hit[ie][ir][ic][ih] << 1) + 1) & gf_mask (15);

#ifdef DEBUG_SVT
                printf ("hit[%d]: %.6x\n", ih, hit[ih]);

#endif /*   */
              }
              else if (ih == HIT_PHI)
              {
                hit[ih] = tf->fep_phi[ie][ir][ic];

#ifdef DEBUG_SVT
                printf
                  ("in gf_fit: tf->wedge[%d] =  %d, SVTSIM_XFTPHIBINS = %d, SVTSIM_NWEDGE = %d\n",
                   ie, tf->wedge[ie], SVTSIM_XFTPHIBINS, SVTSIM_NWEDGE);

#endif /*   */
                hit[ih] -= tf->wedge[ie] * SVTSIM_XFTPHIBINS / SVTSIM_NWEDGE;
                hit[ih] = ((hit[ih] << 3) + (1 << 2)) & gf_mask (15);

#ifdef DEBUG_SVT
                printf ("in gf_fit: phi =  %.6x\n", hit[ih]);

#endif /*   */
              }
              else if (ih == HIT_CRV)
              {
                sign_crv = tf->fep_crv_sign[ie][ir][ic];
                hit[ih] =
                  ((tf->fep_crv[ie][ir][ic] << 8) + (1 << 7)) & gf_mask (15);

#ifdef DEBUG_SVT
                printf
                  ("in gf_fit: curv =  %.6x, sign: %d \n", hit[ih], sign_crv);

#endif /*   */
              }
            }                   /* end for(ih = 0; ih < SVTNHITS; ih++){ */

            /* INTERCEPT */
            theintcp = tf->lfitparfcon[ip][6][iz][which] << 18;

#ifdef DEBUG_SVT
            printf ("in gf_fit: theintcp = %.6x, iz = %d\n", theintcp, iz);

#endif /*   */

#ifdef DEBUG_SVT
            printf ("\n FIT # %d\n", ip);
            for (g = 0; g < SVTNHITS; g++)
            {
              printf ("hit[%d] = %.6x\n", g, hit[g]);
              printf ("coeff[%d][%d] = %.6x\n", ip, g, coeff[ip][g]);
            }

#endif /*   */
            /* printf("in gf_fit: which = %d\n", which); */
            if (gf_fit_proc
                (hit, sign_crv, coeff[ip], theintcp,
                 &tf->fit_fit[ie][ip][ir][ic][0],
                 &tf->fit_err[ie][ir][ic][0]))
            {
              printf ("gf_fit_proc has an error in gf_fit\n");
              return SVTSIM_GF_ERR;
            }

#ifdef DEBUG_SVT
            printf ("in gf_fit: trk->err[%d][%d] = %.6x\n", ir,
                    ic, tf->fit_err[ie][ir][ic]);
            printf
              ("result of gf_fit_proc, before fit_format, trk->fit[%d][%d][%d][0]: %.6x \n",
               ip, ir, ic, tf->fit_fit[ie][ip][ir][ic][0]);

#endif /*   */
          }                     /* end for (ip=0; ip<NFITTER; ip++) { */
        }                       /*  end  if(tf->fep_hitmap[ie][ir][ic] != 0x1f) { */

        else
        {                       /* 5/5 track transformed in 5 4/5 tracks */

#ifdef DEBUG_SVT
          printf ("5/5 track\n");

#endif /*   */
          for (ip = 0; ip < NFITTER; ip++)
          {
            for (ih = 0; ih < NSVX_PLANE; ih++)
            {

#ifdef DEBUG_SVT
              printf ("Track # %d of a 5/5 track\n", ih);

#endif /*   */
              for (il = 0; il < NSVX_PLANE; il++)
              {                 /* one call to gf_fit_proc  
                                   for each ih value  
                                 */

                /* let's calculate the new hitmap */
                if (il != ih)
                {
                  switch (ih)
                  {
                  case 0:      /*  11110 */
                    newhitmap = 0x1e;
                    break;
                  case 1:      /*  11101 */
                    newhitmap = 0x1d;
                    break;
                  case 2:      /*  11011 */
                    newhitmap = 0x1b;
                    break;
                  case 3:      /*  10111 */
                    newhitmap = 0x17;
                    break;
                  case 4:      /*  01111 */
                    newhitmap = 0x0f;
                    break;
                  }
                  if (gf_mkaddr
                      (tf, newhitmap,
                       tf->fep_lcl[ie][ir][ic],
                       tf->fep_zid[ie][ir], &coe_addr,
                       &int_addr, &mka_addr, tf->fit_err_sum))
                  {
                    printf ("gf_mkaddr has an error in tf_fit\n");
                    return SVTSIM_GF_ERR;
                  }

#ifdef DEBUG_SVT
                  printf
                    ("coe_addr = %d, int_addr = %d, mka_addr = %.6x \n",
                     coe_addr, int_addr, mka_addr);
                  printf
                    ("in gf_fit: tf->fep_hitZ[%d][%d][%d] = %d\n",
                     ir, ic, ih, tf->fep_hitZ[ie][ir][ic][ih]);

#endif /*   */
                  if (ih == 0)
                  {
                    iz = tf->fep_hitZ[ie][ir][ic][1];
                  }

                  else
                  {
                    iz = tf->fep_zid[ie][ir] & 7;
                  }
                  which = coe_addr / 6;

#ifdef DEBUG_SVT
                  printf ("in gf_fit: iz = %d\n", iz);

#endif /*   */
                  lwhich = which;
                  svtsim_assert (which >= 0
                                 && which < SVTSIM_NEL (tf->whichFit[iz]));
                  which = tf->whichFit[iz][which];
                  svtsim_assert (which >= 0
                                 && which <
                                 SVTSIM_NEL (tf->lfitpar[ip]
                                             [map[lwhich][ih]][iz]));
                  coeff[ip][il] =
                    map[lwhich][il] <
                    0 ? 0 : (tf->lfitparfcon[ip][map[lwhich][il]][iz][which]);

#ifdef DEBUG_SVT
                  printf
                    ("in gf_fit: coeff[%d][%d] = %.6x, %d\n",
                     ip, il, coeff[ip][il], coeff[ip][il]);

#endif /*   */
                  hit[il] =
                    ((tf->fep_hit[ie][ir][ic][il] << 1) + 1) & gf_mask (15);
                }

                else
                {
                  hit[il] = 0;
                  coeff[ip][il] = 1;
                }

#ifdef DEBUG_SVT
                printf ("hit[%d]: %.6x\n", il, hit[il]);

#endif /*   */
              }                 /* end for(il = 0; il <  NSVX_PLANE; il++) { */

              /* now we have to retrieve the curv and phi and the  
                 corresponding coefficients and the intercept 
               */
              hit[HIT_PHI] = tf->fep_phi[ie][ir][ic];
              hit[HIT_PHI] -=
                tf->wedge[ie] * SVTSIM_XFTPHIBINS / SVTSIM_NWEDGE;
              hit[HIT_PHI] = ((hit[HIT_PHI] << 3) + (1 << 2)) & gf_mask (15);
              coeff[ip][HIT_PHI] =
                map[lwhich][HIT_PHI] <
                0 ? 0 : (tf->lfitparfcon[ip][map[lwhich][HIT_PHI]]
                         [iz][which]);
              sign_crv = tf->fep_crv_sign[ie][ir][ic];
              hit[HIT_CRV] =
                ((tf->fep_crv[ie][ir][ic] << 8) + (1 << 7)) & gf_mask (15);
              coeff[ip][HIT_CRV] =
                map[lwhich][HIT_CRV] <
                0 ? 0 : (tf->lfitparfcon[ip][map[lwhich][HIT_CRV]]
                         [iz][which]);

              /* INTERCEPT */
              theintcp = tf->lfitparfcon[ip][6][iz][which] << 18;

#ifdef DEBUG_SVT
              printf ("in gf_fit: theintcp = %.6x\n", theintcp);

#endif /*   */
#ifdef DEBUG_SVT
              printf ("\n FIT # %d of track %d of a 5/5\n", ip, ih);

#endif /*   */
              if (gf_fit_proc
                  (hit, sign_crv, coeff[ip], theintcp,
                   &tf->fit_fit[ie][ip][ir][ic][ih],
                   &tf->fit_err[ie][ir][ic][ih]))
              {
                printf ("gf_fit_proc has an error in gf_fit\n");
                return SVTSIM_GF_ERR;
              }

#ifdef DEBUG_SVT
              printf ("in gf_fit: trk->err[%d][%d] = %.6x\n",
                      ir, ic, tf->fit_err[ie][ir][ic][ih]);
              printf
                ("after gf_fit_proc: comb # %d,tf->fit_fit[%d][%d][%d][%d] =  %d\n",
                 il, ip, ir, ic, ih, tf->fit_fit[ie][ip][ir][ic][ih]);

#endif /*   */
              tf->fit_err_sum[ie] |= tf->fit_err[ie][ir][ic][ih];
            }                   /*end for(ih = 0; ih < NSVX_PLANE; ih++){ */
          }                     /* end for (ip=0; ip<NFITTER; ip++) { */
        }                       /* end if(tf->fep_hitmap[ie][ir][ic] != 0x1f) { */
      }                         /* end for (ic=0; ic<tf->fep_ncmb[ie][ir]; ic++) { Loop for all combi */
    }                           /* end for (ir=0; ir<tf->fep_nroads; ir++) { Loop for all roads */
  }

#ifdef DUMP
  dump_fit (tf);

#endif /*   */
  gf_fit_format (tf);
  return rc;
}

int
gf_mkaddr (tf_arrays_t tf, int hitmap, int lclmap, int zmap, int *coe_addr,
           int *int_addr, int *addr, int *err)
{
  int iaddr;
  unsigned int datum = 0;

  //uint2  

  /* --------- Executable starts here ------------ */
  if ((hitmap < 0) || (hitmap > gf_mask (NSVX_PLANE + 1)) ||    /* + XFT_LYR */
      (lclmap < 0) || (lclmap > gf_mask (NSVX_PLANE)) ||
      (zmap < 0) || (zmap > gf_mask (GF_ZID_WIDTH)))
    *err |= (1 << SVTSIM_GF_MKADDR_INVALID);
  iaddr = ((zmap & gf_mask (GF_SUBZ_WIDTH))
           + (lclmap << MADDR_NCLS_LSB) + (hitmap << MADDR_HITM_LSB));
  if ((iaddr < 0) || (iaddr >= MAXMKA))
  {
    printf ("gf_mkaddr; invalid MKADDR address\n");
    return SVTSIM_GF_ERR;
  }
  int ldat = 0;
  svtsim_get_gfMkAddr (tf, &ldat, 1, iaddr);
  datum = ldat;
  *int_addr = datum & gf_mask (OFF_SUBA_WIDTH);
  *coe_addr = (datum >> OFF_SUBA_WIDTH) & gf_mask (PAR_ADDR_WIDTH);
  *addr = iaddr;
  return SVTSIM_GF_OK;
}



  /* 
   * fit (no shifts) 
   */
int
gf_fit_proc (int hit[], int sign_crv, long long int coeff[],
             long long int intcp, long long int *result, int *err)
{

  /* 
     the fitter calculates 6 scalar products in parallel: 3 chisq, 3 parameters 
     in uscita ho una parola del tipo: 
     d=10bit+segno, c=7+segno, phi=13, chi1,2,3=14 
   */

  /* 
     Simulate the fit processor in the board. 
     INPUT: hit[]; the hits from SVX and XFT 
     coeff[]; Coefficienct parameters for the linear fit. 
     intcp;   Intercept parameter for the linear fit. 
     OUTPUT: *result; The result of the fit. 
     *err; the error condition. 
   */
  int samesign = 0;
  long long int temp = 0;
  int i = 0;
  int sign_coeff;
  long long int partRes = 0;
  int g = 0;

  /* --------- Executable starts here ------------ */
  *result = 0;
  *err = 0;

#ifdef DEBUG_SVT
  printf ("\n GF_FIT_PROC\n");
  for (g = 0; g < SVTNHITS; g++)
  {
    printf ("hit[%d] = %.6x\n", g, hit[g]);
    printf ("coeff[%d] = %.6x\n", g, coeff[g]);
  }

#endif /*   */
  for (i = 0; i < SVTNHITS; i++)
  {
    if (i < NSVX_PLANE)
    {
      temp += hit[i] * coeff[i];
      partRes = hit[i] * coeff[i];

#ifdef DEBUG_SVT
      printf
        ("coeff[%d] = %.6llx, hit[%d] = %.6x, partRes=%.6llx, temp=%.6llx\n",
         i, coeff[i], i, hit[i], partRes, temp);

#endif /*   */
    }                           /*  end if (i < NSVX_PLANE) { */

    else if (i == HIT_PHI)
    {                           /* XFT phi */
      hit[i] = (hit[i] & 0x400) ? -((~hit[i] & 0x3ff) + 1) : (hit[i] & 0x3ff);
      temp += hit[i] * coeff[i];
      partRes = hit[i] * coeff[i];

#ifdef DEBUG_SVT
      printf
        ("coeff[%d] = %.6llx,  hit[%d] = %.6x, partRes=%.6llx, temp=%.6llx\n",
         i, coeff[i], i, hit[i], partRes, temp);

#endif /*   */
    }

    else if (i == HIT_CRV)
    {                           /* XFT curvature (curv already with sign in  
                                   fep */
      if (sign_crv == 1)
      {                         /* if negative bit is set */
        temp -= hit[i] * coeff[i];
        partRes = (hit[i] * coeff[i]);
      }

      else
      {
        temp += hit[i] * coeff[i];
        partRes = (hit[i] * coeff[i]);
      }

#ifdef DEBUG_SVT
      printf
        ("coeff[%d] = %.6llx, hit[%d] = %.6x, partRes=%.6llx, temp=%.6llx\n",
         i, coeff[i], i, hit[i], partRes, temp);

#endif /*   */
    }
  }

#ifdef DEBUG_SVT
  printf ("Intercept: %.6x \n", intcp);

#endif /*   */
  *result += temp;

#ifdef DEBUG_SVT
  printf ("result (w/o intcp) = %.6llx\n", *result);

#endif /*   */
  samesign = (*result < 0) == (intcp < 0);
  *result += intcp;

#ifdef DEBUG_SVT
  printf ("result (w intcp) = %.12llx\n", *result);

#endif /*   */
  *result = *result < 0 ? -((-*result) >> 17) : *result >> 17;

#ifdef DEBUG_SVT
  printf ("result (after shift >> 17) = %.12llx\n", *result);

#endif /*   */

  /* 
     if (samesign && abs(*result)>gf_mask3(FIT_DWIDTH)) { 
     if (1) { 
     printf(stderr, "gf_fit_proc: result overflow!\n"); 
     printf(stderr, "  intcp=%d result=%d DWIDTH=%d\n", 
     intcp, *result, FIT_DWIDTH); 
     } 
     *err |= (1<<FIT_RESULT_OFLOW_BIT); 
     printf("in gf_fit_proc: fit result overflow!\n"); 
     } 
   */
  if (*result > 0)
    *result &= gf_mask3 (FIT_DWIDTH);

  else
    *result = -(abs (*result) & gf_mask3 (FIT_DWIDTH));

#ifdef DEBUG_SVT
  printf ("in gf_fit_proc: result = %.6x\n", *result);

#endif /*   */
  return SVTSIM_GF_OK;
}

int
gf_fit_format (tf_arrays_t tf)
{
  /* 
     Format the result of the fit 
     - az. angle phi = 13 bits,  no  sign,  phi[0:12] = value, shift = 6 bits  
     - curvature   c =  8 bits,  c[7] =  sign, c[0:7] = value, shift = 6 bits  
     - imp. param. d = 11 bits, d[10] =  sign, d[0:9] = value, shift = 6 bits   
     - chi1,2,3  chi = 14 bits,  no  sign,  chi[0:13] = value, shift = 6 bits  
   */
  /* SA for the moment I'm using the structure opp */
  int ip, ie;
  long long int temp = 0;
  long int r = 0;
  long long int rr = 1;
  int ir, ic, ich;
  int ichi;

#ifdef DEBUG_SVT
  printf ("\n***************** Formatting fit results ****************\n");

#endif /*   */
  for (ie = 0; ie < tf->totEvts; ie++)
  {
    for (ir = 0; ir < tf->fep_nroads[ie]; ir++)
    {                           /* Loop for all roads */
      for (ic = 0; ic < tf->fep_ncmb[ie][ir]; ic++)
      {                         /* Loop for all combinations */

#ifdef DEBUG_SVT
        printf
          ("\n-------------> Processing road %d (ID = %.6x), combination %d\n",
           ir, tf->evt_road[ie][ir], ic);

#endif /*   */
        for (ich = 0; ich < tf->fep_ncomb5h[ie][ir][ic]; ich++)
        {

#ifdef DEBUG_SVT
          printf ("Track # %d of a %d/5\n", ich, tf->fep_ncomb5h[ie][ir][ic]);

#endif /*   */
          /* phi */
          temp = tf->fit_fit[ie][0][ir][ic][ich];

#ifdef DEBUG_SVT
          printf ("phi : %.6x\n", temp);

#endif /*   */
          if (tf->fit_fit[ie][0][ir][ic][ich] > 0)
          {
            temp = (tf->fit_fit[ie][0][ir][ic][ich] + 1);
            temp = temp >> 1;
          }

          else
          {
            temp = (tf->fit_fit[ie][0][ir][ic][ich] - 1);
            temp = -((-temp) >> 1);
          }

#ifdef DEBUG_SVT
          printf ("phi after rounding: %.6x\n", temp);

#endif /*   */
          /* overflow check: 
             if the fit results requires > N to be represented,  
             we set the overflow bit (N is the number of bits to be  
             considered after in output). 
           */
          if (abs (temp) > gf_mask (OPHI_WIDTH))
          {

#ifdef DEBUG_SVT
            printf ("gf_fit_format:phi overflow!\n");

#endif /*   */
            tf->fit_err[ie][ir][ic][ich] |= (1 << FIT_RESULT_OFLOW_BIT);
          }
          temp = (temp < 0 ? -temp : temp) & gf_mask (OPHI_WIDTH);

#ifdef DEBUG_SVT
          printf ("phi fixed width: %.6x\n", temp);
          printf ("\n\n");

#endif /*   */
          tf->fit_fit[ie][0][ir][ic][ich] = temp;

          /* impact parameter */
          temp = tf->fit_fit[ie][1][ir][ic][ich];

#ifdef DEBUG_SVT
          printf ("imp : %.6x\n", temp);

#endif /*   */
          if (tf->fit_fit[ie][1][ir][ic][ich] > 0)
          {
            temp = (tf->fit_fit[ie][1][ir][ic][ich] + 1);
            temp = temp >> 1;
          }

          else
          {
            temp = (tf->fit_fit[ie][1][ir][ic][ich] - 1);
            temp = -((-temp) >> 1);
          }

#ifdef DEBUG_SVT
          printf ("imp after rounding: %.6x\n", temp);

#endif /*   */
          /*overflow check */
          if (abs (temp) > gf_mask (OIMP_WIDTH))
          {
            printf ("gf_fit_format:impact parameter overflow!\n");
            tf->fit_err[ie][ir][ic][ich] |= (1 << FIT_RESULT_OFLOW_BIT);
          }
          temp = (temp < 0 ? -temp : temp) & gf_mask (OIMP_WIDTH);

#ifdef DEBUG_SVT
          printf ("imp fixed width: %.6x\n", temp);

#endif /*   */

          /* now add a bit for the sign  */
          if (tf->fit_fit[ie][1][ir][ic][ich] < 0)
          {

#ifdef DEBUG_SVT
            printf ("negative number\n");

#endif /*   */
            temp += (1 << OIMP_SIGN);
          }

#ifdef DEBUG_SVT
          printf ("imp w sign: %.6x\n", temp);
          printf ("\n\n");

#endif /*   */
          tf->fit_fit[ie][1][ir][ic][ich] = temp;

          /* curvature */
          temp = tf->fit_fit[ie][2][ir][ic][ich];

#ifdef DEBUG_SVT
          printf ("crv : %.6x\n", temp);

#endif /*   */
          if (tf->fit_fit[ie][2][ir][ic][ich] > 0)
          {
            temp = (tf->fit_fit[ie][2][ir][ic][ich] + 1);
            temp = temp >> 1;
          }

          else
          {
            temp = (tf->fit_fit[ie][2][ir][ic][ich] - 1);
            temp = -((-temp) >> 1);
          }

#ifdef DEBUG_SVT
          printf ("crv after rounding: %.6x\n", temp);

#endif /*   */
          /*overflow check */
          if (abs (temp) > gf_mask (OCVR_WIDTH))
          {

/*      	 printf("gf_fit_format:curvature overflow!\n"); */
            tf->fit_err[ie][ir][ic][ich] |= (1 << FIT_RESULT_OFLOW_BIT);
          }
          temp = (temp < 0 ? -temp : temp) & gf_mask (OCVR_WIDTH);

#ifdef DEBUG_SVT
          printf ("crv fixed width: %.6x\n", temp);

#endif /*   */
          /*  now add a bit for the sign  */
          if (tf->fit_fit[ie][2][ir][ic][ich] < 0)
          {

#ifdef DEBUG_SVT
            printf ("negative number\n");

#endif /*   */
            temp += (1 << OCVR_SIGN);
          }

#ifdef DEBUG_SVT
          printf ("crv w sign: %.6x\n", temp);
          printf ("\n\n");

#endif /*   */
          tf->fit_fit[ie][2][ir][ic][ich] = temp;

          /* chi 1,2,3 */
          for (ichi = 3; ichi < 6; ichi++)
          {
            temp = tf->fit_fit[ie][ichi][ir][ic][ich];

#ifdef DEBUG_SVT
            printf ("chi %d : %.6x\n", ichi - 2, temp);

#endif /*   */
            /*      temp = temp & gf_mask(OCHI_WIDTH); */
#ifdef DEBUG_SVT
            printf ("chi %d fixed width: %.6x\n", ichi - 2, temp);
            printf ("\n\n");

#endif /*   */
            tf->fit_fit[ie][ichi][ir][ic][ich] = temp;
          }
        }                       /* end for(ich = 0; ich < fep->ncomb5h[ir][ic]; ich++){ */
      }                         /* end for (ic=0; ic<fep->ncmb[ir]; ic++) {  Loop for all combi */
    }                           /* end for (ic=0; ic<fep->ncmb[ir]; ic++) {  Loop for all roads */
  }                             //end loop on events 

  return 0;
}

int
gf_init (tf_arrays_t * ptr_tf)
{

  //GF 
  struct tf_arrays *tf = malloc (sizeof (struct tf_arrays));
  tf->out = svtsim_cable_new ();
  tf->gf_emsk = 0;
  tf->chi2cut = GF_TOTCHI2_CUTVAL;
  tf->svt_emsk = GF_ERRMASK_SVT;
  tf->cdf_emsk = GF_ERRMASK_CDF;
  tf->eoe_emsk = GF_ERRMASK_EOE;        /* MASK for the errors */
  *ptr_tf = tf;

  return 0;
} 
