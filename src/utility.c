#include "private_includes.h"

void
dump_fout (tf_arrays_t tf)
     /*  
      *   Dump FORMATTER output 
      *     INPUT: fout; FORMATTER output 
      *          */
{
  int it, i, ie;
  for (ie = 0; ie < tf->totEvts; ie++)
  {
    printf ("============== dump_fout ==================\n");
    printf ("#SVT tracks = %d \n", tf->fout_ntrks[ie]);
    for (it = 0; it < tf->fout_ntrks[ie]; it++)
    {
      printf (" ----------------------------------------------\n");
      printf (" road %4X, comb %2d :\n GF OUTPUTs: ",
              tf->fout_iroad[ie][it], tf->fout_icmb[ie][it]);
      printf ("\t");
      for (i = 0; i < 4; i++)
        printf ("%8lX, ", tf->fout_gfword[ie][it][i]);
      printf ("\n\t \t");
      for (i = 4; i < 7; i++)
        printf ("%8lX, ", tf->fout_gfword[ie][it][i]);
    }                           /* for (it = 0; it<opp->ntrks; it++) { */
    printf ("\n ----------------------------------------------\n");
    printf ("EE word : %X,  SVT Error : %d,  OPP Error summary : %X\n",
            tf->fout_ee_word[ie], tf->fout_svterr[ie], tf->fout_err_sum[ie]);
    printf ("\n");
  }                             //end loop on events 
}


/*  
 * Calculate the parity of the word 
 */
int
cal_parity (int word)
{
  int i, par = 0;
  for (i = 0; i < SVT_WORD_WIDTH; i++)
    par ^= ((word >> i) & gf_mask (1));
  return par;
}


void
dump_gf_ievt (tf_arrays_t tf)
 /*  
    Dump the GF input words after decoding. 
    INPUT: gf; Decoded input words 
  */
{
  int ie, ir, id, ih;

  /* --------- Executable starts here ------------ */
  fprintf (stdout, "============== dump_gf_ievt ====================\n");
  for (ie = 0; ie < tf->totEvts; ie++)
  {
    fprintf (stdout, "#road = %d,  EE word: 0x%X,  ERR: %d\n",
             tf->evt_nroads[ie], tf->evt_ee_word[ie], tf->evt_err_sum[ie]);
    for (ir = 0; ir < tf->evt_nroads[ie]; ir++)
    {
      fprintf (stdout, " ----------------------------------------------\n");
      fprintf (stdout,
               " road %4d: ID = %.6x, SECT = %d, ERR = %.6x, Z(in,out) = (%d,%d)\n",
               ir + 1, tf->evt_road[ie][ir], tf->evt_sect[ie][ir],
               tf->evt_err[ie][ir], tf->evt_zid[ie][ir] & gf_mask (3),
               (tf->evt_zid[ie][ir] >> 3) & gf_mask (3));
      for (id = 0; id < XFT_LYR; id++)
      {
        fprintf (stdout, "  Layer   %d: %d hits => ", id,
                 tf->evt_nhits[ie][ir][id]);
        for (ih = 0; ih < tf->evt_nhits[ie][ir][id]; ih++)
        {
          if (tf->evt_lcl[ie][ir][id][ih] == 1)
          {
            fprintf (stdout, "-%.6x, ", tf->evt_hit[ie][ir][id][ih]);
          }
          else
          {
            fprintf (stdout, "%.6x, ", tf->evt_hit[ie][ir][id][ih]);
          }
        }                       /* for (ih=0;ih<tf->evt_nhits[ie][ir][id];ih++) { */
        fprintf (stdout, "\n");
      }                         /* for (id=0; id< XFT_LYR; id++) { */
      fprintf (stdout, "  Layer XFT: %d hits => ",
               tf->evt_nhits[ie][ir][XFT_LYR]);
      for (ih = 0; ih < tf->evt_nhits[ie][ir][XFT_LYR]; ih++)
      {

        //printf("tf->evt_crv_sign[%d][%d][%d] = %d\n", ie,ir,ih,tf->evt_crv_sign[ie][ir][ih]); 
        //printf("tf->h_CRV[%d][%d][%d] = %d\n", ie,ir,ih,tf->h_CRV[ie][ir][ih]); 
        //printf("tf->evt_phi[%d][%d][%d] = %d\n", ie,ir,ih,tf->evt_phi[ie][ir][ih]); 
        if (tf->evt_crv_sign[ie][ir][ih] == 0)
        {
          fprintf (stdout, "CRV=+%.6x PHI=%.6x, ",
                   tf->evt_crv[ie][ir][ih], tf->evt_phi[ie][ir][ih]);
        }
        else if (tf->evt_crv_sign[ie][ir][ih] == 1)
        {
          fprintf (stdout, "CRV=-%.6x PHI=%.6x, ",
                   tf->evt_crv[ie][ir][ih], tf->evt_phi[ie][ir][ih]);
        }
      }
      fprintf (stdout, "\n");
    }
    fprintf (stdout, "\n");
  }                             //end for(ie = 0; ie < NEVTS; ie++) { 
}

void
dump_fep_out (tf_arrays_t tf)
 /*  
    Dump FEP output 
    INPUT: fep; FEP output 
  */
{
  int ir, id, ic;
  int ie;

  /* --------- Executable starts here ------------ */
  for (ie = 0; ie < tf->totEvts; ie++)
  {
    printf ("============== dump_fep_out ==================\n");
    printf ("#road = %d,  EE word: 0x%X,  ERR: %d\n", tf->fep_nroads[ie],
            tf->fep_ee_word[ie], tf->fep_err_sum[ie]);
    for (ir = 0; ir < tf->fep_nroads[ie]; ir++)
    {
      printf (" ----------------------------------------------\n");
      printf (" road %4d: #COMB=%d, ID = %.6x, SECT = %d, Z = %d\n",
              ir + 1, tf->fep_ncmb[ie][ir], tf->fep_road[ie][ir],
              tf->fep_sect[ie][ir], tf->fep_zid[ie][ir]);
      for (ic = 0; ic < tf->fep_ncmb[ie][ir]; ic++)
      {
        printf ("  +++++++++++++++++++++++++++++++\n");
        printf ("  Comb %2d: ", ic + 1);
        printf ("  hits => ");
        for (id = 0; id < XFT_LYR; id++)
          printf (" %.6x, ", tf->fep_hit[ie][ir][ic][id]);
        printf
          ("\n             CRV=%.6x, PHI=%.6x, LCLS=%.6x, HITMAP=%.6x, ERR=%d\n",
           tf->fep_crv[ie][ir][ic], tf->fep_phi[ie][ir][ic],
           tf->fep_lcl[ie][ir][ic], tf->fep_hitmap[ie][ir][ic],
           tf->fep_err[ie][ir][ic]);
      }
    }                           /* for (ir = 0; ir<tf->evt->nroads; ir++) { */
    printf ("\n");
  }
  return;
}

void
dump_fit (tf_arrays_t tf)
     /*  
        Dump FIT output 
        INPUT: fep; FEP output (We need the number of roads and combinations.) 
        fit; FIT output 
      */
{
  int ir, ic, i, ich, ie;
  for (ie = 0; ie < tf->totEvts; ie++)
  {
    printf ("============== dump_fit ==================\n");
    printf ("#road = %d, #comb=", tf->fep_nroads[ie]);
    for (ir = 0; ir < tf->fep_nroads[ie]; ir++)
      printf ("%d ", tf->fep_ncmb[ie][ir]);
    printf (": \n");
    printf (" ----------------------------------------------\n");
    for (ir = 0; ir < tf->fep_nroads[ie]; ir++)
    {
      for (ic = 0; ic < tf->fep_ncmb[ie][ir]; ic++)
      {
        for (ich = 0; ich < tf->fep_ncomb5h[ie][ir][ic]; ich++)
        {
          printf
            (" road %4d, roadID %.6x, comb %2d, which 4/5 %2d : FIT = ",
             ir + 1, tf->evt_road[ie][ir], ic + 1, ich);
          for (i = 0; i < NFITTER; i++)
          {
            if (tf->fit_fit[ie][i][ir][ic][ich] >= 0)
              printf ("%4X, ", tf->fit_fit[ie][i][ir][ic][ich]);

            else
              printf ("%4X, ",
                      abs (tf->fit_fit[ie][i][ir][ic][ich]) + 0x2000);
          }
          printf ("\n");
        }
      }
    }                           /* for (ir = 0; ir<tf->evt->nroads; ir++) { */
    printf ("\n");
  }                             //end loop on events 
}


  /*** END: ***/
int
word_decode (int word, int xft, int *ee_word, int *parity_in, int *sector,
             int *amroad, int *crv, int *crv_sign, int *phi, int *zid,
             int *lcl, int *hit, int *err)
{
  int ee, ep, lyr;

  /* --------- Executable starts here ------------ */
  lyr = -999;                   /* Any invalid numbers != 0-7 */
  if (word > gf_mask (SVT_WORD_WIDTH))
  {
    printf
      ("gf_iword_decode: Input data is larger than the maximum SVT word");
    return SVTSIM_GF_ERR;
  }

  /* check if this is a EP or EE word */
  ee = (word >> SVT_EE_BIT) & gf_mask (1);
  ep = (word >> SVT_EP_BIT) & gf_mask (1);

#ifdef DEBUG_SVT
  printf ("in word decode: word = %.6x\n", word);

#endif /*   */
  if (ee && ep)
  {                             /* End of Event word */

#ifdef DEBUG_SVT
    printf ("in word decode: it's a ee word\n");

#endif /*   */
    *ee_word = word;
    *parity_in = (word >> SVT_PAR_BIT) & gf_mask (1);
    lyr = EE_LYR;
  }
  else if (ee)
  {                             /* only EE bit ON is error condition */
    *err |= (1 << UNKNOWN_ERR);
    lyr = EE_LYR;               /* We have to check */
  }
  else if (ep)
  {                             /* End of Packet word */

#ifdef DEBUG_SVT
    printf ("in word decode: it's a ep word\n");

#endif /*   */
    lyr = EP_LYR;
    *sector = 6;

    /*   *sector = (word >> SVT_SECT_LSB)  & gf_mask(SVT_SECT_WIDTH); */
    *amroad = word & gf_mask (AMROAD_WORD_WIDTH);
  }
  else if (xft)
  {                             /* Second XFT word */

#ifdef DEBUG_SVT
    printf ("in word decode: it's second XFT word\n");

#endif /*   */
    *crv = (word >> SVT_CRV_LSB) & gf_mask (SVT_CRV_WIDTH);
    *crv_sign = (word >> (SVT_CRV_LSB + SVT_CRV_WIDTH)) & gf_mask (1);
    *phi = word & gf_mask (SVT_PHI_WIDTH);
    lyr = XFT_LYR;

#ifdef DEBUG_SVT
    printf ("in word decode: lyr = %d\n", lyr);

#endif /*   */
  }
  else
  {                             /* SVX hits or the first XFT word */

#ifdef DEBUG_SVT
    printf ("in word decode: it's  SVX hits or the first XFT word\n");

#endif /*   */
    lyr = (word >> SVT_LYR_LSB) & gf_mask (SVT_LYR_WIDTH);
    if (lyr == XFT_LYRID)
      lyr = XFT_LYR;
    *zid = (word >> SVT_Z_LSB) & gf_mask (SVT_Z_WIDTH);
    *lcl = (word >> SVT_LCLS_BIT) & gf_mask (1);
    *hit = word & gf_mask (SVT_HIT_WIDTH);
  }
  return lyr;
}
