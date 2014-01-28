#include "private_includes.h"

/* 
 * Compute GF make-address map 
 */
int
svtsim_get_gfMkAddr (tf_arrays_t tf, int *d, int nd, int d0)
{

  /*  
     d0 = iaddr 

   */
  int j = 0;
  int md = 0x4000;
  int iz, lcl, hit;
  int nparam = 6;
  int ipar;
  int dimspa = 6;
  int idim;
  int nbarrel = 6;
  int ibar;
  int nfitBlock = 32;
  int ifitBlock;
  ibar = 0;
  ifitBlock = 0;
  if (d0 + nd > md)
    nd = md - d0;
  for (j = 0; j < nd; j++)
  {
    int i = j + d0;
    int word = 0xffff, intcp = 0, coeff = 0;
    int which;
    iz = i & 7, lcl = i >> 3 & 0x1f, hit = i >> 8 & 0x3f;
    which = svtsim_whichFit (tf, iz, hit, lcl);
    coeff = iz + which * 6;     /* poor choice for illegal iz=6,7, but compatible */
    intcp = which;

#ifdef DEBUG_SVT
    printf ("in svtsim_get_gfMkAddr: nd = %d, d0 = %d\n", nd, d0);
    printf ("in svtsim_get_gfMkAddr: which = %d, coeff = %d, intcp = %d\n",
            which, coeff, intcp);

#endif /*   */
    word = coeff << 3 | intcp;
    d[j] = word;
  } return nd;
}


svtsim_cable_t *
svtsim_cable_new (void)
{
  svtsim_cable_t *cable = 0;
  cable = svtsim_malloc (sizeof (*cable));
  cable->data = 0;
  cable->ndata = 0;
  cable->mdata = 0;
  return cable;
}

void
svtsim_cable_addwords (svtsim_cable_t * cable, unsigned int *word, int nword)
{
  const int minwords = 8;
  int nnew = 0;
  svtsim_assert (cable);
  nnew = cable->ndata + nword;
  if (nnew > cable->mdata)
  {
    cable->mdata = SVTSIM_MAX (minwords, 2 * nnew);
    cable->data =
      svtsim_realloc (cable->data, cable->mdata * sizeof (cable->data[0]));
  }
  if (word)
  {
    memcpy (cable->data + cable->ndata, word, nword * sizeof (word[0]));
  }
  else
  {
    memset (cable->data + cable->ndata, 0, nword * sizeof (word[0]));
  }
  cable->ndata += nword;

  //  printf("in svtsim_cable_addwords: cable->ndata = %d\n", cable->ndata); 
}


void
svtsim_cable_addword (svtsim_cable_t * cable, unsigned int word)
{
  svtsim_cable_addwords (cable, &word, 1);
} 


void
svtsim_cable_copywords (svtsim_cable_t * cable, unsigned int *word, int nword)
{
  svtsim_assert (cable);
  cable->ndata = 0;
  svtsim_cable_addwords (cable, word, nword);
} 


/* 
 * Map layerMask, lcMask into fit block (0..31) for fully general case 
 * in which all fit constants are available.  Note that blocks 25..30 
 * are set up for some exotic cases that may never really occur.  Note 
 * also that before 2002-12-15 or so, only fit block 0 was actually 
 * used, and that as of 2002-12-15, no fit constants yet exist (or are 
 * planned) that make use of the LC bits. 

 * Note even further that the upgraded TF has NO room for fit constants 
 * that use the LC bits, but we do use them  
 * to select the constants when we have 5/5 roads 
 */
int
svtsim_whichFit_full (tf_arrays_t tf, int zin, int layerMask, int lcMask)
{

  /* 
   * Key (*=LC, 5=XFT): 
   *   0  =>  0  1  2  3  5 
   *   1  =>  0  1  2  3* 5 
   *   2  =>  0  1  2* 3  5 
   *   3  =>  0  1* 2  3  5 
   *   4  =>  0* 1  2  3  5 
   *   5  =>  0  1  2  4  5 
   *   6  =>  0  1  2  4* 5 
   *   7  =>  0  1  2* 4  5 
   *   8  =>  0  1* 2  4  5 
   *   9  =>  0* 1  2  4  5 
   *  10  =>  0  1  3  4  5 
   *  11  =>  0  1  3  4* 5 
   *  12  =>  0  1  3* 4  5 
   *  13  =>  0  1* 3  4  5 
   *  14  =>  0* 1  3  4  5 
   *  15  =>  0  2  3  4  5 
   *  16  =>  0  2  3  4* 5 
   *  17  =>  0  2  3* 4  5 
   *  18  =>  0  2* 3  4  5 
   *  19  =>  0* 2  3  4  5 
   *  20  =>  1  2  3  4  5 
   *  21  =>  1  2  3  4* 5  ?? and  1  2  3  5 ?? 
   *  22  =>  1  2  3* 4  5  ?? and  1  2  4  5 ?? 
   *  23  =>  1  2* 3  4  5  ?? and  1  3  4  5 ?? 
   *  24  =>  1* 2  3  4  5  ?? and  2  3  4  5 ?? 
   *  25  =>  ?? 0  1  2  5 ?? 
   *  26  =>  ?? 0  1  3  5 ?? 
   *  27  =>  ?? 0  1  4  5 ?? 
   *  28  =>  ?? 0  2  3  5 ?? 
   *  29  =>  ?? 0  2  4  5 ?? 
   *  30  =>  ?? 0  3  4  5 ?? 
   *  31  =>  invalid/reserved/undecided/dunno 
   */

  /* Key for TF upgrade - this is all we've got! 
   * Just make sure to select if it's 5/5 that we  
   * make a good choice for which 4 layers to use 

   0   =>  0  1  2  3  5 
   1   =>  0  1  2  4  5 
   2   =>  0  1  3  4  5 
   3   =>  0  2  3  4  5 
   4   =>  1  2  3  4  5 */

  /*printf("in svtsim_whichFit_full\n"); 
     printf("zin = %.6x, layerMask = %.6x, lcmask = %.6x \n",zin,layerMask,lcMask); 
   */
  int bogus = 0, l = 0, m = 0;
  int lcUsedMask = 0;
  int nsvx = 0, ngoodsvx = 0, firstlc = 4;
  if (zin < 0 || zin >= SVTSIM_NBAR)
    zin = 0;
  for (l = 0; l < 5; l++)
  {
    if (layerMask >> l & 1)
    {
      nsvx++;
      if (lcMask >> l & 1)
      {
        lcUsedMask |= 1 << m;
        if (firstlc > m)
          firstlc = m;
      }
      else
      {
        ngoodsvx++;
      }
      m++;
    }
  }
  switch (layerMask & 0x1f)
  {
  case 0x0f:                   /* 0123 */
    return 0;
  case 0x17:                   /* 0124 */
    return 1;
  case 0x1b:                   /* 0134 */
    return 2;
  case 0x1d:                   /* 0234 */
    return 3;
  case 0x1e:                   /* 1234 */
    return 4;
  case 0x1f:                   /* 01234 - this is the fun one to be careful with */
    if (lcMask == 0)
      return 2;                 /* use 0134 if we have no LC */

    else if (lcMask == 0x1)
      return 4;

    else if (lcMask == 0x2)
      return 3;

    else if (lcMask == 0x3)
      return 3;

    else if (lcMask == 0x4)
      return 2;

    else if (lcMask == 0x5)
      return 2;

    else if (lcMask == 0x6)
      return 2;

    else if (lcMask == 0x7)
      return 2;

    else if (lcMask == 0x8)
      return 1;

    else if (lcMask == 0x9)
      return 1;

    else if (lcMask == 0xa)
      return 1;

    else if (lcMask == 0xb)
      return 1;

    else if (lcMask == 0xc)
      return 2;

    else if (lcMask == 0xd)
      return 2;

    else if (lcMask == 0xe)
      return 2;

    else if (lcMask == 0xf)
      return 2;

    else                        /* If we have LC on outer layer just use 0123 */
      return 0;
  default:
    return bogus;
  }
}


/* 
 * Map layerMask, lcMask into fit block (0..31), allowing for 
 * degeneracy in TF mkaddr map 
 */
int
svtsim_whichFit (tf_arrays_t tf, int zin, int layerMask, int lcMask)
{

#ifdef DEBUG_SVT
  printf ("in svtsim_whichFit: zin = %d, layerMask = %x, lcMask = %x\n",
          zin, layerMask, lcMask);

#endif /*   */
  int which0 = 0, which = 0;
  if (zin < 0 || zin >= SVTSIM_NBAR)
    zin = 0;
  which0 = svtsim_whichFit_full (tf, zin, layerMask, lcMask);
  which = tf->whichFit[zin][which0];

#ifdef DEBUG_SVT
  printf ("in svtsim_whichFit: which0 = %d, which = %x\n", which0, which);

#endif /*   */
  return which;
}


int
svtsim_fconread (tf_arrays_t tf)
{
  int i = 0, j = 0, k = 0, which = 0;
  char fconFile[SVTSIM_NBAR][100];
  char fcf[SVTSIM_NBAR][100];
  char str[100];

  //initialization of variables 
  tf->mkaddrBogusValue = 1000;

  /* 
     from svtsim_maps.c 
     //  value = svtsim_dict_query(p->mapSetDict, "phiOffsetBins"); 
     //p->dphiDenom = value ? atoi(value) : SVTSIM_NWEDGE; 
     //  sprintf(str, "w%2.2d_phiOffset", wedge); 
     // value = svtsim_dict_query(p->mapSetDict, str); 
     //p->dphiNumer = value ? atoi(value) : wedge; 
   */
  tf->dphiDenom = 12;
  tf->dphiNumer = WEDGE;
  for (i = 0; i < NEVTS; i++)
  {
    tf->wedge[i] = WEDGE;
  }
  for (i = 0; i < SVTSIM_NBAR; i++)
  {
    sprintf (str, "offln_220050_20100810_w%dz%d.fcon", WEDGE, i);
    strncpy (fcf[i], str, sizeof (fcf[i]) - 1);

    //     printf("fcf[%d] = %s\n", i,  fcf[i]); 
  }

  /* 
     Get layer map 
   */
  for (i = 0; i < SVTSIM_NBAR; i++)
  {
    int j = 0;
    for (j = 0; j < SVTSIM_NPL; j++)
    {

      //       int defMap[] = { -1, 0, 1, 2, 3, -1, 5 }; 
      int defMap[] = { -1, 0, 1, 2, 3, 4, 5 };
      tf->lMap[i][j] = defMap[j];

#ifdef DEBUG_READFCON
      printf ("layer map --> tf->lMap[%d][%d]=%d\n", i, j, tf->lMap[i][j]);

#endif /*   */
  }}
  /* 
   * Read fit constants for all barrels 
   */
  for (k = 0; k < SVTSIM_NBAR; k++)
  {
    int n = 0, version = 0;
    int xorg = 0, yorg = 0;
    float x = 0, ftmp = 0;
    char str[80], s[80];
    int h[DIMSPA] = { 0, 0, 0, 0, 0, 0 };
    float dx[DIMSPA] = { 0, 0, 0, 0, 0, 0 };

    /* 
     * Initialize all whichFit entries to bogus values 
     */
    for (i = 0; i < SVTSIM_NEL (tf->whichFit[k]); i++)
      tf->whichFit[k][i] = -1;

    /* 
     * Open .fcon 
     */
    FILE *fp = 0;
    fp = fopen (fcf[k], "r");

    //fp = blobopen(fcf[k], "r"); 
#ifdef DEBUG_READFCON
    printf ("fconread: opening %s for input, barrel %d\n", fcf[k], k);

#endif /*   */
    if (!fp)
    {
      fprintf (stderr, "svtsim: error opening %s for input\n", fcf[k]);
      return -1;
    }

    do
    {
      int l = 0, layerMask = 0, bogusLayers = 0;
      char layers[7] = { 0, 0, 0, 0, 0, 0, 0 };
      if (!fgets (str, sizeof (str), fp))
        break;
      if (strncmp (str, "Version ", strlen ("Version ")))
        continue;
      if ((1 != sscanf (str, "Version %d", &version)) ||
          (version < 3) || (version > 4))
      {
        fprintf (stdout, "svtsim_fconread: %s has wrong format\n", fcf[k]);
        fprintf (stdout, "  expect Version 3 or 4, found %s\n", str);
        fclose (fp);
        return -2;
      }
      ExpectStr ("//NPAR");
      ExpectInt (3);
      ExpectStr ("//Layers");
      Get;
      str[strlen (str) - 1] = 0;
      str[SVTSIM_NEL (layers)] = 0;
      for (l = 0; l < strlen (str); l++)
      {
        int ll = str[l];
        if (ll >= '0' && ll <= '5')
        {
          ll = ll - '0';
        }
        else if (ll == 'F')
        {
          ll = 6;
        }
        else
        {
          ll = -1;
        }
        layers[l] = 'X';
        if (ll >= 0 && ll < SVTSIM_NPL)
        {

#ifdef DEBUG_READFCON
          printf (" str = %s, l = %d, strlen(%s) = %d\n", str, l,
                  str, strlen (str));
          printf ("ll(before assignement) = %d\n", ll);

#endif /*   */
          ll = tf->lMap[k][ll];
          if (ll >= 0)
          {
            layers[l] = '0' + ll;
            layerMask |= 1 << ll;
          }
        }

#ifdef DEBUG_READFCON
        printf ("layers[%d] = %d\n", l, layers[l]);

#endif /*   */
        if (layers[l] == 'X')
          bogusLayers = 1;
      }
      if (bogusLayers)
        continue;
      which = svtsim_whichFit_full (tf, k, layerMask, 0);
      tf->whichFit[k][which] = which;
      if (0)
        printf ("w=%d z=%d: found detLayers=%s => "
                "amLayers=%s mask=%x which=%d\n", WEDGE, k, str, layers,
                layerMask, which);
      ExpectStr ("//Zin,");
      Get;                      /* Zin and longhit_mask data */
      ExpectStr ("//Coordinate");
      Get;
      svtsim_assert (DIMSPA ==
                     sscanf (str, "%f %f %f %f %f %f", dx, dx + 1,
                             dx + 2, dx + 3, dx + 4, dx + 5));
      Get;
      svtsim_assert (DIMSPA ==
                     sscanf (str, "%d %d %d %d %d %d", h, h + 1, h + 2,
                             h + 3, h + 4, h + 5));
      if (version >= 4)
      {
        ExpectStr ("//Origin:");
        Get;
        svtsim_assert (2 == sscanf (str, "%d %d", &xorg, &yorg));
      }
      tf->oX[k] = xorg;
      tf->oY[k] = yorg;
      ExpectStr ("//Parameters:");
      for (i = 0; i < NFITPAR; i++)
      {
        float pstep = 0;
        int fc[DIMSPA + 1] = { 0, 0, 0, 0, 0, 0 };
        int ll[NFITPAR] = { 2, 1, 0, 3, 4, 5 }; /* internal=f,d,c, ext=c,d,f */
        int l = ll[i];
        if (i == 3)
          ExpectStr ("//Kernel");
        Get;
        svtsim_assert (8 ==
                       sscanf (str, "%f %d %d %d %d %d %d %d", &pstep,
                               fc, fc + 1, fc + 2, fc + 3, fc + 4,
                               fc + 5, fc + 6));
        switch (l)
        {
        case 0:
          tf->phiUnit = pstep;
          break;
        case 1:
          tf->dvxUnit = pstep;
          break;
        case 2:
          tf->crvUnit = pstep;
          break;
        case 3:
          tf->k0Unit = pstep;
          break;
        case 4:
          tf->k1Unit = pstep;
          break;
        case 5:
          tf->k2Unit = pstep;
          break;
        }
        for (j = 0; j < DIMSPA; j++)
        {

#ifdef DEBUG_READFCON
          printf ("fc[%d] = %.6x, %.d\n", j, fc[j], fc[j]);

#endif /*   */
          tf->gcon[l][j][k][which] =
            (pstep * fc[j]) / (pow (2.0, h[j]) * dx[j]);
          tf->lfitpar[l][j][k][which] = fc[j];
          tf->lfitpar[l][j][k][which] <<= 30;   /* lfitpar coeffs are <<30 */
          tf->lfitpar[l][j][k][which] = (tf->lfitpar[l][j][k][which] >> h[j]);

          /* fill lfitparfcon array with constants as read from fcon file SA */
          tf->lfitparfcon[l][j][k][which] = fc[j];

#ifdef DEBUG_READFCON
          printf ("lfitparfcon[%d][%d][%d][%d] = %.6x\n", l, j, k,
                  which, tf->lfitparfcon[l][j][k][which]);

#endif /*   */
        }
        tf->gcon[l][DIMSPA][k][which] = fc[DIMSPA] * pstep;
        tf->lfitpar[l][DIMSPA][k][which] = fc[DIMSPA];

        /* fill lfitparfcon array with constants as read from fcon file SA */
        tf->lfitparfcon[l][DIMSPA][k][which] = fc[DIMSPA];

#ifdef DEBUG_READFCON
        printf ("fc[%d] = %.6x, %.d\n", DIMSPA, fc[DIMSPA], fc[DIMSPA]);

#endif /*   */
        if (i == 2)
        {                       /* phi */
          int twopi = floor (0.5 + 2 * M_PI / tf->phiUnit);
          tf->lfitpar[l][DIMSPA][k][which] +=
            twopi * tf->dphiNumer / tf->dphiDenom;

          /* fill lfitparfcon array with constants as read from fcon file SA */
          tf->lfitparfcon[l][DIMSPA][k][which] +=
            twopi * tf->dphiNumer / tf->dphiDenom;
          if (tf->lfitpar[l][DIMSPA][k][which] < 0)
            tf->lfitpar[l][DIMSPA][k][which] += twopi;

          /* fill tf->lfitparfcon array with constants as read from fcon file SA */
          if (tf->lfitparfcon[l][DIMSPA][k][which] < 0)
            tf->lfitparfcon[l][DIMSPA][k][which] += twopi;
          tf->gcon[l][DIMSPA][k][which] +=
            2 * M_PI * tf->dphiNumer / tf->dphiDenom;
        }
        tf->lfitpar[l][DIMSPA][k][which] <<= 16;        /* tf->lfitpar intcp is <<16 */
      }
      GetReal (ftmp);
      svtsim_assert (ftmp == 1.0);      /* chisq scale factor */
    }
    while (1);
    fclose (fp);
  }

  /* 
   * Compute shift values 
   */
  for (i = 0; i < NFITPAR; i++)
  {
    int shftmax_svx = -99, shftmax_phi = -99, shftmax_crv = -99;
    int ishft_svx = 0, ishft_phi = 0, ishft_crv = 0;
    for (k = 0; k < SVTSIM_NBAR; k++)
    {
      for (which = 0; which < SVTSIM_NEL (tf->whichFit[k]); which++)
      {
        int ishft = 0;
        if (tf->whichFit[k][which] != which)
          continue;
        for (j = 0; j < 4; j++)
        {
          ishft = ilog (tf->lfitpar[i][j][k][which]) - 30;
          if (ishft > shftmax_svx)
            shftmax_svx = ishft;
        }
        ishft = ilog (tf->lfitpar[i][4][k][which]) - 30;
        if (ishft > shftmax_phi)
          shftmax_phi = ishft;
        ishft = ilog (tf->lfitpar[i][5][k][which]) - 30;
        if (ishft > shftmax_crv)
          shftmax_crv = ishft;
      }
    }
    ishft_svx = 8 - shftmax_svx;
    if (ishft_svx > 12)
      ishft_svx = 12;
    ishft_phi = shftmax_phi + ishft_svx - 10;
    if (ishft_phi < -3)
      ishft_phi = -3;
    ishft_crv = shftmax_crv + ishft_svx - 10;
    if (ishft_crv < -3)
      ishft_crv = -3;
    tf->result_shiftbits[i] = ishft_svx;
    tf->xftphi_shiftbits[i] = ishft_phi;
    tf->xftcrv_shiftbits[i] = ishft_crv;
  }

  /* 
   * Calculate TF reduced-precision coefficients 
   */
  for (k = 0; k < SVTSIM_NBAR; k++)
  {
    for (i = 0; i < NFITPAR; i++)
    {
      for (j = 0; j < DIMSPA + 1; j++)
      {
        int ishift = 0;
        switch (j)
        {
        case 4:
          ishift = 30 - (tf->result_shiftbits[i] - tf->xftphi_shiftbits[i]);
          break;
        case 5:
          ishift = 30 - (tf->result_shiftbits[i] - tf->xftcrv_shiftbits[i]);
          break;
        case 6:
          ishift = 16;
          break;
        default:
          ishift = 30 - (tf->result_shiftbits[i]);
        }
        for (which = 0; which < SVTSIM_NEL (tf->whichFit[k]); which++)
        {
          if (tf->whichFit[k][which] != which)
            continue;
          tf->ifitpar[i][j][k][which] = tf->lfitpar[i][j][k][which] >> ishift;
        }
      }
    }

    /* 
     * Now propagate existent fit constants into gaps left by  
     * nonexistent fit constants 
     */
    for (i = 0; i <= 0x3f; i++)
    {
      int which0 = -1;
      which = svtsim_whichFit_full (tf, k, i, 0);

      /* 
       * If layer combination has no entry, use a bogus entry 
       */
      if (tf->whichFit[k][which] < 0)
        tf->whichFit[k][which] = tf->mkaddrBogusValue;
      which0 = tf->whichFit[k][which];
      for (j = 0; j <= 0x1f; j++)
      {
        which = svtsim_whichFit_full (tf, k, i, j);

        /* 
         * If long-cluster combination has no entry, use the no-LC entry 
         */
        if (tf->whichFit[k][which] < 0)
          tf->whichFit[k][which] = which0;
      }
    }
    if (0)
    {
      printf ("whichfit(%d):", k);
      for (i = 0; i < SVTSIM_NEL (tf->whichFit[k]); i++)
        printf (" %d", tf->whichFit[k][i]);
      printf ("\n");
    }
  }

#ifdef DEBUG_READFCON
  /*check that constant are correctly saved in the array */
  int nparam = 6;
  int ipar;
  int dimspa = 6;
  int idim;
  int nbarrel = 6;
  int ibar;
  int fitBlock = 5;
  int ifitBlock;
  for (ipar = 0; ipar < nparam; ipar++)
  {
    for (idim = 0; idim < dimspa; idim++)
    {
      for (ibar = 0; ibar < nbarrel; ibar++)
      {
        for (ifitBlock = 0; ifitBlock < fitBlock; ifitBlock++)
        {
          printf
            ("ipar=%d, idim=%d, ibar=%d, ifitBlock=%d, tf->ifitpar = %.6x, lfitparfcon = %.6llx\n",
             ipar, idim, ibar, ifitBlock,
             tf->ifitpar[ipar][idim][ibar][ifitBlock],
             tf->lfitparfcon[ipar][idim][ibar][ifitBlock]);
        }
      }
    }
  }

#endif /*   */
  return 0;
}

int
svtsim_fconread_GPU (tf_arrays_t tf)
{
  int i = 0, j = 0, k = 0, which = 0;
  char fconFile[SVTSIM_NBAR][100];
  char fcf[SVTSIM_NBAR][100];
  char str[100];

  //initialization of variables 
  mkaddrBogusValue = 1977;

  /* 
     from svtsim_maps.c 
     //  value = svtsim_dict_query(p->mapSetDict, "phiOffsetBins"); 
     //p->dphiDenom = value ? atoi(value) : SVTSIM_NWEDGE; 
     //  sprintf(str, "w%2.2d_phiOffset", wedge); 
     // value = svtsim_dict_query(p->mapSetDict, str); 
     //p->dphiNumer = value ? atoi(value) : wedge; 
   */
  dphiDenom = 12;
  dphiNumer = WEDGE;
  for (i = 0; i < NEVTS; i++)
  {
    tf->wedge[i] = WEDGE;
  }
  for (i = 0; i < SVTSIM_NBAR; i++)
  {
    sprintf (str, "offln_220050_20100810_w%dz%d.fcon", WEDGE, i);
    strncpy (fcf[i], str, sizeof (fcf[i]) - 1);

    //     printf("fcf[%d] = %s\n", i,  fcf[i]); 
  }

  /* 
     Get layer map 
   */
  for (i = 0; i < SVTSIM_NBAR; i++)
  {
    int j = 0;
    for (j = 0; j < SVTSIM_NPL; j++)
    {

      //       int defMap[] = { -1, 0, 1, 2, 3, -1, 5 }; 
      int defMap[] = { -1, 0, 1, 2, 3, 4, 5 };
      lMap[i][j] = defMap[j];

#ifdef DEBUG_READFCON
      printf ("layer map -->  lMap[%d][%d]=%d\n", i, j, lMap[i][j]);

#endif /*   */
  }}
  /* 
   * Read fit constants for all barrels 
   */
  for (k = 0; k < SVTSIM_NBAR; k++)
  {
    int n = 0, version = 0;
    int xorg = 0, yorg = 0;
    float x = 0, ftmp = 0;
    char str[80], s[80];
    int h[DIMSPA] = { 0, 0, 0, 0, 0, 0 };
    float dx[DIMSPA] = { 0, 0, 0, 0, 0, 0 };

    /* 
     * Initialize all whichFit entries to bogus values 
     */
    for (i = 0; i < SVTSIM_NEL (whichFit[k]); i++)
      whichFit[k][i] = -1;

    /* 
     * Open .fcon 
     */
    FILE *fp = 0;
    fp = fopen (fcf[k], "r");

    //fp = blobopen(fcf[k], "r"); 
#ifdef DEBUG_READFCON
    printf ("fconread: opening %s for input, barrel %d\n", fcf[k], k);

#endif /*   */
    if (!fp)
    {
      fprintf (stderr, "svtsim: error opening %s for input\n", fcf[k]);
      return -1;
    }

    do
    {
      int l = 0, layerMask = 0, bogusLayers = 0;
      char layers[7] = { 0, 0, 0, 0, 0, 0, 0 };
      if (!fgets (str, sizeof (str), fp))
        break;
      if (strncmp (str, "Version ", strlen ("Version ")))
        continue;
      if ((1 != sscanf (str, "Version %d", &version)) ||
          (version < 3) || (version > 4))
      {
        fprintf (stdout, "svtsim_fconread: %s has wrong format\n", fcf[k]);
        fprintf (stdout, "  expect Version 3 or 4, found %s\n", str);
        fclose (fp);
        return -2;
      }
      ExpectStr ("//NPAR");
      ExpectInt (3);
      ExpectStr ("//Layers");
      Get;
      str[strlen (str) - 1] = 0;
      str[SVTSIM_NEL (layers)] = 0;
      for (l = 0; l < strlen (str); l++)
      {
        int ll = str[l];
        if (ll >= '0' && ll <= '5')
        {
          ll = ll - '0';
        }
        else if (ll == 'F')
        {
          ll = 6;
        }
        else
        {
          ll = -1;
        }
        layers[l] = 'X';
        if (ll >= 0 && ll < SVTSIM_NPL)
        {

#ifdef DEBUG_READFCON
          printf (" str = %s, l = %d, strlen(%s) = %d\n", str, l,
                  str, strlen (str));
          printf ("ll(before assignement) = %d\n", ll);

#endif /*   */
          ll = lMap[k][ll];
          if (ll >= 0)
          {
            layers[l] = '0' + ll;
            layerMask |= 1 << ll;
          }
        }

#ifdef DEBUG_READFCON
        printf ("layers[%d] = %d\n", l, layers[l]);

#endif /*   */
        if (layers[l] == 'X')
          bogusLayers = 1;
      }
      if (bogusLayers)
        continue;
      which = svtsim_whichFit_full (tf, k, layerMask, 0);
      whichFit[k][which] = which;
      if (0)
        printf ("w=%d z=%d: found detLayers=%s => "
                "amLayers=%s mask=%x which=%d\n", WEDGE, k, str, layers,
                layerMask, which);
      ExpectStr ("//Zin,");
      Get;                      /* Zin and longhit_mask data */
      ExpectStr ("//Coordinate");
      Get;
      svtsim_assert (DIMSPA ==
                     sscanf (str, "%f %f %f %f %f %f", dx, dx + 1,
                             dx + 2, dx + 3, dx + 4, dx + 5));
      Get;
      svtsim_assert (DIMSPA ==
                     sscanf (str, "%d %d %d %d %d %d", h, h + 1, h + 2,
                             h + 3, h + 4, h + 5));
      if (version >= 4)
      {
        ExpectStr ("//Origin:");
        Get;
        svtsim_assert (2 == sscanf (str, "%d %d", &xorg, &yorg));
      }
      oX[k] = xorg;
      oY[k] = yorg;
      ExpectStr ("//Parameters:");
      for (i = 0; i < NFITPAR; i++)
      {
        float pstep = 0;
        int fc[DIMSPA + 1] = { 0, 0, 0, 0, 0, 0 };
        int ll[NFITPAR] = { 2, 1, 0, 3, 4, 5 }; /* internal=f,d,c, ext=c,d,f */
        int l = ll[i];
        if (i == 3)
          ExpectStr ("//Kernel");
        Get;
        svtsim_assert (8 ==
                       sscanf (str, "%f %d %d %d %d %d %d %d", &pstep,
                               fc, fc + 1, fc + 2, fc + 3, fc + 4,
                               fc + 5, fc + 6));
        switch (l)
        {
        case 0:
          phiUnit = pstep;
          break;
        case 1:
          dvxUnit = pstep;
          break;
        case 2:
          crvUnit = pstep;
          break;
        case 3:
          k0Unit = pstep;
          break;
        case 4:
          k1Unit = pstep;
          break;
        case 5:
          k2Unit = pstep;
          break;
        }
        for (j = 0; j < DIMSPA; j++)
        {

#ifdef DEBUG_READFCON
          printf ("fc[%d] = %.6x, %.d\n", j, fc[j], fc[j]);

#endif /*   */
          gcon[l][j][k][which] = (pstep * fc[j]) / (pow (2.0, h[j]) * dx[j]);
          lfitpar[l][j][k][which] = fc[j];
          lfitpar[l][j][k][which] <<= 30;       /* lfitpar coeffs are <<30 */
          lfitpar[l][j][k][which] = (lfitpar[l][j][k][which] >> h[j]);

          /* fill lfitparfcon array with constants as read from fcon file SA */
          lfitparfcon[l][j][k][which] = fc[j];

#ifdef DEBUG_READFCON
          printf ("lfitparfcon[%d][%d][%d][%d] = %.6x\n", l, j, k,
                  which, lfitparfcon[l][j][k][which]);

#endif /*   */
        }
        gcon[l][DIMSPA][k][which] = fc[DIMSPA] * pstep;
        lfitpar[l][DIMSPA][k][which] = fc[DIMSPA];

        /* fill lfitparfcon array with constants as read from fcon file SA */
        lfitparfcon[l][DIMSPA][k][which] = fc[DIMSPA];

#ifdef DEBUG_READFCON
        printf ("fc[%d] = %.6x, %.d\n", DIMSPA, fc[DIMSPA], fc[DIMSPA]);

#endif /*   */
        if (i == 2)
        {                       /* phi */
          int twopi = floor (0.5 + 2 * M_PI / phiUnit);
          lfitpar[l][DIMSPA][k][which] += twopi * dphiNumer / dphiDenom;

          /* fill lfitparfcon array with constants as read from fcon file SA */
          lfitparfcon[l][DIMSPA][k][which] += twopi * dphiNumer / dphiDenom;
          if (lfitpar[l][DIMSPA][k][which] < 0)
            lfitpar[l][DIMSPA][k][which] += twopi;

          /* fill  lfitparfcon array with constants as read from fcon file SA */
          if (lfitparfcon[l][DIMSPA][k][which] < 0)
            lfitparfcon[l][DIMSPA][k][which] += twopi;
          gcon[l][DIMSPA][k][which] += 2 * M_PI * dphiNumer / dphiDenom;
        }
        lfitpar[l][DIMSPA][k][which] <<= 16;    /*  lfitpar intcp is <<16 */
      }
      GetReal (ftmp);
      svtsim_assert (ftmp == 1.0);      /* chisq scale factor */
    }
    while (1);
    fclose (fp);
  }

  /* 
   * Compute shift values 
   */
  for (i = 0; i < NFITPAR; i++)
  {
    int shftmax_svx = -99, shftmax_phi = -99, shftmax_crv = -99;
    int ishft_svx = 0, ishft_phi = 0, ishft_crv = 0;
    for (k = 0; k < SVTSIM_NBAR; k++)
    {
      for (which = 0; which < SVTSIM_NEL (whichFit[k]); which++)
      {
        int ishft = 0;
        if (whichFit[k][which] != which)
          continue;
        for (j = 0; j < 4; j++)
        {
          ishft = ilog (lfitpar[i][j][k][which]) - 30;
          if (ishft > shftmax_svx)
            shftmax_svx = ishft;
        }
        ishft = ilog (lfitpar[i][4][k][which]) - 30;
        if (ishft > shftmax_phi)
          shftmax_phi = ishft;
        ishft = ilog (lfitpar[i][5][k][which]) - 30;
        if (ishft > shftmax_crv)
          shftmax_crv = ishft;
      }
    }
    ishft_svx = 8 - shftmax_svx;
    if (ishft_svx > 12)
      ishft_svx = 12;
    ishft_phi = shftmax_phi + ishft_svx - 10;
    if (ishft_phi < -3)
      ishft_phi = -3;
    ishft_crv = shftmax_crv + ishft_svx - 10;
    if (ishft_crv < -3)
      ishft_crv = -3;
    result_shiftbits[i] = ishft_svx;
    xftphi_shiftbits[i] = ishft_phi;
    xftcrv_shiftbits[i] = ishft_crv;
  }

  /* 
   * Calculate TF reduced-precision coefficients 
   */
  for (k = 0; k < SVTSIM_NBAR; k++)
  {
    for (i = 0; i < NFITPAR; i++)
    {
      for (j = 0; j < DIMSPA + 1; j++)
      {
        int ishift = 0;
        switch (j)
        {
        case 4:
          ishift = 30 - (result_shiftbits[i] - xftphi_shiftbits[i]);
          break;
        case 5:
          ishift = 30 - (result_shiftbits[i] - xftcrv_shiftbits[i]);
          break;
        case 6:
          ishift = 16;
          break;
        default:
          ishift = 30 - (result_shiftbits[i]);
        }
        for (which = 0; which < SVTSIM_NEL (whichFit[k]); which++)
        {
          if (whichFit[k][which] != which)
            continue;
          ifitpar[i][j][k][which] = lfitpar[i][j][k][which] >> ishift;
        }
      }
    }

    /* 
     * Now propagate existent fit constants into gaps left by  
     * nonexistent fit constants 
     */
    for (i = 0; i <= 0x3f; i++)
    {
      int which0 = -1;
      which = svtsim_whichFit_full (tf, k, i, 0);

      /* 
       * If layer combination has no entry, use a bogus entry 
       */
      if (whichFit[k][which] < 0)
        whichFit[k][which] = mkaddrBogusValue;
      which0 = whichFit[k][which];
      for (j = 0; j <= 0x1f; j++)
      {
        which = svtsim_whichFit_full (tf, k, i, j);

        /* 
         * If long-cluster combination has no entry, use the no-LC entry 
         */
        if (whichFit[k][which] < 0)
          whichFit[k][which] = which0;
      }
    }
    if (0)
    {
      printf ("whichfit(%d):", k);
      for (i = 0; i < SVTSIM_NEL (whichFit[k]); i++)
        printf (" %d", whichFit[k][i]);
      printf ("\n");
    }
  }

#ifdef DEBUG_READFCON
  /*check that constant are correctly saved in the array */
  int nparam = 6;
  int ipar;
  int dimspa = 6;
  int idim;
  int nbarrel = 6;
  int ibar;
  int fitBlock = 5;
  int ifitBlock;
  for (ipar = 0; ipar < nparam; ipar++)
  {
    for (idim = 0; idim < dimspa; idim++)
    {
      for (ibar = 0; ibar < nbarrel; ibar++)
      {
        for (ifitBlock = 0; ifitBlock < fitBlock; ifitBlock++)
        {
          printf
            ("ipar=%d, idim=%d, ibar=%d, ifitBlock=%d,  ifitpar = %.6x, lfitparfcon = %.6llx\n",
             ipar, idim, ibar, ifitBlock,
             ifitpar[ipar][idim][ibar][ifitBlock],
             lfitparfcon[ipar][idim][ibar][ifitBlock]);
        }
      }
    }
  }

#endif /*   */
  return 0;
}


void
svtsim_assert_set (char *filename, int line)
{
  /*char buf[254]; 
     sprintf(buf,"SVTSIMMODULE::ERROR %s: %d\n",filename,line); 
     if(strlen(buf)+strlen(svtsim_err_str)<svtsim_str_size){ 
     strcat(svtsim_err_str,buf); 
     } */
  if (svtsim_n_err < 10)
    sprintf (svtsim_err_str[svtsim_n_err], "SVTSIM::ASSERT %s: %d\n",
             filename, line);
  svtsim_n_err++;
}

int
svtsim_memDebug (int nChainDump)
{
  int i = 0, t = 0, tb = 0;
  memDebug_t *p = chain_alloc;
  fprintf (stderr, "svtsim_memDebug: n_alloc=%d t_alloc=%d chain+1=%p\n",
           n_alloc, t_alloc, chain_alloc + 1);
  for (; p; p = p->next, i++)
  {
    svtsim_assert (p->deadBeef == 0xdeadBeef);
    svtsim_assert (p->next == 0 || p->next->last == p);
    t++;
    tb += p->userBytes;
    if (p->lineNum < 0)
      continue;                 /* allocated "forever" */
    if (nChainDump < 0 || i < nChainDump)
      fprintf (stderr, "  p=%p sz=%d line=%s:%d\n", p + 1, p->userBytes,
               p->fileName, p->lineNum);
  }
  svtsim_assert (t == n_alloc);
  svtsim_assert (tb == t_alloc);
  return t_alloc;
}

void *
svtsim_malloc1 (size_t size, const char *filename, int linenum)
{
  memDebug_t *p = 0;
  p = (void *) malloc (size + sizeof (*p));
  if (p == NULL)
    printf ("FAILED MALLOC %s %d!\n", filename, linenum);
  svtsim_assert (p != NULL);
  p->deadBeef = 0xdeadBeef;
  p->userBytes = size;
  p->fileName = filename;
  p->lineNum = linenum;
  p->last = 0;
  p->next = chain_alloc;
  chain_alloc = p;
  if (p->next)
  {
    svtsim_assert (p->next->last == 0);
    p->next->last = p;
  }
  memset (p + 1, 0, size);
  n_alloc++;
  t_alloc += p->userBytes;
  if (t_alloc > m_alloc)
    m_alloc = t_alloc;

#if 0
  for (i = 0; i < SVTSIM_NEL (badaddr); i++)
    if (p + 1 == badaddr[i])
    {
      svtsim_memDebug (2);
      svtsim_assert (0);
    }

#endif /*   */
  return ((void *) (p + 1));
} void

svtsim_free1 (void *p, const char *filename, int linenum)
{
  int nbytes = 0;
  memDebug_t *q = ((memDebug_t *) p) - 1;
  if (!p)
    return;
  svtsim_assert (p != (void *) 0xffffffff);
  if (q->deadBeef != 0xdeadBeef)
  {
    fprintf (stderr, "%p->deadBeef==0x%x (%s:%d)\n", q, q->deadBeef,
             filename, linenum);
    free (p);
    return;
  }
  svtsim_assert (q->deadBeef == 0xdeadBeef);
  svtsim_assert (q->lineNum >= 0);      /* don't free "forever" mallocs */
  if (q->last)
  {
    svtsim_assert (q->last->next == q);
    q->last->next = q->next;
  }
  else
  {
    svtsim_assert (chain_alloc == q);
    chain_alloc = q->next;
  }
  if (q->next)
  {
    svtsim_assert (q->next->last == q);
    q->next->last = q->last;
  }
  n_alloc--;
  t_alloc -= q->userBytes;
  nbytes = sizeof (*q) + q->userBytes;
  memset (q, -1, nbytes);
  free (q);
}

void *
svtsim_realloc1 (void *inptr, size_t size, const char *filename, int linenum)
{
  memDebug_t *p = ((memDebug_t *) inptr) - 1;
  if (!inptr)
    return svtsim_malloc1 (size, filename, linenum);
  if (!size)
  {
    svtsim_free1 (inptr, filename, linenum);
    return 0;
  }
  if (0)
  {                             /* debug */
    svtsim_free1 (inptr, filename, linenum);
    return svtsim_malloc1 (size, filename, linenum);
  }
  svtsim_assert (p->deadBeef == 0xdeadBeef);
  if (p->last)
  {
    svtsim_assert (p->last->next == p);
  }
  else
  {
    svtsim_assert (p == chain_alloc);
  }
  if (p->next)
    svtsim_assert (p->next->last == p);
  t_alloc -= p->userBytes;
  p = (memDebug_t *) realloc (p, size + sizeof (*p));
  if (p->last)
    p->last->next = p;
  else
    chain_alloc = p;
  if (p->next)
    p->next->last = p;
  p->userBytes = size;
  t_alloc += p->userBytes;
  if (t_alloc > m_alloc)
    m_alloc = t_alloc;
  svtsim_assert (p != 0);

#if 0
  for (i = 0; i < SVTSIM_NEL (badaddr); i++)
    if (p + 1 == badaddr[i])
    {
      svtsim_memDebug (2);
      svtsim_assert (0);
    }

#endif /*   */
  return p + 1;
}
