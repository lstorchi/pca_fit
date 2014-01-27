#include "svtsim_functions.h"

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


  /* 
   * Default BLOB loader: read file from online directory 
   */
static int
getBlobFromFile (const char *name, char **datap, void *extra)
{
  int dlen = 0;
  FILE *fp = 0;
  char str[256];
  char *dlPath = (char *) extra;
  if (!dlPath)
    dlPath = "/cdf/onln/code/cdfprod/svt_config";
  *datap = 0;
  sprintf (str, "%s/%s", dlPath, name);
  svtsim_assert (strlen (str) < sizeof (str));
  fp = fopen (str, "r");
  if (!fp)
  {
    fprintf (stderr, "svtsim::getBlobFromFile: can't open %s\n", str);
    return -2;
  }
  svtsim_assert (0 == fseek (fp, 0, SEEK_END));
  dlen = ftell (fp);
  svtsim_assert (0 == fseek (fp, 0, SEEK_SET));
  *datap = svtsim_malloc (dlen + 1);
  svtsim_assert (1 == fread (*datap, dlen, 1, fp));
  (*datap)[dlen] = 0;
  fclose (fp);
  return dlen;
}


  /* 
   * Pointer to function to fetch a BLOB from the database; also keep 
   * a pointer to some data that the fetch function may need (e.g. data path) 
   */
static int (*pGetBlob) (const char *name, char **datap, void *extra) =
  getBlobFromFile;
static void *pGetBlobExtra = 0;
typedef struct blob_s
{
  char *data;
  int dlen;
  int offset;
} blob_t;

  /* 
   * Fetch a BLOB from the database: return length>=0 or error<0 
   */
static int
getBlob (const char *name, char **datap)
{
  int rc = -1;
  if (pGetBlob)
    rc = (*pGetBlob) (name, datap, pGetBlobExtra);
  if (1)
    fprintf (stderr, "getBlob(%s) => %d\n", name, rc);
  return rc;
}

static blob_t *
blobopen (const char *name, const char *dummy)
{
  blob_t *b = 0;
  b = svtsim_malloc (sizeof (*b));
  b->offset = 0;
  b->dlen = getBlob (name, &b->data);
  if (b->dlen < 0)
  {
    svtsim_free (b);
    return 0;
  }
  return b;
}

static char *
blobgets (char *s, int size, blob_t * b)
{
  printf ("in blobgets: s = %s\n", s);
  int np = 0;
  const char *p = 0;
  if (!b || !b->data || b->offset >= b->dlen)
    return 0;
  p = b->data + b->offset;
  while (b->offset + np < b->dlen && np < size - 1)
    if (p[np++] == '\n')
      break;
  b->offset += np;
  svtsim_assert (np < size);
  strncpy (s, p, np);
  s[np] = 0;
  return s;
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


  /* 
   * Integer binary logarithm (ceil(log(fabs(ix+epsilon))/log(2))) 
   */
static int
ilog (long long ix)
{
  long long one = 1;
  if (ix == one << 63)
    return 64;
  if (ix < 0)
    ix = -ix;

  /* 
   * Oops, VxWorks compiler won't shift long long by int 
   */
#define IHTFOS(i) if (ix>>(i)&1) return (i)+1
  IHTFOS (62);
  IHTFOS (61);
  IHTFOS (60);
  IHTFOS (59);
  IHTFOS (58);
  IHTFOS (57);
  IHTFOS (56);
  IHTFOS (55);
  IHTFOS (54);
  IHTFOS (53);
  IHTFOS (52);
  IHTFOS (51);
  IHTFOS (50);
  IHTFOS (49);
  IHTFOS (48);
  IHTFOS (47);
  IHTFOS (46);
  IHTFOS (45);
  IHTFOS (44);
  IHTFOS (43);
  IHTFOS (42);
  IHTFOS (41);
  IHTFOS (40);
  IHTFOS (39);
  IHTFOS (38);
  IHTFOS (37);
  IHTFOS (36);
  IHTFOS (35);
  IHTFOS (34);
  IHTFOS (33);
  IHTFOS (32);
  IHTFOS (31);
  IHTFOS (30);
  IHTFOS (29);
  IHTFOS (28);
  IHTFOS (27);
  IHTFOS (26);
  IHTFOS (25);
  IHTFOS (24);
  IHTFOS (23);
  IHTFOS (22);
  IHTFOS (21);
  IHTFOS (20);
  IHTFOS (19);
  IHTFOS (18);
  IHTFOS (17);
  IHTFOS (16);
  IHTFOS (15);
  IHTFOS (14);
  IHTFOS (13);
  IHTFOS (12);
  IHTFOS (11);
  IHTFOS (10);
  IHTFOS (9);
  IHTFOS (8);
  IHTFOS (7);
  IHTFOS (6);
  IHTFOS (5);
  IHTFOS (4);
  IHTFOS (3);
  IHTFOS (2);
  IHTFOS (1);
  IHTFOS (0);

#undef IHTFOS
  return 0;
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
} void

svtsim_cable_copywords (svtsim_cable_t * cable, unsigned int *word, int nword)
{
  svtsim_assert (cable);
  cable->ndata = 0;
  svtsim_cable_addwords (cable, word, nword);
} int

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
