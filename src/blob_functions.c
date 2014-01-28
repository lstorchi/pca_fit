#include "private_includes.h"

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
