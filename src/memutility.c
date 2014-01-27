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
} 

void
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


