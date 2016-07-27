#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char ** argv)
{
  if (argc != 2)
  {
    fprintf (stderr, "usage: %s eta \n", argv[0]);
    return 1;
  }

  double eta = atof(argv[1]);
  double val = 2.0 * atan (exp (-1.0e0 * eta));
  double cott = cos(val)/sin(val);

  fprintf (stdout, "cot(theta): %f \n", cott);

  return 0;
}
