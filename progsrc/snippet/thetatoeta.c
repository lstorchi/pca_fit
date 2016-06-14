#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main (int argc, char ** argv)
{
  if (argc != 2)
  {
    fprintf (stderr, "usage: %s cot(theta) \n", argv[0]);
    return 1;
  }

  double theta = atan(1.0e0 / atof(argv[1]));
  double etaorig = 0.0e0;
  double tantheta2 = tan (theta/2.0e0);
  if (tantheta2 < 0.0)
    etaorig = 1.0e0 * log (-1.0e0 * tantheta2);
  else
    etaorig = -1.0e0 * log (tantheta2);

  fprintf (stdout, "eta: %f \n", etaorig);

  return 0;
}
