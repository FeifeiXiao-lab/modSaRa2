/* diagnosticValue.c -- get an array of local diagnostic function value */

void diagnosticValue(double *y, int *h, int *n, double *z)
{
    int i, j;
    int hh = *h, nn = *n;

    for(i = 0; i < hh; i++)
        z[0] += (y[i] - y[i + hh])/(1.0 * hh);

    for(j = 1; j < nn; j++)
        z[j] = z[j - 1] + (2.0 * y[j - 1 + hh] - y[j - 1]
               - y[j + 2 * hh - 1])/(1.0 * hh);
}
