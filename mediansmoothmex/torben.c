/*
 * The following code is public domain.
 * Algorithm by Torben Mogensen, implementation by N. Devillard.
 * This code in public domain.
 */
int torben(const void *x, const int n)
{
    const elem_type *m=x;
    elem_type guess;
    int i, less, greater, equal, imin, imax, iequal, imaxltguess, imingtguess;
    
    imin = imax = 0;
    for (i=1; i<n; i++) {
        if (m[i] < m[imin]) imin = i;
        if (m[i] > m[imax]) imax = i;
    }
    while (1) {
        guess = (m[imin] + m[imax]) / 2;
        less = greater = equal = 0;
        imaxltguess = imin;
        imingtguess = imax;
        for (i=0; i<n; i++) {
            if (m[i] < guess) {
                less++;
                if (m[i] > m[imaxltguess]) imaxltguess = i;
            } else if (m[i] > guess) {
                greater++;
                if (m[i] < m[imingtguess]) imingtguess = i;
            } else {
                equal++;
                iequal = i;
            }
        }
        if (less <= (n+1)/2 && greater <= (n+1)/2) break;
        else if (less > greater) imax = imaxltguess;
        else imin = imingtguess;
    }
    if (less >= (n+1)/2) return imaxltguess;
    else if (less+equal >= (n+1)/2) return iequal;
    else return imingtguess;
}
