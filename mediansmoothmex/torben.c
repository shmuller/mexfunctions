/*
 * The following code is public domain.
 * Algorithm by Torben Mogensen, implementation by N. Devillard.
 * This code in public domain.
 */
elem_type torben(const elem_type m[], const int n)
{
    int i, less, greater, equal;
    elem_type min, max, guess, maxltguess, mingtguess;
    min = max = m[0];
    for (i=1; i<n; i++) {
        if (m[i] < min) min = m[i];
        if (m[i] > max) max = m[i];
    }
    while (1) {
        guess = (min + max) / 2;
        less = greater = equal = 0;
        maxltguess = min;
        mingtguess = max;
        for (i=0; i<n; i++) {
            if (m[i] < guess) {
                less++;
                if (m[i] > maxltguess) maxltguess = m[i];
            } else if (m[i] > guess) {
                greater++;
                if (m[i] < mingtguess) mingtguess = m[i];
            } else equal++;
        }
        if (less <= (n+1)/2 && greater <= (n+1)/2) break;
        else if (less > greater) max = maxltguess;
        else min = mingtguess;
    }
    if (less >= (n+1)/2) return maxltguess;
    else if (less+equal >= (n+1)/2) return guess;
    else return mingtguess;
}
