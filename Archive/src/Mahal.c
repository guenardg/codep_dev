/*************************************************************************
 
 (c) 2008-2020 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Mahalanobis Distance**
 
 This file is part of codep
 
 codep is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 codep is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with codep.  If not, see <https://www.gnu.org/licenses/>.
 
 C fonctions definitions
 
 *************************************************************************/

#include"Mahal.h"

/* Computes the Mahalanobis distance between two matrices.
 That function can take two matrices with identical numbers
 of columns and output a rectangular
 distance matrix.
 */

void dist_Mahal(double* from, double* to, int* n, int* tri, int* m,
                double* d,
                double* cov, int* sq) {
  int end, i, j, k, l, ii, jj;
  double acc;
  end = (*tri) ? n[0]*(n[0]-1)/2 : n[0]*n[1];
  for(i = 0, j = (*tri) ? 1 : 0, k = 0; k < end; k++) {
    for(l = 0, ii = i, jj = j; l < *m; l++,
        ii += n[0], jj += n[1]) {
      //
      // mod from here...
      acc = from[ii] - to[jj];
      acc *= acc;
      d[k] += acc;
      // mod to there.
      //
    }
    if(!*sq)
      d[k] = sqrt(d[k]);
    j++;
    if(j == n[1]) {
      i++;
      j = (*tri) ? i+1 : 0;
    }
    if(!(k%INTERRUPT_CHECK))
      R_CheckUserInterrupt();
  }
  return;
}

