int mv1 ( int ns, int *fn, int nVar, int nPts, int *ind, int *cnt, int *sup );
/*
 * DESCRITPION :
 *   Computes the mixed volume of the polytopes spanned by given supports.
 *
 * ON ENTRY :
 *   ns      number of 'characters' in the file name fn;
 *   fn      name of the output file for the mixed-cell configuration;
 *   nVar    ambient dimension, length of the vectors in supports;
 *   nPts    total number of points in the supports;
 *   ind     ind(i) is the start of the i-th support;
 *   cnt     cnt(i) counts the length of the i-th support;
 *   sup     coordinates of the points in the supports. */
