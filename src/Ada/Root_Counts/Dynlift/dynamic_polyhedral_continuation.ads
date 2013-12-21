with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Dynamic_Mixed_Subdivisions;         use Dynamic_Mixed_Subdivisions;

package Dynamic_Polyhedral_Continuation is

-- DESCRIPTION :
--   This package contains the utilities for polyhedral continuation,
--   to apply in combination with the dynamic lifting algorithm.

  procedure Dynamic_Unmixed_Solve
                ( file : in file_type; n : in integer32;
                  L : in List; order,inter : in boolean; maxli : in integer32;
                  lifted,lifted_last : in out List; t : in out Triangulation;
                  q : in Poly_Sys; qsols : in out Solution_List );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the unmixed
  --   polynomial system q.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   n          dimension before lifting;
  --   l          list of integer vectors of length n;
  --   order      if true, then the points are already ordered and will
  --              be processed like they occur in the list,
  --              if false, then a random order will be chosen;
  --   inter      if true, then the list may contain interior points,
  --              i.e. points x in conv(l\{x}),
  --              if false, no interior points occur in the list;
  --   maxli      maximum value of the lifting function,
  --              if = 0, then no flattening will be applied,
  --              if > 0, then no points will be given a lifting value
  --              larger than maxli, eventually flattening will be applied;
  --   lifted     points that already have been lifted;
  --   lifted_last is pointer to the last element of the list lifted;
  --   t          contains a triangulation of the lifted points;
  --   q          a polynomial system with randomly chosen coefficients,
  --              the lists of exponent vectors of q equals the list l.

  -- ON RETURN :
  --   lifted     the lifted points;
  --   t          a regular triangulation of the points in l;
  --   qsols      the solutions of q.

  procedure Dynamic_Cayley_Solve
                ( file : in file_type; n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter : in boolean;
                  maxli : in integer32; lifted : in out Array_of_Lists;
                  mixsub : in out Mixed_Subdivision; numtri : out natural32;
                  q : in Poly_Sys; qsols : in out Solution_List );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to the Cayley trick.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   n          length of the vectors in the supports;
  --   mix        type of mixture;
  --   supports   lists of integer vector of length n;
  --   order      if true, then the points are already ordered and will
  --              be processed like they occur in the list,
  --              if false, then a random order will be chosen;
  --   inter      if true, then the list may contain interior points,
  --              i.e. points x in conv(l\{x}),
  --              if false, no interior points occur in the list;
  --   maxli      maximum value of the lifting function,
  --              if = 0, then no flattening will be applied,
  --              i.e. there is no maximum lifting value,
  --              if > 0, then no points will be given a lifting value
  --              greater than maxli;
  --   q          a polynomial system with randomly chosen coefficients,
  --              the lists of exponent vectors of q equals the given supports.

  -- ON RETURN :
  --   lifted     the lifted supports;
  --   mixsub     the mixed subdivision;
  --   numtri     number of simplices in the triangulation of the
  --              auxiliary polytope,
  --   qsols      the solutions of q.

  procedure Dynamic_Mixed_Solve
                ( file : in file_type; n : in integer32; mix : in Vector;
                  supports : in Array_of_Lists; order,inter,conmv : in boolean;
                  maxli : in integer32; mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures; 
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                  q : in Poly_Sys; qsols : in out Solution_List );

  -- DESCRIPTION :
  --   Application of the dynamic lifting algorithm to a tuple of supports,
  --   by the lift-and-prune algorithm.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   n          length of the vectors in the supports;
  --   mix        type of mixture;
  --   supports   lists of integer vector of length n;
  --   order      if true, then the points are already ordered and will
  --              be processed like they occur in the list,
  --              if false, then a random order will be chosen;
  --   inter      if true, then the list may contain interior points,
  --              i.e. points x in conv(l\{x}),
  --              if false, no interior points occur in the list;
  --   conmv      if true, then online checks on zero contributions will
  --              be made, otherwise this will not happen;
  --   maxli      maximum value of the lifting function,
  --              if = 0, then no flattening will be applied,
  --              if > 0, then no points will be given a lifting value
  --              larger than maxli, eventually flattening will be applied;
  --   mixsub     regular mixed subdivision of the lifted points;
  --   fs         the face structures of the lower hulls of the lifted
  --              points, necessary to re-start the algorithm efficiently;
  --   nbsucc     number of successful face-face combinations;
  --   nbfail     number of unsuccessful face-face combinations.
  --   q          a polynomial system with randomly chosen coefficients,
  --              the lists of exponent vectors of q equals the given supports.

  -- ON RETURN :
  --   mixsub     regular mixed subdivision of the lifted points;
  --   fs         face structures of the updated lower hulls;
  --   nbsucc     updated number of successful face-face combinations;
  --   nbfail     updated number of unsuccessful face-face combinations;
  --   qsols      the solutions of q.

end Dynamic_Polyhedral_Continuation;
