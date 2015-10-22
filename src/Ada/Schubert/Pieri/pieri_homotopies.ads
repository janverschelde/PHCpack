with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;
with Localization_Posets;                use Localization_Posets;

package Pieri_Homotopies is

-- DESCRIPTION :
--   This package provides the homotopy constructors for the poset-oriented
--   Pieri homotopy algorithm for four cases of increasing complexity :
--     1) hypersurface intersection conditions
--     2) general co-dimension intersections
--     3) q-curves satisfying interpolation-intersection conditions
--     4) q-curves the meet general planes of varying dimensions
--        at specified interpolation points.
--   The prefixes One_ and Two_ refer to the cases of respectively one
--   and two moving equations in the homotopy.

  function Moving_Parameter ( n,xk,tk : integer32;
                              start,target : Complex_Number ) return Poly;

  -- DESCRIPTION :
  --   Returns the equation c*(x-start)*(1-t) + (x-target)*t = 0 for
  --   the motion of x from start to target as t goes from 0 to 1.
  --   This version uses a constant c, randomly generated from within.

  -- ON ENTRY :
  --   n         total number of variables, continuation parameter t included;
  --   xk        index of the moving variable x;
  --   tk        index of the continuation parameter t;
  --   start     starting value for x;
  --   target    target value for x.

  function One_Hypersurface_Pieri_Homotopy
                  ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the Pieri homotopy for the hypersurface case,
  --   when the type of the node is either top or bottom,
  --   which means that only one intersection condition is folded in.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   nd           node in the localization poset, must be top or bottom;
  --   expbp        general format of the intersection condition;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions.

  function Two_Hypersurface_Pieri_Homotopy
                  ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the Pieri homotopy for the hypersurface case,
  --   when the type of the node is mixed,
  --   which means that two intersection conditions are folded in.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   nd           node in the localization poset, must be mixed;
  --   expbp        general format of the intersection condition;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions.

  function One_General_Pieri_Homotopy
                  ( n,ind : integer32; nd : Node; bs : Bracket_System;
                    start,target : Standard_Complex_Matrices.Matrix;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat ) return Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the Pieri homotopy to satisfy one general linear subspace
  --   intersection, when the type of the node is either top or bottom.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   ind          indicates the plane planes(ind) that is folded in;
  --   nd           node in the localization poset, must be top or bottom;
  --   bs           collects the structure of the equations;
  --   start        specialized plane that is met at the start of the homotopy;
  --   target       plane that has to be met at the end of the homotopy;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions.

  function Two_General_Pieri_Homotopy
                  ( n,ind : integer32; nd : Node;
                    top_bs,bot_bs : Bracket_System;
                    top_start,top_target,bot_start,bot_target
                      : Standard_Complex_Matrices.Matrix;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat ) return Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the Pieri homotopy to satisfy two general linear subspace
  --   intersections, when the type of the node is mixed.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   ind          indicates the plane planes(ind) that is folded in;
  --   nd           node in the localization poset, must be mixed;
  --   top_bs       the structure of the equations for top pivots;
  --   bot_bs       the structure of the equations for bottom pivots;
  --   top_start    special plane met at the start for top pivots;
  --   top_target   plane to be met at the end for top pivots;
  --   bot_start    special plane met at the start for bottom pivots;
  --   bot_target   plane to be met at the end for bottom pivots;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions.

  function One_Quantum_Pieri_Homotopy
                  ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat; s : Vector ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the Pieri homotopy to compute q-curves for one interpolation-
  --   intersection condition, when the node is either top or bottom.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   nd           node in the localization poset, must be top or bottom;
  --   expbp        general format of the intersection condition;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions;
  --   s            interpolation points where the planes are sampled.

  function Two_Quantum_Pieri_Homotopy
                  ( n : integer32; nd : Node; expbp : Bracket_Polynomial;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat; s : Vector ) return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the Pieri homotopy to compute q-curves for one interpolation-
  --   intersection condition, when the node is either top or bottom.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   nd           node in the localization poset, must be mixed;
  --   expbp        general format of the intersection condition;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions;
  --   s            interpolation points where the planes are sampled.

  function One_General_Quantum_Pieri_Homotopy
                  ( n,ind : integer32; nd : Node; s_mode : natural32;
                    bs : Bracket_System;
                    start,target : Standard_Complex_Matrices.Matrix;
                    xpm : Standard_Complex_Poly_Matrices.Matrix;
                    planes : VecMat; s : Vector ) return Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the quantum Pieri homotopy to satisfy one general linear
  --   subspace intersection, when the type is either top or bottom.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   ind          indicates the plane planes(ind) that is folded in;
  --   nd           node in the localization poset, must be top or bottom;
  --   s_mode       = 0 : s goes from 0 to 1,
  --                = 1 : s remains constant at 1,
  --                = 2 : s goes from 1 to target value s(ind);
  --   bs           collects the structure of the equations;
  --   start        specialized plane that is met at the start of the homotopy;
  --   target       plane that has to be met at the end of the homotopy;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions;
  --   s            interpolation points where the planes are sampled.

 -- function Two_General_Quantum_Pieri_Homotopy
 --               ( n,ind : natural; nd : Node; top_bs,bot_bs : Bracket_System;
 --                 top_start,top_target,bot_start,bot_target
 --                   : Standard_Complex_Matrices.Matrix;
 --                 xpm : Standard_Complex_Poly_Matrices.Matrix;
 --                 planes : VecMat; s : Vector ) return Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Returns the quantum Pieri homotopy to satisfy two general linear
  --   subspace intersections, when the type of the node is mixed.

  -- ON ENTRY :
  --   n            dimension of the working space, equals m+p;
  --   ind          indicates the plane planes(ind) that is folded in;
  --   nd           node in the localization poset, must be mixed;
  --   top_bs       the structure of the equations for top pivots;
  --   bot_bs       the structure of the equations for bottom pivots;
  --   top_start    special plane met at the start for top pivots;
  --   top_target   plane to be met at the end for top pivots;
  --   bot_start    special plane met at the start for bottom pivots;
  --   bot_target   plane to be met at the end for bottom pivots;
  --   xpm          localization pattern corresponding to the pivots in nd;
  --   planes       the planes that form the intersection conditions;
  --   s            interpolation points where the planes are sampled.

end Pieri_Homotopies;
