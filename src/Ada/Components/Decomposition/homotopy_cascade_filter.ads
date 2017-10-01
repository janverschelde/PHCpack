with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Irreducible_Decomp;
with Multprec_Irreducible_Decomp;

package Homotopy_Cascade_Filter is

-- DESCRIPTION :
--   This package provides routines to perform the witness generate and
--   classify, using a cascade of homotopies to filter solution lists.

-- INITIALIZE, UPDATE, AND MAIN LOOP :

  procedure Standard_Initialize
              ( soco : in out Standard_Irreducible_Decomp.Solution_Components;
                n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Initializes the solution components, placing the original system p
  --   at the 0-th component.

  procedure Multprec_Initialize
              ( soco : in out Multprec_Irreducible_Decomp.Solution_Components;
                n : in integer32;
                p : in Multprec_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   Initializes the solution components, placing the original system p
  --   at the 0-th component.

  procedure Standard_Update_Hypersurfaces
              ( file : in file_type;
                soco : in out Standard_Irreducible_Decomp.Solution_Components;
                n,top,k,itp : in natural32; skewproj : in boolean;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List; tol_sing : in double_float;
                timings : out Duration; fp,fp_last : in out List );

  -- DESCRIPTION :
  --   Updates the hypersurfaces and filter statistics.

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   soco       hypersurface equations for solution components;
  --   n          original dimension;
  --   top        top level;
  --   k          level, number of slices and added variables;
  --   itp        interpolation type :
  --               = 1 : massive interpolation with full grid of points,
  --               = 2 : incremental interpolation, one point after the other,
  --               = 3 : subspace restriction and projection from point;
  --   skewproj   if true, skew line projections will be used;
  --   embsys     embedding of the original polynomial system;
  --   sols       solutions to the embedded system;
  --   tol_sing   solutions with rcond lower than tol_sing are singular;
  --   fp         filter statistics;
  --   fp_last    pointer to the last element in the list fp.

  -- ON RETURN :
  --   soco       updated structure of hypersurface equations;
  --   timings    user cpu time for sampling and interpolating.
  --   fp,fp_last are updated filter statistics.

  procedure Multprec_Update_Hypersurfaces
              ( file : in file_type;
                soco : in out Multprec_Irreducible_Decomp.Solution_Components;
                n,top,k,size,itp : in natural32; skewproj : in boolean;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List; tol_sing : in double_float;
                timings : out Duration; fp,fp_last : in out List );

  -- DESCRIPTION :
  --   This is the multi-precision analogue for updating the hypersurface
  --   equations.  All parameters have the same meaning as the standard one.
  --   As additional parameter we now have "size" to denote the size of
  --   the numbers and the original system "orgsys" with multi-precision
  --   coefficients.

  procedure Standard_Cascade_Loop
              ( file : in file_type; n,k,itp : in natural32;
                skewproj : in boolean;
                embp : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                sols : in out Solution_List;
                soco : in out Standard_Irreducible_Decomp.Solution_Components;
                tol,tol_sing,tol_eval : in double_float;
                npa,ns0,ns1,div : in out Standard_Natural_Vectors.Vector;
                fp,fp_last : in out List;
                gentimes,clatimes : in out Array_of_Duration );

  -- DESCRIPTION :
  --   Performs a sequence of homotopies, starting at the top level k,
  --   down to level 0, using embeddings of the original polynomial system.
 
  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   n          dimension of the original polynomial system;
  --   k          top level, number of added slices and slack variables;
  --   itp        interpolation type :
  --               = 1 : massive interpolation with full grid of points,
  --               = 2 : incremental interpolation, one point after the other,
  --               = 3 : subspace restriction and projection from point;
  --   skewproj   if true, skew line projections will be used;
  --   embp       embedding sequence of the original polynomial system;
  --   sols       solutions at top level k;
  --   soco       initialized solution components;
  --   tol        tolerance for first filter of solution;
  --   tol_sing   tolerance for singularity test;
  --   tol_eval   tolerance for residual in component test;
  --   npa        #paths traced at every level, initialized at k;
  --   ns0        #paths ended with slack variable equal to zero, at level k;
  --   ns1        #paths ended with slack variable nonzero, at level k;
  --   div        #paths diverged to infinity, initialized at level k;
  --   fp         statistics about filtered points;
  --   fp_last    pointer to the last element in the list fp;
  --   gentimes   timings for the generation of generic points, initialized;
  --   clatimes   timings for the classfication of points.

  -- ON RETURN :
  --   sols       last list of computed solutions, sols is used as work space;
  --   soco       hypersurface equations for solution components;
  --   npa        updated #paths traced at each level;
  --   ns0        updated #paths ended at zero slack variables;
  --   ns1        updated #paths ended at nonzero slack variables;
  --   div        updated #paths diverged to infinity;
  --   gentimes   updated timings for computation of generic points;
  --   clatimes   updated timings for classification of points.

  procedure Multprec_Cascade_Loop
              ( file : in file_type; n,k,size,itp : in natural32;
                skewproj : in boolean;
                embp : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in out Solution_List;
                soco : in out Multprec_Irreducible_Decomp.Solution_Components;
                tol,tol_sing,tol_eval : in double_float;
                npa,ns0,ns1,div : in out Standard_Natural_Vectors.Vector;
                fp,fp_last : in out List;
                gentimes,clatimes : in out Array_of_Duration );

  -- DESCRIPTION :
  --   Performs a sequence of homotopies, starting at the top level k,
  --   down to level 0, using embeddings of the original polynomial system.
 
  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   n          dimension of the original polynomial system;
  --   k          top level, number of added slices and slack variables;
  --   size       size of the multi-precision numbers;
  --   itp        interpolation type :
  --               = 1 : massive interpolation with full grid of points,
  --               = 2 : incremental interpolation, one point after the other,
  --               = 3 : subspace restriction and projection from point;
  --   skewproj   if true, skew line projections will be used;
  --   embp       embedding sequence of the original polynomial system;
  --   orgsys     original system with multi-precision coefficients;
  --   sols       solutions at top level k;
  --   soco       initialized solution components;
  --   tol        tolerance for first filter of solution;
  --   tol_sing   tolerance for singularity test;
  --   tol_eval   tolerance for residual in component test;
  --   npa        #paths traced at every level, initialized at k;
  --   ns0        #paths ended with slack variable equal to zero, at level k;
  --   ns1        #paths ended with slack variable nonzero, at level k;
  --   div        #paths diverged to infinity, initialized at level k;
  --   fp         statistics about filtered points;
  --   fp_last    pointer to the last element in the list fp;
  --   gentimes   timings for the generation of generic points, initialized;
  --   clatimes   timings for the classfication of points.

  -- ON RETURN :
  --   sols       last list of computed solutions, sols is used as work space;
  --   soco       hypersurface equations for solution components;
  --   npa        updated #paths traced at each level;
  --   ns0        updated #paths ended at zero slack variables;
  --   ns1        updated #paths ended at nonzero slack variables;
  --   div        updated #paths diverged to infinity;
  --   gentimes   updated timings for computation of generic points;
  --   clatimes   updated timings for classification of points.

end Homotopy_Cascade_Filter;
