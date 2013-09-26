with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems; 
with Standard_Complex_Solutions; 
with Multprec_Complex_VecVecs;           use Multprec_Complex_VecVecs;
with Multprec_Complex_Poly_Systems;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Multprec_Irreducible_Decomp is

-- DESCRIPTION :
--   To represent a component of solutions we use three pieces of data:
--     1) list of generic points;
--     2) random hyperplanes to define a projection operator;
--     3) the hypersurface equations.
--   The generic points are solutions of the original polynomial system
--   sliced with as many random hyperplanes as the dimension of the
--   solution component.  The number of generic points is the degree
--   of the component.  We want to use only one equation to describe the
--   component.  Therefore we project a d-dimensional component to a space
--   of dimension d+1, using the random projection operator.  The projection
--   of a point is defined by the first d+1 entries in the evaluation of
--   the point in the random hyperplane.  The hypersurface equation is
--   obtained by interpolation in sample points.
--   Since the projection operator must be the same for all components
--   of the same dimension, we collect all d-dimensional components in
--   one record structure: Hypersurfaces.
--   The solution components are then stored as an array of hypersurfaces,
--   where the d-th entry holds the component of dimension d.
--   For d=0, the equations are the original polynomial equations and
--   the list pts contains the isolated solutions.

--   We represent the interpolating equations in multi-precision.
--   The generic points are originally given in standard precision,
--   and then later refined up to the required precision.
--   The random hyperplanes are represented in multi-precision format,
--   although they are only generated as double float vectors.

  type Hypersurfaces ( n,d : integer32 ) is record
    -- n is the dimension of the larger, working space
    -- d is the dimension of the solution component
    pts : Standard_Complex_Solutions.Solution_List;     -- generic points
    hyp : VecVec(1..n);                                 -- random hyperplanes 
    equ : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys; 
                                                    -- hypersurface equations
  end record;

  type Link_to_Hypersurfaces is access Hypersurfaces;

  type Solution_Components is
    array ( integer32 range <> ) of Link_to_Hypersurfaces;

-- CREATOR :

  procedure Create_Hypersurfaces
              ( file : in file_type; n,d,size,itp : in natural32;
                skewproj : in boolean;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                genpts : in Standard_Complex_Solutions.Solution_List;
                surf : out Link_to_Hypersurfaces; timing : out Duration );

  -- DESCRIPTION :
  --   With the generic points genpts as solutions of the embedded system,
  --   the corresponding hypersurfaces are constructed.
  --   Note that there is no sharing, the generic points are duplicated.

  -- ON ENTRY :
  --   file     to write intermediate results;
  --   n        working dimension of the system and solutions;
  --   d        level, n-d is the dimension of the original system;
  --   size     size of the floating-point numbers in interpolator;
  --   itp       interpolation type :
  --              = 1 : massive interpolation with full grid of points,
  --              = 2 : incremental interpolation, one point after the other,
  --              = 3 : subspace restriction and projection from point;
  --   embsys   the embedding of the original system;
  --   orgsys   original system with multi-precision coefficients;
  --   genpts   generic points are solutions to the embedded system;

  -- ON RETURN :
  --   surf     created hypersurface equations;
  --   timing   user cpu time needed for sampling and interpolating.

-- SELECTORS :

  procedure Component_Test
              ( file : in file_type; level : in natural32;
                surf : in Hypersurfaces; tol : in double_float;
                sols : in Standard_Complex_Solutions.Solution_List;
                solslvl : in natural32; num : out natural32;
                filpts : in out List;
                rem_first,rem_last
                  : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Evaluates the solutions sols at level solslvl after projection on
  --   the hyperplane of the hypersurfaces in the hypersurface equations
  --   to determine whether they lie on the solution component.
  --   Results are written to a log file.
  --   The number "num" on return counts the number of solutions that
  --   satisfy the interpolating polynomials of the Hypersurfaces.
  --   The filtered points statistical data is updated.
  --   The remainder of the solutions that fail the component test are
  --   stored in the list with head rem_first and tail rem_last.

  procedure Component_Test
              ( file : in file_type; soco : in Solution_Components;
                tol : in double_float;
                sols : in Standard_Complex_Solutions.Solution_List;
                level : in natural32; num : out natural32;
                filpts : in out List;
                rem_sols : out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs the component test on the solution list sols, along all
  --   components in soco.  Results are written to a log file.
  --   The solutions that do not lie on the components are in rem_sols.

-- DESTRUCTORS :

  procedure Clear ( h : in out Hypersurfaces );
  procedure Clear ( h : in out Link_to_Hypersurfaces );
  procedure Clear ( s : in out Solution_Components );

  -- DESCRIPTION :
  --   Deallocation of the memory occupied by the data structures.

end Multprec_Irreducible_Decomp;
