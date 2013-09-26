with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;     use Standard_Complex_Polynomials;
with Partitions_of_Sets_of_Unknowns;   use Partitions_of_Sets_of_Unknowns;

package Drivers_to_Intrinsic_Solvers is

-- DESCRIPTION :
--   This package provides drivers to the intrinsic implementation
--   of diagonal homotopies.

-- DRIVERS TO SOLVERS OF TWO POLYNOMIALS :

  procedure Total_Degree_Hypersurface_Intersection
               ( n : in natural32; p,q : in Poly );
  procedure Total_Degree_Hypersurface_Intersection
               ( file : in file_type; n : in natural32; p,q : in Poly );

  -- DESCRIPTION :
  --   Computes witness points on the hypersurfaces defined by p and q,
  --   and then calls a total degree diagonal homotopy to compute witness
  --   points on the intersection of the two hyperfaces.

  procedure Multihomogeneous_Hypersurface_Intersection
               ( n : in natural32; p,q : in Poly; z : in Partition );
  procedure Multihomogeneous_Hypersurface_Intersection
               ( file : in file_type;
                 n : in natural32; p,q : in Poly; z : in Partition );

  -- DESCRIPTION :
  --   Computes witness points on hypersurfaces with a multi-homogeneous
  --   structure, given by the partition z.

-- MAIN INTERACTIVE DRIVERS :

  procedure Intersection_of_Two_Hypersurfaces;

  -- DESCRIPTION :
  --   Driver routine to intersect two hypersurfaces defined by
  --   polynomials in the same number of variables and with standard
  --   complex coefficients.

  procedure Driver_to_Hypersurface_Solvers;

  -- DESCRIPTION :
  --   Interactive driver to solve polynomial systems incrementally.

end Drivers_to_Intrinsic_Solvers;
