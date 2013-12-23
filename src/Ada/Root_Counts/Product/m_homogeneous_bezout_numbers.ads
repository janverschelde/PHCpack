with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package m_Homogeneous_Bezout_Numbers is

-- DESCRIPTION :
--   This package allows the computation of m-homogeneous Bezout numbers.
--   It provides various enumeration strategies for computing a minimal
--   m-homogeneous Bezout number.

  function Cardinalities ( z : Partition ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector which contains the cardinality of each set.

  function Total_Degree ( p : Poly_Sys ) return natural32;
  function Total_Degree ( p : Poly_Sys ) return natural64;
  function Total_Degree ( p : Poly_Sys ) return Natural_Number;

  -- DESCRIPTION :
  --   Returns the 1-homogeneous Bezout number of the system.

  function Bezout_Number ( p : Poly_Sys; z : Partition ) return natural64;
  function Bezout_Number ( p : Poly_Sys; z : Partition; max : natural64 )
                         return natural64;

  -- DESCRIPTION :
  --   Returns the m-homogeneous Bezout number w.r.t. the given partition.
  --   When max is given as parameter, the computation stops when the
  --   result becomes larger than or equal to max.

  function Bezout_Number ( p : Poly_Sys ) return natural64;
  function Bezout_Number ( max : natural64; p : Poly_Sys ) return natural64;
  function Bezout_Number ( p : Poly_Sys; min : natural64 ) return natural64;
  function Bezout_Number ( max : natural64; p : Poly_Sys; min : natural64 )
                         return natural64;
  -- DESCRIPTION :
  --   A minimal m-homogeneous Bezout number of the polynomial system 
  --   p is computed, by generating partitions of the sets of unknowns.

  -- ON ENTRY :
  --   p         a polynomial system;
  --   max       a maximum number of partitions to be evaluated,
  --             if max=0, then the total degree is returned;
  --   min       the procedure stops when the Bezout number becomes
  --             smaller than min.

  procedure Bezout_Number 
               ( p : in Poly_Sys; b : out natural64;
                 m : out natural32; z : in out Partition );
  procedure Bezout_Number 
               ( max : in natural64; p : in Poly_Sys; b : out natural64;
                 m : out natural32; z : in out Partition );
  procedure Bezout_Number
               ( p : in Poly_Sys; min : in natural64; b : out natural64;
                 m : out natural32; z : in out Partition );
  procedure Bezout_Number
               ( max : in natural64; p : in Poly_Sys; min : in natural64;
                 b : out natural64; m : out natural32; z : in out Partition );

  -- DESCRIPTION :
  --   A minimal m-homogeneous Bezout number of the polynomial system
  --   is computed.  The partition with the calculated Bezout number
  --   is returned.

  -- ON ENTRY :
  --   p         a polynomial system;
  --   max       a maximum number of partition to be evaluated,
  --             if max=0, then the total degree is returned;
  --   min       the procedure stops when the Bezout number becomes
  --             smaller than min.

  -- ON RETURN :
  --   b         a minimal m-homogeneous Bezout number;
  --   m         the number of sets in the partition z;
  --   z         the partition with the computed Bezout number b.

  procedure PB ( p : in Poly_Sys; b : out natural64;
                 m : out natural32; z : in out Partition );

  -- DESCRIPTION :
  --   This is a fast heuristic for computing `the' Bezout number of p.
  --   It is fast because it does not generate all partitions,
  --   it is a heuristic because it does not always give then
  --   minimal partition.

  -- ON ENTRY :
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         an m-homogeneous Bezout number;
  --   m         the number of sets in the partition z;
  --   z         the partition corresponding to b.

  procedure Patch ( p : in Poly_Sys; z : in out Partition;
                    nz : in natural32; bz : in out natural64 );

  -- DESCRIPTION :
  --   Because of this bug :
  --    3 
  --     x**2 + x*y + z;
  --     y**2 + x*y + z;
  --     x + y + z;
  --   The heuristic partitioner only returns {z} instead of {x y z}.

end m_Homogeneous_Bezout_Numbers;
