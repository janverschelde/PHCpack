with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with TripDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with PentDobl_Complex_Solutions;
with OctoDobl_Complex_Solutions;
with DecaDobl_Complex_Solutions;
with HexaDobl_Complex_Solutions;
with Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

package Multi_Projective_Transformations is

-- DESCRIPTION :
--   A multi-projective space is defined by a partition of the set of unknowns.
--   After a multi-projective transformation, the degree of a term in a set
--   of the partition is the same for all terms in the multi-homogenous form
--   of the polynomial.  Transformations are supported for polynomials with
--   coefficients in double, double double, triple double, quad double,
--   penta double, octo double, deca double, and hexa double precision.

  function Multiset_Degrees
             ( p : in Standard_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in DoblDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in TripDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in QuadDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in PentDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in OctoDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in DecaDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;
  function Multiset_Degrees
             ( p : in HexaDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the degrees of the polynomial p in the m sets
  --   of the partition z, as a vector of range 1..m.

  function Make_Homogeneous
             ( t : Standard_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return Standard_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : DoblDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return DoblDobl_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : TripDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return TripDobl_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : QuadDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return QuadDobl_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : PentDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return PentDobl_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : OctoDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return OctoDobl_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : DecaDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return DecaDobl_Complex_Polynomials.Term;
  function Make_Homogeneous
             ( t : HexaDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return HexaDobl_Complex_Polynomials.Term;

  -- DESCRIPTION :
  --   Returns the term with m variables added,
  --   with degrees to match the given multiset degrees in d.

  function Make_Homogeneous
             ( p : Standard_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return Standard_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : DoblDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return DoblDobl_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : TripDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return TripDobl_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : QuadDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return QuadDobl_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : PentDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return PentDobl_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : OctoDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return OctoDobl_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : DecaDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return DecaDobl_Complex_Polynomials.Poly;
  function Make_Homogeneous
             ( p : HexaDobl_Complex_Polynomials.Poly; 
               m : natural32; z : Partition )
             return HexaDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the polynomial p with m variables added,
  --   with degrees to match the given multiset degrees in d.

  function Make_Homogeneous
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Homogeneous
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system p with m variables added,
  --   with degrees to match the multiset degrees of each polynomial.

  function Standard_Random_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term;
  function DoblDobl_Random_Linear_Term
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Term;
  function TripDobl_Random_Linear_Term
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Term;
  function QuadDobl_Random_Linear_Term
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Term;
  function PentDobl_Random_Linear_Term
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Term;
  function OctoDobl_Random_Linear_Term
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Term;
  function DecaDobl_Random_Linear_Term
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Term;
  function HexaDobl_Random_Linear_Term
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Term;

  -- DESCRIPTION :
  --   Returns a term in the i-th variable, with random coefficient,
  --   as a term in n variables.

  function Standard_Start_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term;
  function DoblDobl_Start_Linear_Term
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Term;
  function TripDobl_Start_Linear_Term
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Term;
  function QuadDobl_Start_Linear_Term
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Term;
  function PentDobl_Start_Linear_Term
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Term;
  function OctoDobl_Start_Linear_Term
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Term;
  function DecaDobl_Start_Linear_Term
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Term;
  function HexaDobl_Start_Linear_Term
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Term;

  -- DESCRIPTION :
  --   Returns the i-th variable as a term in n variables,
  --   with coefficient equal to one.

  function Standard_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return Standard_Complex_Polynomials.Poly;
  function DoblDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return DoblDobl_Complex_Polynomials.Poly;
  function TripDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return TripDobl_Complex_Polynomials.Poly;
  function QuadDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return QuadDobl_Complex_Polynomials.Poly;
  function PentDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return PentDobl_Complex_Polynomials.Poly;
  function OctoDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return OctoDobl_Complex_Polynomials.Poly;
  function DecaDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return DecaDobl_Complex_Polynomials.Poly;
  function HexaDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return HexaDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a linear polynomial in the variables in s,
  --   with nonzero constant coefficient as a polynomial in n variables.

  function Standard_Start_Linear_Polynomial
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Poly;
  function DoblDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Poly;
  function TripDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Poly;
  function QuadDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Poly;
  function PentDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Poly;
  function OctoDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Poly;
  function DecaDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Poly;
  function HexaDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the start polynomial Zi - 1,
  --   as a polynomial in n variables.

  function Standard_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function DoblDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function TripDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function QuadDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function PentDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;
  function OctoDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;
  function DecaDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;
  function HexaDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns m random linear polynomials in n+m variables in the sets
  --   of the partition z, with the constant added and the random term
  --   for the extra homogeneous m-th variable.

  function Standard_Start_Linear_Polynomials
             ( n,m : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function DoblDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function TripDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function QuadDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function PentDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;
  function OctoDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;
  function DecaDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;
  function HexaDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns m start polynomials in n+m variables Zi - 1,
  --   for i in range 1..m.

  function Add_Ones ( s : Standard_Complex_Solutions.Solution;
                      m : natural32 )
                    return Standard_Complex_Solutions.Solution;
  function Add_Ones ( s : DoblDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return DoblDobl_Complex_Solutions.Solution;
  function Add_Ones ( s : TripDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return TripDobl_Complex_Solutions.Solution;
  function Add_Ones ( s : QuadDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return QuadDobl_Complex_Solutions.Solution;
  function Add_Ones ( s : PentDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return PentDobl_Complex_Solutions.Solution;
  function Add_Ones ( s : OctoDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return OctoDobl_Complex_Solutions.Solution;
  function Add_Ones ( s : DecaDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return DecaDobl_Complex_Solutions.Solution;
  function Add_Ones ( s : HexaDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return HexaDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Returns a solution with the same coordinates as in s,
  --   but with m added coordinates, with values all equal to one.

  function Add_Ones ( sols : Standard_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return Standard_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : DoblDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return DoblDobl_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : TripDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return TripDobl_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : QuadDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return QuadDobl_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : PentDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return PentDobl_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : OctoDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return OctoDobl_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : DecaDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return DecaDobl_Complex_Solutions.Solution_List;
  function Add_Ones ( sols : HexaDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return HexaDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the solutions with the same coordinates as in sols,
  --   but with m coordinates, all equal to one, added to each solution.

  procedure Add_Ones ( sols : in out Standard_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out TripDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out PentDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out OctoDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out DecaDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );
  procedure Add_Ones ( sols : in out HexaDobl_Complex_Solutions.Solution_List;
                       m : in natural32 );

  -- DESCRIPTION :
  --   Replaces every solution in sols by a solution with the same
  --   coordinates, and with m ones added to each solution vector.

  function Make_Affine ( sol : Standard_Complex_Solutions.Solution;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return Standard_Complex_Solutions.Solution;
  function Make_Affine ( sol : DoblDobl_Complex_Solutions.Solution;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return DoblDobl_Complex_Solutions.Solution;
  function Make_Affine ( sol : QuadDobl_Complex_Solutions.Solution;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Divides every coordinate in sol.v by the value of the added
  --   homogeneous coordinate of each set, as defined by idz,
  --   the index representation of the partition of the variables.

  function Make_Affine ( sols : Standard_Complex_Solutions.Solution_List;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return Standard_Complex_Solutions.Solution_List;
  function Make_Affine ( sols : DoblDobl_Complex_Solutions.Solution_List;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return DoblDobl_Complex_Solutions.Solution_List;
  function Make_Affine ( sols : QuadDobl_Complex_Solutions.Solution_List;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Returns the affine version of the solutions in sols,
  --   dividing every coordinate in a solution by the value of the
  --   added homogeneous coordinate of each set, as defined by idz,
  --   the index representation of the partition of the variables.

  procedure Make_Affine
              ( sols : in out Standard_Complex_Solutions.Solution_List;
                m : in natural32;
                idz : in Standard_Natural_Vectors.Vector );
  procedure Make_Affine
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                m : in natural32;
                idz : in Standard_Natural_Vectors.Vector );
  procedure Make_Affine
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                m : in natural32;
                idz : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Converts the list of solutions sols to affine coordinates,
  --   dividing every coordinate in a solution by the value of the
  --   added homogeneous coordinate of each set, as defined by idz,
  --   the index representation of the partition of the variables.

  function Multi_Projective_Transformation
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return TripDobl_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return PentDobl_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys;
  function Multi_Projective_Transformation
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system p with m variables added,
  --   with degrees to match the multiset degrees of each polynomial.
  --   If start, then start linear polynomials are added,
  --   otherwise, random linear polynomial are added.

end Multi_Projective_Transformations;
