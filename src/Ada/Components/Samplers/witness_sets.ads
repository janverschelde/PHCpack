with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Witness_Sets is

-- DESCRIPTION :
--   To compute generic points on components of solutions, we add
--   a random hyperplane to the system and apply an embedding with
--   an additional variable to square the overdetermined systems.
--   This gives rise to the "witness set" data structure which consists
--   of an embedded system and a set of its isolated solutions.
--   The number of hyperplanes added to the system equals the dimension
--   and the number of solutions equals the degree of the algebraic set
--   represented by the witness set.  The operations in this package
--   provide tools to define witness sets.

-- GENERATE RANDOM SLICES :

  function Random_Hyperplanes ( k,n : natural32 )
                              return Standard_Complex_VecVecs.VecVec;
  function Random_Hyperplanes ( k,n : natural32 )
                              return DoblDobl_Complex_VecVecs.VecVec;
  function Random_Hyperplanes ( k,n : natural32 )
                              return QuadDobl_Complex_VecVecs.VecVec;
  function Random_Hyperplanes ( k,n,size : natural32 )
                              return Multprec_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns k random hyperplanes, represented as vectors of range 0..n,
  --   either with coefficients as standard complex numbers,
  --   or with multi-precision complex numbers of the given size.

-- MANIPULATION OF POLYNOMIALS :

  function Add_Dummy ( n,k,i : natural32 )
                     return Standard_Complex_Polynomials.Poly;
  function Add_Dummy ( n,k,i : natural32 )
                     return DoblDobl_Complex_Polynomials.Poly;
  function Add_Dummy ( n,k,i : natural32 )
                     return QuadDobl_Complex_Polynomials.Poly;
  function Add_Dummy ( n,k,i : natural32 )
                     return Standard_Complex_Laurentials.Poly;
  function Add_Dummy ( n,k,i : natural32 )
                     return DoblDobl_Complex_Laurentials.Poly;
  function Add_Dummy ( n,k,i : natural32 )
                     return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns one monomial: "zzi", where i is the number of slack
  --   variable ranging between 1 and k.  The number of variables
  --   in the polynomial on return is n+k.

  function Add_Dummies 
             ( n,k : natural32 ) return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Adds k linear terms in the k extra slack variables to an empty
  --   polynomial in n variables,
  --   so the polynomial on return has n+k variables.

  function Add_Embedding ( p : Standard_Complex_Polynomials.Poly;
                           k : natural32 )
                         return Standard_Complex_Polynomials.Poly;
  function Add_Embedding ( p : Standard_Complex_Laurentials.Poly;
                           k : natural32 )
                         return Standard_Complex_Laurentials.Poly;
  function Add_Embedding ( p : DoblDobl_Complex_Polynomials.Poly;
                           k : natural32 )
                         return DoblDobl_Complex_Polynomials.Poly;
  function Add_Embedding ( p : DoblDobl_Complex_Laurentials.Poly;
                           k : natural32 )
                         return DoblDobl_Complex_Laurentials.Poly;
  function Add_Embedding ( p : QuadDobl_Complex_Polynomials.Poly;
                           k : natural32 )
                         return QuadDobl_Complex_Polynomials.Poly;
  function Add_Embedding ( p : QuadDobl_Complex_Laurentials.Poly;
                           k : natural32 )
                         return QuadDobl_Complex_Laurentials.Poly;

  function Add_Embedding ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                           k : natural32 )
                         return Standard_Complex_Poly_Systems.Poly_Sys;
  function Add_Embedding ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                           k : natural32 )
                         return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Does Add_Variables(p,k) and adds to the result k linear terms
  --   with random coefficients in the k extra variables.
  --   On a system p, the system on return has the same range.

-- OPERATIONS ON POLYNOMIAL SYSTEMS :

  function Make_Square ( f : Standard_Complex_Poly_Systems.Poly_Sys;
                         k : natural32 )
                       return Standard_Complex_Poly_Systems.Poly_Sys;
  function Make_Square ( f : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                         k : natural32 )
                       return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Make_Square ( f : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                         k : natural32 )
                       return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   If the number of equations in f is more than the co-dimension k,
  --   the extra equations in f will be multiplied with a random
  --   constant and added to the first k equation in f.
  --   The system on return is always a copy, independent of f.

  function Add_Slice ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                       hyp : Standard_Complex_Vectors.Vector )
                     return Standard_Complex_Poly_Systems.Poly_Sys;
  function Add_Slice ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                       hyp : DoblDobl_Complex_Vectors.Vector )
                     return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Add_Slice ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                       hyp : QuadDobl_Complex_Vectors.Vector )
                     return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Adds the hyperplane to the polynomial system p.  The coefficients of
  --   the hyperplane are indexed from 0 up to n = p'last, to represent
  --     hyp(0) + hyp(1)*x(1) + hyp(2)*x(2) + .. + hyp(n)*x(n) = 0, 
  --   The last equation of the system on return is the representation
  --   of the hyperplane as a polynomial.
  --   The equations of p are shared with those of the system on return.

  function Remove_Slice ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                        return Standard_Complex_Poly_Systems.Poly_Sys;
  function Remove_Slice ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                        return Standard_Complex_Laur_Systems.Laur_Sys;
  function Remove_Slice ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                        return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Remove_Slice ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                        return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Remove_Slice ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                        return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Remove_Slice ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                        return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Replaces the last equation of p by x_n = 0, where n = p'last. 

  function Eliminate_Slice ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                             k,i : natural32 )
                           return Standard_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the polynomial system where the i-th unknown has been
  --   eliminated, using the k-th slice in the equations of p.

  -- REQUIRED : Degree(p(k),i) = 1.

  function Embed ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                 return Standard_Complex_Poly_Systems.Poly_Sys;
  function Embed ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                 return Standard_Complex_Laur_Systems.Laur_Sys;
  function Embed ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                 return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Embed ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                 return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Embed ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                 return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Embed ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                 return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Augments the number of variables in the polynomials with one,
  --   in standard double, double double, and quad double precision,
  --   every monomial is multiplied with z^0.

  function Embed ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                   gamma : Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Poly_Systems.Poly_Sys;
  function Embed ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                   gamma : DoblDobl_Complex_Vectors.Vector )
                 return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Embed ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                   gamma : QuadDobl_Complex_Vectors.Vector )
                 return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   To every polynomial in p the term gamma(i)*z is added,
  --   in standard double, double double, and quad double precision,
  --   where z is a new variable.

  procedure Store_Random_Hyperplanes
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                n,k : in natural32 );
  procedure Store_Random_Hyperplanes
              ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                n,k : in natural32 );
  procedure Store_Random_Hyperplanes
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                n,k : in natural32 );
  procedure Store_Random_Hyperplanes
              ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                n,k : in natural32 );
  procedure Store_Random_Hyperplanes
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                n,k : in natural32 );
  procedure Store_Random_Hyperplanes
              ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                n,k : in natural32 );

  -- DESCRIPTION :
  --   In the last k entries of p, random hyperplanes are stored,
  --   where n is the last index in p of an input polynomial.

  function Slice_and_Embed
              ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys;
  function Slice_and_Embed
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return Standard_Complex_Laur_Systems.Laur_Sys;
  function Slice_and_Embed
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Slice_and_Embed
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Slice_and_Embed
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Slice_and_Embed
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Adds k slices and k slack variables to the system p,
  --   in standard double, double double, and quad double precision.

  -- REQUIRED : p is a square system.

  function Embed_with_Dummies
              ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys;
  function Embed_with_Dummies
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return Standard_Complex_Laur_Systems.Laur_Sys;
  function Embed_with_Dummies
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Embed_with_Dummies
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Embed_with_Dummies
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                k : natural32 )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Embed_with_Dummies
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                k : natural32 )
              return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Adds k dummy slack variables with k extra equations setting those
  --   dummy variables to zero.  The original last k hyperplanes of p
  --   (see the requirement) are replaced with the equations that set
  --   the dummy variables to zero.

  -- REQUIRED :
  --   p is a square system and contains at the end at least k hyperplanes.

  function Slice_and_Embed 
              ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : natural32 ) 
              return Standard_Complex_Poly_Systems.Array_of_Poly_Sys;

  -- DESCRIPTION :
  --   Returns a sequence of embedded systems of range 0..k.
  --   The 0-th entry is the original system p, the i-th entry is obtained
  --   by adding i random hyperplanes and i additional variables to p.
  --   If k > 1, then the slices are orthogonal with respect to each other.

  function Slice_and_Embed
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; k,l : natural32 )
             return Standard_Complex_Poly_Systems.Array_of_Poly_Sys;

  -- DESCRIPTION :
  --   In constructing the sequence of embedded systems, like above,
  --   the last l equations of p are copied and used as slices.
  --   This is useful in the treatment of underdetermined systems.

  -- REQUIRED : l <= k.

  function Slices ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                    k : natural32 ) return Standard_Complex_VecVecs.VecVec;
  function Slices ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                    k : natural32 ) return Standard_Complex_VecVecs.VecVec;
  function Slices ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    k : natural32 ) return DoblDobl_Complex_VecVecs.VecVec;
  function Slices ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    k : natural32 ) return DoblDobl_Complex_VecVecs.VecVec;
  function Slices ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    k : natural32 ) return QuadDobl_Complex_VecVecs.VecVec;
  function Slices ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    k : natural32 ) return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the vector representation of the k slices added to p.
  --   Each vector in the vector of vectors of return contains the
  --   coefficients cff:  cff(0) + cff(1)*x1 + .. + cff(n)*xn.

  function Square ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Square ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys;
  function Square ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Square ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Square ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return QUadDobl_Complex_Poly_Systems.Poly_Sys;
  function Square ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return QUadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   The system on return is square, has as many equations as unknowns.
  --   If p is square, then Square(p) = p, otherwise slices or additional
  --   unknowns are added depending on whether p is has too few equations
  --   or too few unknowns.

  function Remove_Embedding ( p : Standard_Complex_Polynomials.Poly;
                              dim : natural32 )
                            return Standard_Complex_Polynomials.Poly;
  function Remove_Embedding ( p : Standard_Complex_Laurentials.Poly;
                              dim : natural32 )
                            return Standard_Complex_Laurentials.Poly;
  function Remove_Embedding ( p : DoblDobl_Complex_Polynomials.Poly;
                              dim : natural32 )
                            return DoblDobl_Complex_Polynomials.Poly;
  function Remove_Embedding ( p : DoblDobl_Complex_Laurentials.Poly;
                              dim : natural32 )
                            return DoblDobl_Complex_Laurentials.Poly;
  function Remove_Embedding ( p : QuadDobl_Complex_Polynomials.Poly;
                              dim : natural32 )
                            return QuadDobl_Complex_Polynomials.Poly;
  function Remove_Embedding ( p : QuadDobl_Complex_Laurentials.Poly;
                              dim : natural32 )
                            return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Removes all dim added slack  added in the embedding.

  function Remove_Embedding1 ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                               dim : natural32 )
                             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Remove_Embedding1 ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                               dim : natural32 )
                             return Standard_Complex_Laur_Systems.Laur_Sys;
  function Remove_Embedding1 ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                               dim : natural32 )
                             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Remove_Embedding1 ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                               dim : natural32 )
                             return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Remove_Embedding1 ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                               dim : natural32 )
                             return QuadDobl_Complex_Poly_Systems.Poly_Sys;
  function Remove_Embedding1 ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                               dim : natural32 )
                             return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Removes all dim added variables and the terms added to p in
  --   the embedding.  On return is a system of range 1..p'last-dim,
  --   but more equations at the end may have become zero.

  function Number_of_Zero_Equations 
             ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return natural32;
  function Number_of_Zero_Equations 
             ( p : Standard_Complex_Laur_Systems.Laur_Sys ) return natural32;
  function Number_of_Zero_Equations 
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return natural32;
  function Number_of_Zero_Equations 
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys ) return natural32;
  function Number_of_Zero_Equations 
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return natural32;
  function Number_of_Zero_Equations 
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of zero equations that may occur at the 
  --   end after removing an embedding.
  
  function Complete ( n,k : natural32;
                      p : Standard_Complex_Poly_Systems.Poly_Sys )
                    return Standard_Complex_Poly_Systems.Poly_Sys;
  function Complete ( n,k : natural32;
                      p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                    return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Complete ( n,k : natural32;
                      p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                    return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   When p is not a complete intersection, random multiples of the
  --   extra polynomials are added to have n-k equations on return.

  -- ON ENTRY :
  --   n         number of variables;
  --   k         dimension of the solution component;
  --   p         polynomial system.

  -- ON RETURN :
  --   n-k equations which also define the solution component.
  --   In case p'length is already n-k, a copy of p is returned.
  
  function Complete ( n,k : natural32;
                      p : Standard_Complex_Laur_Systems.Laur_Sys )
                    return Standard_Complex_Laur_Systems.Laur_Sys;
  function Complete ( n,k : natural32;
                      p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                    return DoblDobl_Complex_Laur_Systems.Laur_Sys;
  function Complete ( n,k : natural32;
                      p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                    return QuadDobl_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   When p is not a complete intersection, random multiples of the
  --   extra polynomials are added to have n-k equations on return.

  -- ON ENTRY :
  --   n         number of variables;
  --   k         dimension of the solution component;
  --   p         polynomial system.

  -- ON RETURN :
  --   n-k equations which also define the solution component.
  --   In case p'length is already n-k, a copy of p is returned.

-- OPERATIONS ON SOLUTION LISTS :

  function Add_Component
             ( s : Standard_Complex_Solutions.Solution;
               c : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Solutions.Solution;
 
  -- DESCRIPTION :
  --   Add c as last component to the solution vector.

  function Add_Component
             ( sols : Standard_Complex_Solutions.Solution_List;
               c : Standard_Complex_Numbers.Complex_Number )
             return Standard_Complex_Solutions.Solution_List;
 
  -- DESCRIPTION :
  --   Add c as last component to the solution vectors in the list.

  function Add_Embedding
             ( s : Standard_Complex_Solutions.Solution; k : natural32 )
             return Standard_Complex_Solutions.Solution;
  function Add_Embedding
             ( s : DoblDobl_Complex_Solutions.Solution; k : natural32 )
             return DoblDobl_Complex_Solutions.Solution;
  function Add_Embedding
             ( s : QuadDobl_Complex_Solutions.Solution; k : natural32 )
             return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Adds k zeros to every solution,
  --   in standard double, double double, or quad double precision.

  function Add_Embedding
             ( sols : Standard_Complex_Solutions.Solution_List;
               k : natural32 )
             return Standard_Complex_Solutions.Solution_List;
  function Add_Embedding
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List;
  function Add_Embedding
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Adds k zeros to every solution,
  --   in standard double, double double, or quad double precision.

  function Remove_Component
             ( s : Standard_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution;
  function Remove_Component
             ( s : DoblDobl_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution;
  function Remove_Component
             ( s : QuadDobl_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution;
 
  -- DESCRIPTION :
  --   Removes the last component from the solution vector,
  --   in standard double, double double, or quad double precision.

  procedure Remove_Component
              ( sols : in out Standard_Complex_Solutions.Solution_List );
 
  -- DESCRIPTION :
  --   Removes the last component from the solution vectors in the list.

  function Remove_Component
              ( sols : Standard_Complex_Solutions.Solution_List )
              return Standard_Complex_Solutions.Solution_List;
  function Remove_Component
              ( sols : DoblDobl_Complex_Solutions.Solution_List )
              return DoblDobl_Complex_Solutions.Solution_List;
  function Remove_Component
              ( sols : QuadDobl_Complex_Solutions.Solution_List )
              return QuadDobl_Complex_Solutions.Solution_List;
 
  -- DESCRIPTION :
  --   Removes the last component from the solution vectors in the list,
  --   in standard double, double double, or quad double precision.

  function Remove_Embedding
             ( s : Standard_Complex_Solutions.Solution; k : natural32 )
             return Standard_Complex_Solutions.Solution;
  function Remove_Embedding
             ( s : DoblDobl_Complex_Solutions.Solution; k : natural32 )
             return DoblDobl_Complex_Solutions.Solution;
  function Remove_Embedding
             ( s : QuadDobl_Complex_Solutions.Solution; k : natural32 )
             return QuadDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Removes the last k components from the solution,
  --   in standard double, double double, or quad double precision.

  function Remove_Embedding
             ( sols : Standard_Complex_Solutions.Solution_List;
               k : natural32 )
             return Standard_Complex_Solutions.Solution_List;
  function Remove_Embedding
             ( sols : DoblDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List;
  function Remove_Embedding
             ( sols : QuadDobl_Complex_Solutions.Solution_List;
               k : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Removes the last k components from each solution in the list,
  --   in standard double, double double, or quad double precision.

end Witness_Sets;
