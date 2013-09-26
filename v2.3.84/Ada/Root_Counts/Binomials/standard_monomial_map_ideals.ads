with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Lists;       use Standard_Complex_Laur_Lists;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;

package Standard_Monomial_Map_Ideals is

-- DESCRIPTION :
--   A monomial map defines each variable as a monomial of parameters.
--   This parametric representation allows to set up the defining
--   equations of the solution set represented by the monomial map.

  function One_Variable_Equation ( n,i : integer32 ) return Poly;

  -- DESCRIPTION :
  --   Returns the equation defined by setting the i-th variable to zero,
  --   represented as one term of a polynomial in n variables.

  function Is_Constant ( map : Monomial_Map; i : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the i-th variable as defined by the map is constant,
  --   returns false otherwise.

  function One_Variable_Constant_Equation
             ( n,i : integer32; c : Complex_Number ) return Poly;

  -- DESCRIPTION :
  --   Returns the polynomial x(i) - c as a polynomial in n variables.

  function Variable_Difference_Equation
              ( n,i,j : integer32; ci,cj : Complex_Number ) return Poly;

  -- DESCRIPTION :
  --   Returns cj*x(i) - ci*x(j) as a polynomial in n variables.

  function Same_Exponents
              ( map : Monomial_Map; i,j : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x(i) and x(j) share the same exponents.

  function Variable_Differences ( map : Monomial_Map ) return Poly_List;

  -- DESCRIPTION :
  --   Returns the equations of scaled differences of variables,
  --   as defined by the monomial map.

  function Share_One_Same_Parameter
              ( map : Monomial_Map; i,j : integer32 ) return integer32;

  -- DESCRIPTION :
  --   If the monomials for x(i) and x(j) differ only in the power of one
  --   and the same parameter, then the index of that parameter is returned,
  --   otherwise 0 is returned.

  function Variable_Relation_Same_Sign
              ( n,i,j,a,b,e : integer32;
                c,d : Complex_Number ) return Poly;

  -- DESCRIPTION :
  --   Computes the relation between variables i and j
  --   with x(i) = c*t^a, x(j) = d*t^b, e = lcm(|a|,|b|)
  --   where a and b have the same sign.

  function Variable_Relation_Opposite_Sign
              ( n,i,j,a,b,e : integer32;
                c,d : Complex_Number ) return Poly;

  -- DESCRIPTION :
  --   Computes the relation between variables i and j
  --   with x(i) = c*t^a, x(j) = d*t^b, e = lcm(|a|,|b|)
  --   where a and b have the same sign.

  function Variable_Relation
              ( n,i,j,a,b : integer32;
                c,d : Complex_Number ) return Poly;

  -- DESCRIPTION :
  --   Computes the relation between variables i and j
  --   with x(i) = c*t^a and x(j) = d*t^b.

  function Paired_Variables ( map : Monomial_Map ) return Poly_List;

  -- DESCRIPTION :
  --   Returns the relations between all variables which share one
  --   and the same parameter.

  function Number_of_Nonzeroes
             ( v : Standard_integer_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nonzero elements in v.

  function Basis_Parameter
             ( map : Monomial_Map; i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   A variable is a basis variable if its coefficient equals one
  --   and if the support of its monomial corresponds to one parameter.
  --   If the i-th variable is not a basis parameter, then 0 is returned,
  --   otherwise the index of the parameter is returned.

  procedure Exponents_of_Basis
               ( map : in Monomial_Map; i,m : in integer32;
                 exp : out Standard_Integer_Matrices.Matrix;
                 pars,vars : out Standard_Integer_Vectors.Vector );
 
  -- DESCRIPTION :
  --   Returns the exponents of the parameters for the i-th variable
  --   in the top column of the m-by-m matrix on return,
  --   with in the other rows the other basis elements.

  -- REQUIRED : m = Number_of_Nonzeroes(map.v(i)) > 1,
  --   exp is an (m+1)-by-m matrix, pars'range = 1..m = vars'range.

  -- ON ENTRY :
  --   map       a monomial map;
  --   i         current variable under consideration;
  --   m         the number of nonzero elements in map.v(i),
  --             which equals the number of parameters used for x(i).

  -- ON RETURN :
  --   exp       first row contains the nonzero elements of map.v(i),
  --             the next rows contains the monomial definitions for
  --             the variables that are basis elements and with their
  --             parameter occurring in map.v(i);
  --   pars      if pars(k) = 1, then the k-th parameter was selected;
  --   vars      if vars(k) = 1, then the k-th variable was selected.

  procedure Multiply ( e : in out Standard_Integer_Matrices.Matrix;
                       f : out integer32 );

  -- DESCRIPTION :
  --   Given an (m+1)-by-m matrix e, the first row is multiplied to
  --   insure that every element on the first row is a multiple of
  --   the corresponding nonzero element in the same column.
  --   On return in f is the accumulated multiplication factor,
  --   as product of the multipliers for the first row.

  function Basis_Equation ( map : Monomial_Map; i : integer32 ) return Poly;

  -- DESCRIPTION :
  --   If the monomial for the i-th variable has more than one parameter
  --   with nonzero power and if there are sufficiently many basis monomials,
  --   then the defining equation is returned.

  function Equations ( map : Monomial_Map ) return Poly_List;
  function Equations ( map : Monomial_Map ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the list of equations for the solution set
  --   given by the parametric representation defined by the map.

end Standard_Monomial_Map_Ideals;
