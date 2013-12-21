with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
--with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;

package Affine_Binomial_Iterator is

-- DESCRIPTION :
--   The focus of this package is to enumerate all factors that
--   contribute to a generalized permanent of the incidence matrix
--   of the bipartite graph linking monomials to variables.

  procedure Initialize_Subsets ( n,max : in integer32 );

  -- DESCRIPTION :
  --   Initializes the state of the iterator to enumerate all subsets of
  --   a set of n variables.  The parameter max determines the maximum 
  --   number of variables to be set to zero.

  procedure Initialize_Iterator 
              ( A : in Standard_Integer_Matrices.Matrix;
                max : in integer32 );

  -- DESCRIPTION :
  --   The incidence matrix is stored in the iterator to enumerate all
  --   candidate affine solution sets.  The parameter max limits the
  --   cardinality of the set of selected variables to be set to zero.

  function Set_to_Zero
              ( A : Standard_Integer_Matrices.Matrix; i : integer32;
                s : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the selection of variables set to zero
  --   make the i-th monomial defined by the i-th row of A zero.

  -- ON ENTRY :
  --   A        the incidence matrix of the bipartite graph:
  --            A[i,j] > 0 means that x[j] occurs with positive
  --            powers in the i-th monomial;
  --   i        if i is even, then monomial i occurs as first monomial
  --            in the binomial i/2 in the system, otherwise the
  --            monomial is second in the binomial equation i/2;
  --   s        the current selection of variables to be set to zero.

  function All_Set_to_Zero 
             ( A : Standard_Integer_Matrices.Matrix;
               s : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the selection of the variables in s causes
  --   all monomials in A to vanish.

  function Set_Nonzero
              ( A : Standard_Integer_Matrices.Matrix; i : integer32;
                s : Standard_Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the selection of variables that may not be zero
  --   prevents the i-th monomial from every becoming zero.

  -- ON ENTRY :
  --   A        the incidence matrix of the bipartite graph:
  --            A[i,j] > 0 means that x[j] occurs with positive
  --            powers in the i-th monomial;
  --   i        if i is even, then monomial i occurs as first monomial
  --            in the binomial i/2 in the system, otherwise the
  --            monomial is second in the binomial equation i/2;
  --   s        the current selection of variables which may not be zero.

  procedure Update_Present_Variables
              ( s : in out Standard_Integer_Vectors.Vector;
                A : in Standard_Integer_Matrices.Matrix; i : in integer32 );

  -- DESCRIPTION :
  --   Sets s(j) to 1 if A(i,j) > 0.  The purpose of s is to record
  --   which variables should not be set to zero because they occur in
  --   a monomial of an equation that will be skipped in the enumeration.
 
  generic

    with procedure Report
           ( setzero : in Standard_Integer_Vectors.Vector;
             cntzero : in integer32;
             nonzero : in Standard_Integer_Vectors.Vector;
            -- resteqs : in Standard_Integer_Vectors.Vector;
             cntrest : in integer32;
             continue : out boolean );

    -- DESCRIPTION :
    --   Each time a new choice is found, Report is called.
   
    -- ON ENTRY :
    --   setzero    indices of variables selected to be zero;
    --   cntzero    number of variables selected to be zero;
    --   nonzero    indices of variables that should not be zero,
    --              because of skipped binomial equations.
    --   resteqs    the k-th equation is remaining if resteqs(k) = 1,
    --              if resteqs(k) = 0, then the k-th equation vanishes;
    --   cntrest    number of remaining equations.
  
    -- ON RETURN :
    --   continue   true if the enumeration should go on,
    --              false if the enumerate must terminate.

  procedure Enumerate
              ( A : in Standard_Integer_Matrices.Matrix;
                s0_max : in integer32 );

  -- DESCRIPTION :
  --   Enumerates choices of variables to be set to zero,
  --   going over the rows of A.

  -- ON ENTRY :
  --   A        rows of A are indexed by the monomials,
  --            the columns by the variables 
  --   A[i,j] = 1 : the j-th variable has positive power in monomial i,
  --          = 0 : zeroing the j-th variable does nothing to monomial i;
  --   s0_max   maximum number of variables which may be set to zero,
  --            to limit the depth of the enumeration tree.

  procedure Next_Subset
              ( setzero : out Standard_Integer_Vectors.Vector;
                cntzero : out integer32 );

  -- DESCRIPTION :
  --   Returns the next choice of variables to be set to zero,
  --   considering all subsets of cardinality <= max,
  --   where max was given to the Initialize_Subsets procedure.

  -- REQUIRED : Initialize_Subsets was executed.
   
  -- ON ENTRY :
  --   setzero  indices of variables selected to be zero;
  --   cntzero  number of variables selected to be zero;
  --            if -1, then the iteration stops.

  procedure Next_Selection
              ( setzero : out Standard_Integer_Vectors.Vector;
                cntzero : out integer32;
                nonzero : out Standard_Integer_Vectors.Vector;
                cntrest : out integer32 );

  -- DESCRIPTION :
  --   Returns the next choice of variables to be set to zero,
  --   taking into account the incidence matrix given to the
  --   procedure Initialize_Iterator.
   
  -- REQUIRED : Initialize_Iterator was executed.

  -- ON ENTRY :
  --   setzero  indices of variables selected to be zero;
  --   cntzero  number of variables selected to be zero;
  --            if -1, then the iteration stops.
  --   nonzero  indices of variables that should not be zero,
  --            because of skipped equations;
  --   cntrest  number of remaining equations.

  -- function Generated_Selections return List;

  -- DESCRIPTION :
  --   Returns the list of generated selections of subsets of variables,
  --   This list is maintained to avoid duplicate reporting and is only
  --   complete at the end of the iterator.

  -- NOTE : for large enumerations searching each time through that list
  --   could be prohibitively inefficient, although with a hash table,
  --   the generated selections could work well.

  procedure Clear_Subsets;

  -- DESCRIPTION :
  --   Clears the state of the iterator to enumerate all subsets.
 
  procedure Clear_Iterator;

  -- DESCRIPTION :
  --   Clears the iterator to enumerate all candidate affine solution sets.

end Affine_Binomial_Iterator;
