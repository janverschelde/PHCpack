with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;

package QuadDobl_Predictor_Convolutions is

-- DESCRIPTION :
--   Exports data structures and procedures to construct a rational prediction
--   for the next solution on a path defined by a polynomial homotopy,
--   represented by convolution circuits, in quad double precision.

-- DATA STRUCTURE :
--   The predictor type stores a solution in a vector of dimension dim,
--   of power series of degree deg, for a rational approximation of degrees
--   numdeg and dendeg, respectively of numerator and denominator.
--   The work space for the pivoting information of the linear system
--   solving of the series computation and the rational approximation.
--   The work space (wrk below) for Newton's method on power series
--   has a dimension equal to the number of equations in the homotopy.

  type LU_Predictor ( dim,deg,numdeg,dendeg : integer32 ) is record
    sol : QuadDobl_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector;        -- Newton work space
    numcff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    newtpiv : Standard_Integer_Vectors.Vector(1..dim);    -- pivots for Newton
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for rational approximation
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
  end record;
  type Link_to_LU_Predictor is access LU_Predictor;

-- DATA STRUCTURE FOR SVD :
--   The predictor with the Singular Value Decomposition works for systems
--   that are overdetermined, with more equations than variables.
--   The number of equations is stored in the neq attribute,
--   the number of variables is stored in the dim attribute.
--   The number of variables plus one is the dim1 attribute.
--   The update in Newton's method is stored twice, once as dx,
--   as a vector of power series, and once as xd, in linearized format,
--   as a series with vector coefficients. 

  type SVD_Predictor ( neq,dim,dim1,deg,numdeg,dendeg : integer32 ) is record
    sol : QuadDobl_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector;        -- of range 1..neq
    dx : QuadDobl_Complex_VecVecs.VecVec(1..dim);         -- vector of series
    xd : QuadDobl_Complex_VecVecs.VecVec(0..deg);         -- series of vectors
    svl : QuadDobl_Complex_Vectors.Vector(1..dim1);       -- singular values
    U : QuadDobl_Complex_Matrices.Matrix(1..neq,1..neq);  -- U in SVD
    V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);  -- V in SVD
    ewrk : QuadDobl_Complex_Vectors.Link_to_Vector;       -- of range 1..dim
    numcff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : QuadDobl_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for the rational approximation
    rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg); -- right hand side
  end record;
  type Link_to_SVD_Predictor is access SVD_Predictor;

  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) return LU_Predictor;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return Link_to_LU_Predictor;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32 ) return SVD_Predictor;
  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return Link_to_SVD_Predictor;

  -- DESCRIPTION :
  --   Given a solution vector and dimensions,
  --   makes the power series for sol and allocate work space.

  -- ON ENTRY :
  --   sol      a solution vector;
  --   neq      number of equations in the homotopy;
  --   deg      degree of the power series;
  --   numdeg   degree of the numerator of the rational approximation;
  --   dendeg   degree of the denominator of the rational approximation.

  -- ON RETURN :
  --   data allocated to run the Newton-Fabry-Pade predictor.

  procedure Clear ( p : in out Link_to_LU_Predictor );
  procedure Clear ( p : in out Link_to_SVD_Predictor );

  -- DESCRIPTION :
  --   Deallocates the memory for the series, the work space,
  --   and the rational approximation.

end QuadDobl_Predictor_Convolutions;
