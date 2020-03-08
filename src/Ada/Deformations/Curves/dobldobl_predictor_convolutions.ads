with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;

package DoblDobl_Predictor_Convolutions is

-- DESCRIPTION :
--   Exports data structures and procedures to construct a rational prediction
--   for the next solution on a path defined by a polynomial homotopy,
--   represented by convolution circuits, in double double precision.

-- DATA STRUCTURE :
--   The predictor type stores a solution in a vector of dimension dim,
--   of power series of degree deg, for a rational approximation of degrees
--   numdeg and dendeg, respectively of numerator and denominator.
--   The work space for the pivoting information of the linear system
--   solving of the series computation and the rational approximation.
--   The work space (wrk below) for Newton's method on power series
--   has a dimension equal to the number of equations in the homotopy.

  type Predictor ( dim,deg,numdeg,dendeg : integer32 ) is record
    sol : DoblDobl_Complex_VecVecs.VecVec(1..dim);        -- degree deg series
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector;        -- Newton work space
    numcff : DoblDobl_Complex_VecVecs.VecVec(1..dim);     -- numerator coeffs
    dencff : DoblDobl_Complex_VecVecs.VecVec(1..dim);     -- denominator coeffs
    newtpiv : Standard_Integer_Vectors.Vector(1..dim);    -- pivots for Newton
    padepiv : Standard_Integer_Vectors.Vector(1..dendeg); -- pivots for Pade
    mat : Matrix(1..dendeg,1..dendeg); -- matrix for rational approximation
    rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
  end record;

  type Link_to_Predictor is access Predictor;

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32 ) return Predictor;
  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32 )
                  return Link_to_Predictor;

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

  procedure Clear ( p : in out Link_to_Predictor );

  -- DESCRIPTION :
  --   Deallocates the memory for the series and rational approximation.

end DoblDobl_Predictor_Convolutions;
