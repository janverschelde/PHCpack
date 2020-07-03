with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;        use Standard_Floating_VecVecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;

package Standard_Coefficient_Convolutions is

-- DESCRIPTION :
--   This package offers more efficient vector computations on 
--   separated real and imaginary parts, splitted from complex vectors.

-- DATA STRUCTURES :

-- A convolution circuit is a data structure for the efficient evaluation
-- and differentiation of polynomials in several variables at the
-- coefficient vectors of power series using the reverse mode of
-- algorithmic differentiation.

  type Circuit ( nbr,dim,dim1,dim2 : integer32 ) is record
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponent vectors
    idx : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponent indices
    fac : Standard_Integer_VecVecs.VecVec(1..nbr); -- factor indices
   -- the complex coefficients of monomials are stored in parts
    rcf : Standard_Floating_VecVecs.VecVec(1..nbr); -- real parts
    icf : Standard_Floating_VecVecs.VecVec(1..nbr); -- imaginary parts
   -- the complex constant coefficients is stored in parts
    rct : Standard_Floating_Vectors.Link_to_Vector; -- real parts
    ict : Standard_Floating_Vectors.Link_to_Vector; -- imaginary parts
   -- workspace for products and coefficient vectors of series are
   -- stored in real and imaginary parts
    rfwd,ifwd : Standard_floating_VecVecs.VecVec(1..dim1); -- forward products
    rbck,ibck : Standard_floating_VecVecs.VecVec(1..dim2); -- backward products
    rcrs,icrs : Standard_floating_VecVecs.VecVec(1..dim2); -- cross products
   -- real and imaginary parts of workspace and accumulator series
    rwrk,iwrk : Standard_Floating_Vectors.Link_to_Vector;  -- workspace
    racc,iacc : Standard_Floating_Vectors.Link_to_Vector;  -- accumulator
  end record;

  type Link_to_Circuit is access Circuit;

  type Circuits is array ( integer32 range <> ) of Link_to_Circuit;

  type Link_to_Circuits is access Circuits;

-- A system stores the sequence of convolution circuits for each polynomial
-- and the work space to hold auxiliary results and the final outcomes.

  type System ( neq,neq1,dim,dim1,deg : integer32 ) is record
    crc : Circuits(1..neq);    -- circuits for the equations
    mxe : Standard_Integer_Vectors.Vector(1..dim); -- exponent maxima
    rpwt : Link_to_VecVecVec;  -- real parts of the power table
    ipwt : Link_to_VecVecVec;  -- imaginary parts of the power table
   -- work space for EvalDiff on one circuit in ryd and iyd
    ryd : Standard_Floating_VecVecs.VecVec(1..dim1); -- real parts
    iyd : Standard_Floating_VecVecs.VecVec(1..dim1); -- imaginary parts
   -- linearized evaluated power series in vy
    vy : Standard_Complex_VecVecs.VecVec(0..deg);
   -- delinearized evaluated power series in yv
    yv : Standard_Complex_VecVecs.VecVec(1..neq);
   -- differentiation result as matrix series in vm
    vm : Standard_Complex_VecMats.VecMat(0..deg);
  end record;

  type Link_to_System is access System;

  type System_Array is array ( integer32 range<> ) of Link_to_System;

-- CONSTRUCTORS AND ALLOCATORS :

  function Exponent_Maxima
             ( c : Circuits; dim : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the maximal exponents of the dim variables in the circuits.
  --   The result of this function is the mxe for the Allocate and the 
  --   Create of the power table.

  function Allocate ( mxe : Standard_Integer_Vectors.Vector;
                      deg : integer32 )
                    return Link_to_VecVecVec;

  -- DESCRIPTION :
  --   Allocates space for the power table, given the exponent maxima
  --   for each variable in mxe and the degrees of the power series in deg.

  procedure Create ( rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                     mxe : in Standard_Integer_Vectors.Vector;
                     deg : in integer32;
                     rpwt,ipwt : out Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Stores all powers x(i)^k for k ranging from 2 to mxe(i),
  --   for i in mxe'range = rx'range, in the power table.
  --   The i-th entry in the power table contains the powers of x(i),
  --   if mxe(i) > 1, starting with x(i)^2 at the first position.
  --   This Create procedure combines the functions Allocate(mxe,deg) above
  --   and the Compute(rpwt,ipwt,mxe,rx,ix) procedure below.

  function Linearized_Allocation
             ( dim,deg : integer32 )
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of the series
  --   truncated to degree deg and initialized to zero, in linearized form.
  --   The vector on return has range 0..deg and represents a series
  --   truncated to degree deg and with vectors of range 1..dim as
  --   its coefficients.

  function Allocate_Coefficients
             ( nbq,nvr,deg : integer32 )
             return Standard_Complex_VecMats.VecMat;

  -- DESCRIPTION :
  --   Returns a vector of matrices, of range 0..deg,
  --   of nbq-by-nvr matrices of coefficients,
  --   where nbq equals the number of equations
  --   and nvr is the number of variables.

  function Create ( c : Circuits; dim,deg : integer32 ) return System;
  function Create ( c : Circuits; dim,deg : integer32 ) return Link_to_System;

  -- DESCRIPTION:
  --   The system on return stores the convolution circuits in crc,
  --   contains the values for the exponent maxima in mxe, and
  --   has allocated space for rpwt, ipwt, ryd, iyd, vy, yv, and vm.

-- COPY PROCEDURES :

  procedure Copy ( c_from : in Circuit; c_to : out Circuit );

  -- DESCRIPTION :
  --   Copies the circuit c_from to the circuit c_to.

  -- REQUIRED :
  --   c_from.nbr = c_to.nbr, c_from.dim = c_to.dim,
  --   c_from.dim1 = c_to.dim1, and c_from.dim2 = c_to.dim2.

  procedure Copy ( c_from : in Link_to_Circuit; c_to : out Link_to_Circuit );

  -- DESCRIPTION :
  --   Copies the circuit c_from to the circuit c_to.
  --   Deallocates c_to before the copy.

  procedure Copy ( c_from : in Circuits; c_to : out Circuits );

  -- DESCRIPTION :
  --   Copies all circuits from c_from to c_to.

  -- REQUIRED : c_from'range = c_to'range.

  procedure Copy ( s_from : in System; s_to : out System );

  -- DESCRIPTION :
  --   Copies the system s_from to the system s_to.

  -- REQUIRED :
  --   s_from.neq = s_to.neq, s_from.neq1 = s_to.neq1,
  --   s_from.dim = s_to.dim, s_from.dim1 = s_to.dim1, and
  --   s_from.deg = s_to.deg.

  procedure Copy ( s_from : in Link_to_System; s_to : out Link_to_System );

  -- DESCRIPTION :
  --   Copies the system s_from to the system s_to.
  --   The system s_to is deallocated before the copy.

-- BASIC COMPUTATIONAL PROCEDURES :

  procedure Update ( rvl,ivl : in Standard_Floating_Vectors.Link_to_Vector;
                     rnc,inc : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Update ( deg : in integer32;
                     rvl,ivl : in Standard_Floating_Vectors.Link_to_Vector;
                     rnc,inc : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds the numbers in rnc to the rvl and 
  --   adds the numbers in inc to the ivl.
  --   If the degree deg is provided,
  --   then the update is restricted to the range 0..deg.

  -- REQUIRED : all vectors have the same range,
  --   or include the range 0..deg.

  procedure Multiply
              ( xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Multiply
              ( deg : in integer32;
                xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficients of the vector x with y,
  --   with real parts in xr, yr, and imaginary parts in xi, yi,
  --   and stores the results in the z, with real and imaginary
  --   parts in zr and zi.  If the degree deg is provided,
  --   then the convolution is restricted up to the degree deg.

  -- REQUIRED :
  --   All vectors have the same range, or include the range 0..deg.

-- COMPUTING THE POWER TABLE :

  procedure Compute ( rpwt,ipwt : in Link_to_VecVecVec;
                      mxe : in Standard_Integer_Vectors.Vector;
                      rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Compute ( deg : in integer32;
                      rpwt,ipwt : in Link_to_VecVecVec;
                      mxe : in Standard_Integer_Vectors.Vector;
                      rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Computes the powers in the allocated power table for the
  --   coefficients in the power series in x.

  -- REQUIRED : rpwt = Allocate(mxe,deg) and ipwt = Allocate(mxe,deg)
  --   have been executed, and rpwt'range = x'range = ipwt'range,
  --   or at least include the range 0..deg.

  -- ON ENTRY :
  --   deg      (optional) degree of the power series, if provided,
  --            then the convolutions are limited to degree deg,
  --            otherwise the convolutions run till the end of the vectors;
  --   rpwt     allocated power table for the exponent maxima in mxe
  --            and for the real parts of the coefficients of series
  --            of the degree of the series in x;
  --   ipwt     allocated power table for the exponent maxima in mxe
  --            and for the imaginary parts of the coefficients of series
  --            of the degree of the series in x;
  --   mxe      exponent maxima for each variable in the system;
  --   rx       real parts of coefficients of power series;
  --   ix       imaginary parts of coefficients of power series.

  -- ON RETURN :
  --   rpwt     real parts of the coefficients of the power table;
  --   ipwt     imaginary parts of the coefficients of the power table.

-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel ( rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( deg : in integer32;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( rx,ix : in Standard_Floating_VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( deg : in integer32;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product
  --   of the series in rx, ix, which are all of the same degree.

  -- REQUIRED :
  --   The range of rx and ix starts at 1 and ends at 2 or higher.
  --   All vectors are allocated and of the range 0..degree,
  --   for some same fixed degree.

  -- ON ENTRY :
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   rx           a vector with range starting at 1, ending at 2 or higher,
  --                with the real parts of the coefficients of series;
  --   ix           a vector with range starting at 1, ending at 2 or higher,
  --                with the imaginary parts of the coefficients of series;
  --   idx          if provided, then only those indices of x with values
  --                in idx participate and dim = idx'last,
  --                otherwise, all values of x participate and dim = x'last;
  --   rfwd,ifwd    space allocated for dim-1 coefficient vectors,
  --                for series of the same fixed degree;
  --   backward     has space reserved for dim-2 coefficient vectors,
  --                for series of the same fixed degree;
  --   cross        has space reserved for dim-2 coefficient vectors,
  --                for series of the same fixed degree.

  -- ON RETURN :
  --   rfwd         accumulates the real parts of the forward products,
  --                rfwd(dim-1) holds the coefficients for the product,
  --                rfwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   ifwd         accumulates the imaginary parts of the forward products,
  --                ifwd(dim-1) holds the coefficients for the product,
  --                ifwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   rbck         accumulates the real parts of the backward products,
  --                rbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   ibck         accumulates the imaginary parts of the backward products,
  --                rbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   rcrs         stores the real parts of the cross products,
  --                rcrs(k) contains the real parts of the coefficients
  --                of the partial derivative w.r.t. k+1;
  --   icrs         stores the imaginary parts of the cross products,
  --                icrs(k) contains the imaginary parts of coefficients
  --                of the partial derivative w.r.t. k+1.

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( deg : in integer32;
                    idx : in Standard_Integer_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Speel ( deg : in integer32;
                    idx : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluation and differentiation of the sum of products,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   idx          indexed representation of a sum of products of variables;
  --   rcff         real parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   icff         imaginary parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   rx           real parts of coefficients of series of same degree;
  --   ix           imaginary parts of coefficients of series of same degree;
  --   rfwd         work space allocated for rx'last-1 coefficient vectors
  --                of the same fixed degree as the series in rx;
  --   ifwd         work space allocated for ix'last-1 coefficient vectors
  --                of the same fixed degree as the series in ix;
  --   rbck         work space allocated for rx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rx;
  --   ibck         work space allocated for ix'last-2 coefficient vectors
  --                of the same fixed degree as the series in ix;
  --   rcrs         work space allocated for rx'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   icrs         work space allocated for ix'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   ryd          vector of range 0..rx'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   iyd          vector of range 0..ix'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rwrk         work space for the coefficients of the same fixed degree;
  --   iwrk         work space for the coefficients of the same fixed degree.

  -- ON RETURN :
  --   ryd          ryd(rx'last+1) contains the real parts of coefficients
  --                of the value of the sum of products, evaluated at x,
  --                ryd(k) is the real part of the k-th partial derivative;
  --   iyd          iyd(ix'last+1) contains the real parts of coefficients
  --                of the value of the sum of products, evaluated at x,
  --                iyd(k) is the real part of the k-th partial derivative.

  procedure Multiply_Factor
                  ( xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_Vectors.Link_to_Vector;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec );
  procedure Multiply_Factor
                  ( deg : in integer32;
                    xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_Vectors.Link_to_Vector;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Multiplies the coefficient with the common factor.

  -- REQUIRED : facidx /= null.

  -- ON ENTRY :
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   xpk          k-th exponent vector;
  --   facidx       factor index of k-th exponents in xpk;
  --   rx           real parts of the values for the variables;
  --   ix           imaginary parts of the values for the variables;
  --   rcff         real parts of the coefficient of the monomial;
  --   icff         imaginary parts coefficient of the monomial;
  --   rwrk         allocated space for the real parts
  --                of the coefficients of series of same degree;
  --   iwrk         allocated space for the imaginary parts
  --                of the coefficients of series of same degree;
  --   racc         allocated space for the real parts 
  --                of coefficients of series of same degree;
  --   iacc         allocated space for the imaginary parts 
  --                of coefficients of series of same degree;
  --   rpwt         power table of the real parts for the values in x;
  --   ipwt         power table of the real imaginary for the values in x.

  -- ON RETURN :
  --   racc         accumulates the product of the real parts
  --                of the coefficients with the evaluated powers of x
  --                as defined by xpk and the factor index in facidx.
  --   iacc         accumulates the product of the imaginary parts
  --                of the coefficients with the evaluated powers of x
  --                as defined by xpk and the factor index in facidx.

  procedure Multiply_Power
                  ( multiplier : in integer32;
                    rcff : in Standard_Floating_Vectors.Link_to_Vector; 
                    icff : in Standard_Floating_Vectors.Link_to_Vector ); 
  procedure Multiply_Power
                  ( deg,multiplier : in integer32;
                    rcff : in Standard_Floating_Vectors.Link_to_Vector; 
                    icff : in Standard_Floating_Vectors.Link_to_Vector ); 

  -- DESCRIPTION :
  --   Multiplies the coefficients of the power series with multiplier.

  -- ON ENTRY:
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   multiplier   is the multiplier exponent;
  --   rcff         real parts of coefficients of a power series;
  --   icff         imaginary parts of coefficients of a power series.

  -- ON RETURN :
  --   rcff         real parts of coefficients multiplied;
  --   icff         imaginary parts of coefficients multiplied.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec );
  procedure Speel ( deg : in integer32;
                    xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    rx,ix : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Evaluation and differentiation of a polynomial,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   xps          exponent vector of the monomials;
  --   idx          indexed representation of the variables in exponents;
  --   fac          factor index of the exponents;
  --   rcff         real parts of coefficients of the products;
  --   icff         imaginary parts of coefficients of the products;
  --   rx           real parts of coefficients of series of same degree;
  --   ix           imaginary parts of coefficients of series of same degree;
  --   rfwd         work space allocated for rx'last-1 coefficient vectors
  --                of the same fixed degree as the series in rx;
  --   ifwd         work space allocated for ix'last-1 coefficient vectors
  --                of the same fixed degree as the series in ix;
  --   rbck         work space allocated for rx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rx;
  --   ibck         work space allocated for ix'last-2 coefficient vectors
  --                of the same fixed degree as the series in ix;
  --   rcrs         work space allocated for rx'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   icrs         work space allocated for ix'last-2 coefficient vectors
  --                of the same fixed degree as the series in x;
  --   ryd          vector of range 0..rx'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   iyd          vector of range 0..ix'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rwrk         work space for the coefficients of the same fixed degree;
  --   iwrk         work space for the coefficients of the same fixed degree;
  --   racc         work space for the coefficients of the same fixed degree;
  --   iacc         work space for the coefficients of the same fixed degree;
  --   rpwt         power table of the real parts for the values in x;
  --   ipwt         power table of the real imaginary for the values in x.

  -- ON RETURN :
  --   ryd          ryd(rx'last+1) contains the real parts of the coefficient
  --                vector of the value of the sum of products evaluated at x,
  --                ryd(k) is the real part of the k-th partial derivative;
  --   iyd          ryd(ix'last+1) contains the real parts of the coefficient
  --                vector of the value of the sum of products evaluated at x,
  --                iyd(k) is the real part of the k-th partial derivative.

-- EVALUATION OF POWER SERIES COEFFICIENTS :

  procedure EvalCoeff ( c : in Circuit; t : in double_float;
                        rct,ict : out double_float;
                        rcf : out Standard_Floating_Vectors.Vector;
                        icf : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   The coefficients of the power series in the circuit c 
  --   are evaluated at t.

  -- REQUIRED : rcf'range = icf'range = 1..c.nbr.

  -- ON ENTRY :
  --   c        a coefficient circuit;
  --   t        value for the continuation parameter.

  -- ON RETURN :
  --   rct      real part of the constant coefficient of c, at t;
  --   ict      imaginary part of the constant coefficient of c, at t;
  --   rcf      real parts of the coefficients of c, evaluated at t;
  --   icf      imaginary parts of the coefficients of c, evaluated at t.

-- EVALUATION AND DIFFERENTIATION ON CIRCUITS :

  procedure EvalDiff ( c : in Circuit;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec );
  procedure EvalDiff ( deg : in integer32; c : in Circuit;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Wraps the Speel procedure for the convolution circuit c, to
  --   evaluate at the series with real and imaginary parts in rx and ix,
  --   with the aid of the power table, with real and imaginary parts
  --   in rpwt and ipwt.  The result are placed in ryd and iyd.

  -- ON ENTRY :
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   c            a circuit properly defined and allocated;
  --   rx           real parts of coefficients of series of same degree;
  --   ix           imaginary parts of coefficients of series of same degree;
  --   ryd          vector of range 0..rx'last+1 with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   iyd          vector of range 0..ix'last+1 with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rpwt         power table of the real parts for the values in x;
  --   ipwt         power table of the real imaginary for the values in x.

  -- ON RETURN :
  --   ryd          ryd(rx'last+1) contains the real parts of the coefficient
  --                vector of the value of the sum of products evaluated at x,
  --                ryd(k) is the real part of the k-th partial derivative;
  --   iyd          ryd(ix'last+1) contains the real parts of the coefficient
  --                vector of the value of the sum of products evaluated at x,
  --                iyd(k) is the real part of the k-th partial derivative.

  procedure EvalDiff ( c : in Circuits;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                       vy : in Standard_Complex_VecVecs.VecVec;
                       vm : in Standard_Complex_VecMats.VecMat );
  procedure EvalDiff ( deg : in integer32; c : in Circuits;
                       rx : in Standard_Floating_VecVecs.VecVec;
                       ix : in Standard_Floating_VecVecs.VecVec;
                       rpwt,ipwt : in Link_to_VecVecVec;
                       ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                       vy : in Standard_Complex_VecVecs.VecVec;
                       vm : in Standard_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Evaluates and differentiates the convolution circuits in c
  --   at the series x.

  -- ON ENTRY :
  --   deg          (optional) degree of the power series,
  --                if provided, then the convolutions will stop at deg,
  --                otherwise continue till the end of the vectors;
  --   c            an array of convolution circuits;
  --   rx           real parts of coefficients of series of same degree;
  --   ix           imaginary parts of coefficients of series of same degree;
  --   rpwt         power table of the real parts for the values in x;
  --   ipwt         power table of the real imaginary for the values in x;
  --   ryd          work space of range 1..rx'last+1 to contain the real
  --                parts of the gradient and the value of the circuits in c;
  --   iyd          work space of range 1..ix'last+1 to contain the imaginary
  --                parts of the gradient and the value of the circuits in c;
  --   vy           allocated space for the values of the circuits at x,
  --                done by the above procedure Linearized_Allocation,
  --   vm           space allocated for a series of some fixed degree
  --                with matrix coefficients.

  -- ON RETURN :
  --   vy           values of the circuits at x, in linearized form;
  --   vm           the evaluated circuits at x as a series 
  --                of some fixe degree with matrix coefficients.

  procedure Delinearize ( vy,yv : in Standard_Complex_VecVecs.VecVec );
  procedure Delinearize ( deg : in integer32;
                          vy,yv : in Standard_Complex_VecVecs.VecVec );

  --  DESCRIPTION :
  --    Converts the linearized representation in vy into the vector yv
  --    of coefficient vectors of the series.
  --    This conversion is convenient for the difference computation
  --    and needed in the application of Newton's method.
  --    If the degree deg is provided, then the conversion stops at deg,
  --    otherwise, the conversion runs over all coefficients in the series.

  -- REQUIRED :
  --   if vy'range = 0..degree and vy(k)'range = 1..dimension,
  --   then yv'range = 1..dimension and yv(k)'range = 0..degree.

  procedure EvalDiff ( s : in System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec );
  procedure EvalDiff ( deg : in integer32; s : in System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec );
  procedure EvalDiff ( s : in Link_to_System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec );
  procedure EvalDiff ( deg : in integer32; s : in Link_to_System;
                       rx,ix : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Wraps the EvalDiff on the convolution circuits in s.crc,
  --   at the power series with coefficients in x,
  --   with real and imaginary parts in rx and ix.
  --   If the degree deg is provided, then the convolutions stop at deg,
  --   otherwise, the convolutions run over all coefficients in the series.

  -- REQUIRED :
  --   All data in s are allocated properly with respect to dimension
  --   and degree, the power table s.pwt is up to data with the given
  --   real and imaginary parts in rx and ix.

-- DEALLOCATORS :

  procedure Clear ( c : in out Circuit );
  procedure Clear ( c : in out Link_to_Circuit );
  procedure Clear ( c : in out Circuits );
  procedure Clear ( c : in out Link_to_Circuits );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the convolution circuits.

  procedure Clear ( s : in out System );
  procedure Clear ( s : in out Link_to_System );
  procedure Clear ( s : in out System_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the system (array) s.

end Standard_Coefficient_Convolutions;
