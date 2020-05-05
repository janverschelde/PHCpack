with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Standard_Coefficient_Convolutions is

-- DESCRIPTION :
--   This package offers more efficient vector computations on 
--   separated real and imaginary parts, splitted from complex vectors.

-- DATA STRUCTURES :

  type VecVecVec is
    array ( integer32 range <> ) of Standard_Floating_VecVecs.Link_to_VecVec;
  -- A three dimensional structure to store the coefficient vectors,
  -- real or imaginary parts, of powers of series.
 
  type Link_to_VecVecVec is access VecVecVec; -- to store the power table

  type VecVecVec_Array is array ( integer32 range <> ) of Link_to_VecVecVec;

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
  --   Stores all powers x(i)^k for k ranging from 2 to d(i),
  --   for i in d'range = x'range, in the power table.
  --   The i-th entry in the power table contains the powers of x(i),
  --   if d(i) > 1, starting with x(i)^2 at the first position.
  --   This Create procedure combines the functions Allocate(mxe,deg) 
  --   and the Compute(rpwt,ipwt,d,rx,ix) procedure below.

-- BASIC COMPUTATIONAL PROCEDURES :

  procedure Update ( rvl,ivl : in Standard_Floating_Vectors.Link_to_Vector;
                     rnc,inc : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds the numbers in rnc to the rvl and 
  --   adds the numbers in inc to the ivl.

  -- REQUIRED : all vectors have the same range.

  procedure Multiply
              ( xr,xi,yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr,zi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficients of the vector x with y,
  --   with real parts in xr, yr, and imaginary parts in xi, yi,
  --   and stores the results in the z, with real and imaginary
  --   parts in zr and zi.

  -- REQUIRED :
  --   All vectors have the same range.

-- COMPUTING THE POWER TABLE :

  procedure Compute ( rpwt,ipwt : in Link_to_VecVecVec;
                      mxe : in Standard_Integer_Vectors.Vector;
                      rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Computes the powers in the allocated power table for the
  --   coefficients in the power series in x.

  -- REQUIRED : rpwt = Allocate(mxe,deg) and ipwt = Allocate(mxe,deg)
  --   have been executed, and rpwt'range = x'range = ipwt'range.

  -- ON ENTRY :
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

  procedure Speel ( rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Speel ( rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rfwd,ifwd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product
  --   of the series in rx, ix, which are all of the same degree.

  -- REQUIRED :
  --   The range of rx and ix starts at 1 and ends at 2 or higher.
  --   All vectors are allocated and of the range 0..degree,
  --   for some same fixed degree.

  -- ON ENTRY :
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
                    rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.Link_to_VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.Link_to_VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluation and differentiation of the sum of products,
  --   given in indexed format at a power series.

  -- ON ENTRY :
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
                    rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rcff,icff : in Standard_Floating_Vectors.Link_to_Vector;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Multiplies the coefficient with the common factor.

  -- REQUIRED : facidx /= null.

  -- ON ENTRY :
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

  -- DESCRIPTION :
  --   Multiplies the coefficients of the power series with multiplier.

  -- ON ENTRY:
  --   multiplier   is the multiplier exponent;
  --   rcff         real parts of coefficients of a power series;
  --   icff         imaginary parts of coefficients of a power series.

  -- ON RETURN :
  --   rcff         real parts of coefficients multiplied;
  --   icff         imaginary parts of coefficients multiplied.

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.Link_to_VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.Link_to_VecVec;
                    rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt,ipwt : in Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Evaluation and differentiation of a polynomial,
  --   given in indexed format at a power series.

  -- ON ENTRY :
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

-- DEALLOCATORS :

  procedure Clear ( pwt : in out VecVecVec );
  procedure Clear ( pwt : in out Link_to_VecVecVec );
  procedure Clear ( pwt : in out VecVecVec_Array );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the power table(s) pwt.

end Standard_Coefficient_Convolutions;
