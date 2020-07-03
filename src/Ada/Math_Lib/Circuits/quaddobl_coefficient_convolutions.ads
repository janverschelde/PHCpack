with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;

package QuadDobl_Coefficient_Convolutions is

-- DESCRIPTION :
--   The convolutions on power series in quad double precision
--   work on complex vectors in a 2-vector representation,
--   defined by the real parts and the imaginary parts,
--   stored in four consecutive doubles, for the highest, second highest,
--   second lowest and lowest double for each quad double number.
--   The floating vectors are thus four times as long as the complex vectors.

-- DATA STRUCTURES :

-- A convolution circuit is a data structure for the efficient evaluation
-- and differentiation of polynomials in several variables at the
-- coefficient vectors of power series using the reverse mode of
-- algorithmic differentiation.

  type Circuit ( nbr,dim,dim1,dim2 : integer32 ) is record
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponent vectors
    idx : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponent indices
    fac : Standard_Integer_VecVecs.VecVec(1..nbr); -- factor indices
   -- the complex coefficients of monomials are stored in 2 vectors
    rcff : Standard_Floating_VecVecs.VecVec(1..nbr); -- real parts
    icff : Standard_Floating_VecVecs.VecVec(1..nbr); -- imaginary parts
   -- the complex constant coefficient is stored in 2 vectors
    rcst : Standard_Floating_Vectors.Link_to_Vector; -- real parts
    icst : Standard_Floating_Vectors.Link_to_Vector; -- imaginary parts
   -- workspace for products and coefficient vectors of series are
   -- stored in real and imaginary parts, in long vectors
   -- the forward products in 2-vector representation :
    rfwd,ifwd : Standard_floating_VecVecs.VecVec(1..dim1);
   -- the backward products in 2-vector representation :
    rbck,ibck : Standard_floating_VecVecs.VecVec(1..dim2);
   -- the cross products in 2-vector representation :
    rcrs,icrs : Standard_floating_VecVecs.VecVec(1..dim2);
   -- real and imaginary parts of workspace and accumulator series
    rwrk,iwrk : Standard_Floating_Vectors.Link_to_Vector;
    racc,iacc : Standard_Floating_Vectors.Link_to_Vector;
  end record;

  type Link_to_Circuit is access Circuit;

  type Circuits is array ( integer32 range <> ) of Link_to_Circuit;

  type Link_to_Circuits is access Circuits;

-- A system stores the sequence of convolution circuits for each polynomial
-- and the work space to hold auxiliary results and the final outcomes.

  type System ( neq,neq1,dim,dim1,deg : integer32 ) is record
    crc : Circuits(1..neq);    -- circuits for the equations
    mxe : Standard_Integer_Vectors.Vector(1..dim); -- exponent maxima
   -- the power table is stored in two parts
    rpwt : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    ipwt : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
   -- work space for EvalDiff on one circuit in ryd and iyd
    ryd : Standard_Floating_VecVecs.VecVec(1..dim1); -- real parts
    iyd : Standard_Floating_VecVecs.VecVec(1..dim1); -- imaginary parts
   -- linearized evaluated power series in vy
    vy : QuadDobl_Complex_VecVecs.VecVec(0..deg);
   -- delinearized evaluated power series in yv
    yv : QuadDobl_Complex_VecVecs.VecVec(1..neq);
   -- differentiation result as matrix series in vm
    vm : QuadDobl_Complex_VecMats.VecMat(0..deg);
  end record;

  type Link_to_System is access System;

  type System_Array is array ( integer32 range<> ) of Link_to_System;

-- ALLOCATORS AND CONSTRUCTORS :

  function Exponent_Maxima
             ( c : Circuits; dim : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the maximal exponents of the dim variables in the circuits.
  --   The result of this function is the mxe for the Power Table computation.

  function Linearized_Allocation
             ( dim,deg : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of the series
  --   truncated to degree deg and initialized to zero, in linearized form.
  --   The vector on return has range 0..deg and represents a series
  --   truncated to degree deg and with vectors of range 1..dim as
  --   its coefficients.

  function Allocate_Coefficients
             ( nbq,nvr,deg : integer32 )
             return QuadDobl_Complex_VecMats.VecMat;

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
  --   has allocated space for rhpwt, ihpwt, rlpwt, ilpwt, rhyd, ihyd,
  --   rlyd, ilyd, vy, yv, and vm.

-- COMPUTING THE POWER TABLE :

  procedure Compute
              ( rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes the powers in the allocated power table for the
  --   coefficients in the power series in x.

  -- REQUIRED : rpwt = Standard_Floating_VecVecVecs.Allocate(mxe,deg)
  --   and ipwt = Standard_Floating_VecVecVecs.Allocate(mxe,deg)
  --   have been executed, and *pwt'range = x'range.

  -- ON ENTRY :
  --   rpwt     allocated power table for the exponent maxima in mxe
  --            and for the real parts of the coefficients of series,
  --            of the degree of the series in x;
  --   ipwt     allocated power table for the exponent maxima in mxe
  --            and for the imaginary parts of the coefficients of series,
  --            of the degree of the series in x;
  --   mxe      exponent maxima for each variable in the system;
  --   xr       real parts of coefficients of power series;
  --   xi       imaginary parts of coefficients of power series;
  --   u,v,w    are work space vector of range 0..3.

  -- ON RETURN :
  --   rpwt     real high parts of the coefficients of the power table;
  --   ipwt     imaginary high parts of the coefficients of the power table.

-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel ( xr,xi : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Speel ( xr,xi : in Standard_Floating_VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product of the
  --   series in xr and xi which are all of the same degree.

  -- REQUIRED :
  --   The range of xr and xi starts at 1 and ends at 2 or higher.
  --   All vectors are allocated and of the range 0..degree,
  --   for some same fixed degree.

  -- ON ENTRY :
  --   xr           a vector with range starting at 1, ending at 2 or higher,
  --                with the real high parts of the coefficients of series;
  --   xi           a vector with range starting at 1, ending at 2 or higher,
  --                with the imaginary high parts of the coefficients;
  --   idx          if provided, then only those indices of x with values
  --                in idx participate and dim = idx'last,
  --                otherwise, all values of x participate and dim = x'last;
  --   rfwd,ifwd    space allocated for dim-1 coefficient vectors,
  --                for real and imaginary parts of series coefficients;
  --   rbck,ibck    has space reserved for dim-2 coefficient vectors,
  --                for real and imaginary parts of series coefficients;
  --   rcrs,icrs    has space reserved for dim-2 coefficient vectors,
  --                for real and imaginary parts of series coefficients;
  --   u,v,w        are work space vectors of range 0..3.

  -- ON RETURN :
  --   rfwd         accumulates the real parts of the forward products,
  --                rfwd(dim-1) holds the coefficients for the product,
  --                rfwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   ifwd         accumulates the imaginary parts of forward products,
  --                ifwd(dim-1) holds the coefficients for the product,
  --                ifwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   rbck         accumulates the real parts of the backward products,
  --                rbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   ibck         accumulates the imaginary parts of backward products,
  --                rbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   rcrs         stores the real parts of the cross products,
  --                rcrs(k) contains the real parts of the coefficients
  --                of the partial derivative w.r.t. k+1;
  --   icrs         stores the imaginary parts of the cross products,
  --                icrs(k) contains the imaginary parts of coefficients
  --                of the partial derivative w.r.t. k+1.

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    xr,xi : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rcff,icff : in Standard_Floating_VecVecs.VecVec;
                    xr,xi : in Standard_Floating_VecVecs.VecVec;
                    rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                    rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                    rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                    ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                    rwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluation and differentiation of the sum of products,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   idx          indexed representation of a sum of products of variables;
  --   rcff         real parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   icff         imaginary parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   xrx          real parts of coefficients of series of same degree;
  --   xix          imaginary parts of coefficients of series;
  --   rfwd         work space allocated for xr'last-1 coefficient vectors
  --                of the same fixed degree as the series in rhx;
  --   ifwd         work space allocated for xi'last-1 coefficient vectors
  --                of the same fixed degree as the series in ihx;
  --   rbck         work space allocated for xr'last-2 coefficient vectors
  --                of the same fixed degree as the series in rhx;
  --   ibck         work space allocated for xi'last-2 coefficient vectors
  --                of the same fixed degree as the series in ihx;
  --   rcrs         work space allocated for rhx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rhx;
  --   icrs         work space allocated for ihx'last-2 coefficient vectors
  --                of the same fixed degree as the series in ihx;
  --   rhyd         vector of range 0..xr'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   ihyd         vector of range 0..xi'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rwrk         work space for the coefficients of the same fixed degree;
  --   iwrk         work space for the coefficients of the same fixed degree;
  --   u,v,w        are work space vectors of range 0..3.

  -- ON RETURN :
  --   ryd          ryd(xrx'last+1) contains the real parts of
  --                coefficients of the value of the sum of products at x;
  --                ryd(k) is the real high part of the k-th partial
  --                derivative at x;
  --   iyd          iyd(xi'last+1) contains the imaginary parts of
  --                coefficients of the value of the sum of products at x,
  --                iyd(k) is the imaginary high part of the k-th partial
  --                derivative at x.

  procedure Multiply_Factor
              ( xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rcff,icff : in Standard_Floating_Vectors.Link_to_Vector;
                rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficient with the common factor.

  -- REQUIRED : facidx /= null.

  -- ON ENTRY :
  --   xpk          k-th exponent vector;
  --   facidx       factor index of k-th exponents in xpk;
  --   xr           real parts of the values for the variables;
  --   xi           imaginary parts of the values for the variables;
  --   rcff         real parts of the coefficient of the monomial;
  --   icff         imaginary parts coefficient of the monomial;
  --   rwrk         allocated space for the real parts
  --                of the coefficients of series of same degree;
  --   iwrk         allocated space for the imaginary high parts
  --                of the coefficients of series of same degree;
  --   racc         allocated space for the real parts 
  --                of coefficients of series of same degree;
  --   iacc         allocated space for the imaginary parts 
  --                of coefficients of series of same degree;
  --   rpwt         power table of the real parts for the values in x;
  --   ipwt         power table of the imaginary parts;
  --   u,v,w        are work space vectors of range 0..3.

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
  --   icff         imaginary parts of coefficients of a power series;

  -- ON RETURN :
  --   rcff         real parts of coefficients multiplied;
  --   icff         imaginary parts of coefficients multiplied.

  procedure Speel
              ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                rcff,icff : in Standard_Floating_VecVecs.VecVec;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rfwd,ifwd : in Standard_Floating_VecVecs.VecVec;
                rbck,ibck : in Standard_Floating_VecVecs.VecVec;
                rcrs,icrs : in Standard_Floating_VecVecs.VecVec;
                ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                rwrk,iwrk : in Standard_Floating_Vectors.Link_to_Vector;
                racc,iacc : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluation and differentiation of a polynomial,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   xps          exponent vector of the monomials;
  --   idx          indexed representation of the variables in exponents;
  --   fac          factor index of the exponents;
  --   rcff         real parts of coefficients of the products;
  --   icff         imaginary parts of coefficients of the products;
  --   xr           real parts of coefficients of series of same degree;
  --   xr           imaginary parts of coefficients of series of same degree;
  --   rfwd         work space allocated for xr'last-1 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   ifwd         work space allocated for xi'last-1 coefficient vectors
  --                of the same fixed degree as the series in xi;
  --   rbck         work space allocated for xr'last-2 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   ibck         work space allocated for xi'last-2 coefficient vectors
  --                of the same fixed degree as the series in xi;
  --   rcrs         work space allocated for xr'last-2 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   icrs         work space allocated for xi'last-2 coefficient vectors
  --                of the same fixed degree as the series in xr;
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
  --   u,v,w        are work space vectors of range 0..3.

  -- ON RETURN :
  --   ryd          ryd(xr'last+1) contains the real parts of the coefficient
  --                vector of the value of the sum of products evaluated at x,
  --                ryd(k) is the real part of the k-th partial derivative;
  --   iyd          ryd(xi'last+1) contains the real parts of the coefficient
  --                vector of the value of the sum of products evaluated at x,
  --                iyd(k) is the real part of the k-th partial derivative.

-- EVALUATION AND DIFFERENTIATION ON CIRCUITS :

  procedure EvalDiff
              ( c : in Circuit;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Wraps the Speel procedure for the convolution circuit c, to
  --   evaluate at the series with 4-vector representations,
  --   with the aid of the power table.

  -- ON ENTRY :
  --   c        a circuit properly defined and allocated;
  --   xr       real parts of coefficients of series of same degree;
  --   xi       imaginary parts of coefficients of series;
  --   rpwt     power table of the real parts for the values in x;
  --   ipwt     power table of the imaginary parts for the values in x.
  --   ryd      vector of range 0..rx'last+1 with space allocated for the
  --            coefficients of power series of the same fixed degree;
  --   iyd      vector of range 0..ix'last+1 with space allocated for the
  --            coefficients of power series of the same fixed degree;
  --   u,v,w    are work space vectors of range 0..3.

  -- ON RETURN :
  --   ryd      ryd(rx'last+1) contains the real parts of the
  --            coefficient vector of the value of the sum of products at x,
  --            ryd(k) is the real part of the k-th partial derivative;
  --   iyd      ryd(ihx'last+1) contains the imaginary parts of the
  --            coefficient vector of the value of the sum of products at x,
  --            iyd(k) is the imaginary high part of the k-th partial
  --            derivative;

  procedure EvalDiff
              ( c : in Circuits;
                xr,xi : in Standard_Floating_VecVecs.VecVec;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ryd,iyd : in Standard_Floating_VecVecs.VecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluates and differentiates the convolution circuits in c
  --   at the series x.

  -- ON ENTRY :
  --   c        an array of convolution circuits;
  --   xr       real parts of coefficients of series of same degree;
  --   xi       imaginary parts of coefficients of series;
  --   rpwt     power table of the real parts for the values in x;
  --   ipwt     power table of the imaginary parts for the values in x;
  --   ryd      work space of range 1..xr'last+1 to contain the real parts
  --            of the gradient and the value of the circuits in c;
  --   iyd      work space of range 1..xi'last+1 to contain the imaginary
  --            parts of the gradient and the value at x;
  --   vy       allocated space for the values of the circuits at x,
  --            done by the above procedure Linearized_Allocation,
  --   vm       space allocated for a series of some fixed degree
  --            with matrix coefficients;
  --   u,v,w    are work space vectors of range 0..3.

  -- ON RETURN :
  --   vy       values of the circuits at x, in linearized form;
  --   vm       the evaluated circuits at x as a series 
  --            of some fixe degree with matrix coefficients.

  procedure Delinearize ( vy,yv : in QuadDobl_Complex_VecVecs.VecVec );

  --  DESCRIPTION :
  --    Converts the linearized representation in vy into the vector yv
  --    of coefficient vectors of the series.
  --    This conversion is convenient for the difference computation
  --    and needed in the application of Newton's method.

  -- REQUIRED :
  --   if vy'range = 0..degree and vy(k)'range = 1..dimension,
  --   then yv'range = 1..dimension and yv(k)'range = 0..degree.

  procedure EvalDiff ( s : in System;
                       xr,xi : in Standard_Floating_VecVecs.VecVec;
                       u,v,w : in Standard_Floating_Vectors.Link_to_Vector );
  procedure EvalDiff ( s : in Link_to_System;
                       xr,xi : in Standard_Floating_VecVecs.VecVec;
                       u,v,w : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Wraps the EvalDiff on the convolution circuits in s.crc,
  --   at the power series with coefficients in x,
  --   represented by the vectors xr and xi.
  --   The work space vectors u, v, and w have range 0..3. 

  -- REQUIRED :
  --   All data in s are allocated properly with respect to dimension
  --   and degree, the power table is up to data with the given parts
  --   in the two vectors xr and xi.

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

end QuadDobl_Coefficient_Convolutions;
