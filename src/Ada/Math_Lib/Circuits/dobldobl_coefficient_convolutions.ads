with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;

package DoblDobl_Coefficient_Convolutions is

-- DESCRIPTION :
--   The convolutions on power series in double double precision work
--   on complex vectors in a 4-vector representation,
--   defined by the high parts of the real and imaginary parts, and 
--   the low parts of the real and imaginary parts, abbreviated 
--   respectively by real high (rh), imaginary high (ih), real low (rl),
--   and imaginary low (il).

-- DATA STRUCTURES :

-- A convolution circuit is a data structure for the efficient evaluation
-- and differentiation of polynomials in several variables at the
-- coefficient vectors of power series using the reverse mode of
-- algorithmic differentiation.

  type Circuit ( nbr,dim,dim1,dim2 : integer32 ) is record
    xps : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponent vectors
    idx : Standard_Integer_VecVecs.VecVec(1..nbr); -- exponent indices
    fac : Standard_Integer_VecVecs.VecVec(1..nbr); -- factor indices
   -- the complex coefficients of monomials are stored in 4 vectors
    rhcf : Standard_Floating_VecVecs.VecVec(1..nbr); -- real high parts
    ihcf : Standard_Floating_VecVecs.VecVec(1..nbr); -- imaginary high parts
    rlcf : Standard_Floating_VecVecs.VecVec(1..nbr); -- real low parts
    ilcf : Standard_Floating_VecVecs.VecVec(1..nbr); -- imaginary low parts
   -- the complex constant coefficients is stored in 4 vectors
    rhct : Standard_Floating_Vectors.Link_to_Vector; -- real high parts
    ihct : Standard_Floating_Vectors.Link_to_Vector; -- imaginary high parts
    rlct : Standard_Floating_Vectors.Link_to_Vector; -- real low parts
    ilct : Standard_Floating_Vectors.Link_to_Vector; -- imaginary low parts
   -- workspace for products and coefficient vectors of series are
   -- stored in real and imaginary parts, in high and low parts
   -- the forward products in 4-vector representation :
    rhfwd,ihfwd,rlfwd,ilfwd : Standard_floating_VecVecs.VecVec(1..dim1);
   -- the backward products in 4-vector representation :
    rhbck,ihbck,rlbck,ilbck : Standard_floating_VecVecs.VecVec(1..dim2);
   -- the cross products in 4-vector representation :
    rhcrs,ihcrs,rlcrs,ilcrs : Standard_floating_VecVecs.VecVec(1..dim2);
   -- real and imaginary parts of workspace and accumulator series
    rhwrk,ihwrk,rlwrk,ilwrk : Standard_Floating_Vectors.Link_to_Vector;
    rhacc,ihacc,rlacc,ilacc : Standard_Floating_Vectors.Link_to_Vector;
  end record;

  type Link_to_Circuit is access Circuit;

  type Circuits is array ( integer32 range <> ) of Link_to_Circuit;

  type Link_to_Circuits is access Circuits;

-- A system stores the sequence of convolution circuits for each polynomial
-- and the work space to hold auxiliary results and the final outcomes.

  type System ( neq,neq1,dim,dim1,deg : integer32 ) is record
    crc : Circuits(1..neq);    -- circuits for the equations
    mxe : Standard_Integer_Vectors.Vector(1..dim); -- exponent maxima
   -- the power table is stored in four parts
    rhpwt : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    ihpwt : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    rlpwt : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    ilpwt : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
   -- work space for EvalDiff on one circuit in ryd and iyd
    rhyd : Standard_Floating_VecVecs.VecVec(1..dim1); -- real high parts
    ihyd : Standard_Floating_VecVecs.VecVec(1..dim1); -- imaginary high parts
    rlyd : Standard_Floating_VecVecs.VecVec(1..dim1); -- real low parts
    ilyd : Standard_Floating_VecVecs.VecVec(1..dim1); -- imaginary low parts
   -- linearized evaluated power series in vy
    vy : DoblDobl_Complex_VecVecs.VecVec(0..deg);
   -- delinearized evaluated power series in yv
    yv : DoblDobl_Complex_VecVecs.VecVec(1..neq);
   -- differentiation result as matrix series in vm
    vm : DoblDobl_Complex_VecMats.VecMat(0..deg);
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
             return DoblDobl_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of the series
  --   truncated to degree deg and initialized to zero, in linearized form.
  --   The vector on return has range 0..deg and represents a series
  --   truncated to degree deg and with vectors of range 1..dim as
  --   its coefficients.

  function Allocate_Coefficients
             ( nbq,nvr,deg : integer32 )
             return DoblDobl_Complex_VecMats.VecMat;

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
              ( rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Computes the powers in the allocated power table for the
  --   coefficients in the power series in x.

  -- REQUIRED : rhpwt = Standard_Floating_VecVecVecs.Allocate(mxe,deg),
  --   ihpwt = Standard_Floating_VecVecVecs.Allocate(mxe,deg),
  --   rlpwt = Standard_Floating_VecVecVecs.Allocate(mxe,deg), and
  --   ilpwt = Standard_Floating_VecVecVecs.Allocate(mxe,deg)
  --   have been executed, and *pwt'range = x'range.

  -- ON ENTRY :
  --   rhpwt    allocated power table for the exponent maxima in mxe
  --            and for the real high parts of the coefficients of series,
  --            of the degree of the series in x;
  --   ihpwt    allocated power table for the exponent maxima in mxe
  --            and for the imaginary high parts of the coefficients of
  --            series, of the degree of the series in x;
  --   rlpwt    allocated power table for the exponent maxima in mxe
  --            and for the real low parts of the coefficients of series,
  --            of the degree of the series in x;
  --   ilpwt    allocated power table for the exponent maxima in mxe
  --            and for the imaginary low parts of the coefficients of
  --            series, of the degree of the series in x;
  --   mxe      exponent maxima for each variable in the system;
  --   rhx      real high parts of coefficients of power series;
  --   ihx      imaginary high parts of coefficients of power series;
  --   rlx      real low parts of coefficients of power series;
  --   ilx      imaginary low parts of coefficients of power series.

  -- ON RETURN :
  --   rhpwt    real high parts of the coefficients of the power table;
  --   ihpwt    imaginary high parts of the coefficients of the power table;
  --   rlpwt    real low parts of the coefficients of the power table;
  --   ilpwt    imaginary low parts of the coefficients of the power table.


-- REVERSE MODE OF ALGORITHMIC DIFFERENTIATION :

  procedure Speel ( rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                    rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                    rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                    rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                    rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                    rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                    rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                    idx : in Standard_Integer_Vectors.Vector;
                    rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                    rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                    rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                    rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                    rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                    rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Inline computation of the coefficients of the product of the
  --   series in rhx, ihx, rlx, ilx, which are all of the same degree.

  -- REQUIRED :
  --   The range of rx and ix starts at 1 and ends at 2 or higher.
  --   All vectors are allocated and of the range 0..degree,
  --   for some same fixed degree.

  -- ON ENTRY :
  --   rhx          a vector with range starting at 1, ending at 2 or higher,
  --                with the real high parts of the coefficients of series;
  --   ihx          a vector with range starting at 1, ending at 2 or higher,
  --                with the imaginary high parts of the coefficients;
  --   idx          if provided, then only those indices of x with values
  --                in idx participate and dim = idx'last,
  --                otherwise, all values of x participate and dim = x'last;
  --   rhfwd,ihfwd  space allocated for dim-1 coefficient vectors,
  --                for high parts of series of the same fixed degree;
  --   rlfwd,ilfwd  space allocated for dim-1 coefficient vectors,
  --                for low parts of series of the same fixed degree;
  --   rhbck,ihbck  has space reserved for dim-2 coefficient vectors,
  --                for high parts series of the same fixed degree;
  --   rlbck,ilbck  has space reserved for dim-2 coefficient vectors,
  --                for low parts series of the same fixed degree;
  --   rhcrs,ihcrs  has space reserved for dim-2 coefficient vectors,
  --                for high parts series of the same fixed degree;
  --   rlcrs,ilcrs  has space reserved for dim-2 coefficient vectors,
  --                for low parts series of the same fixed degree.

  -- ON RETURN :
  --   rhfwd        accumulates the real high parts of the forward products,
  --                rhfwd(dim-1) holds the coefficients for the product,
  --                rhfwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   ihfwd        accumulates the imaginary high parts of forward products,
  --                ihfwd(dim-1) holds the coefficients for the product,
  --                ihfwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   rlfwd        accumulates the real low parts of the forward products,
  --                rlfwd(dim-1) holds the coefficients for the product,
  --                rlfwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   ilfwd        accumulates the imaginary low parts of forward products,
  --                ilfwd(dim-1) holds the coefficients for the product,
  --                ilfwd(dim-2) holds the coefficients for the value
  --                of the last partial derivative in x;
  --   rhbck        accumulates the real high parts of the backward products,
  --                rhbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   ihbck        accumulates the imaginary high parts of backward products,
  --                rhbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   rlbck        accumulates the real low parts of the backward products,
  --                rlbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   ilbck        accumulates the imaginary low parts of backward products,
  --                rlbck(dim-2) holds the coefficients of the first
  --                partial derivative of the product;
  --   rhcrs        stores the real high parts of the cross products,
  --                rhcrs(k) contains the real parts of the coefficients
  --                of the partial derivative w.r.t. k+1;
  --   ihcrs        stores the imaginary high parts of the cross products,
  --                ihcrs(k) contains the imaginary parts of coefficients
  --                of the partial derivative w.r.t. k+1;
  --   rlcrs        stores the real low parts of the cross products,
  --                rhcrs(k) contains the real parts of the coefficients
  --                of the partial derivative w.r.t. k+1;
  --   ilcrs        stores the imaginary low parts of the cross products,
  --                ihcrs(k) contains the imaginary parts of coefficients
  --                of the partial derivative w.r.t. k+1.

  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                    rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                    rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                    rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                    rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                    rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                    rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec;
                    rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                    rlyd,ilyd : in Standard_Floating_VecVecs.VecVec );
  procedure Speel ( idx : in Standard_Integer_VecVecs.VecVec;
                    rhcff,ihcff : in Standard_Floating_VecVecs.VecVec;
                    rlcff,ilcff : in Standard_Floating_VecVecs.VecVec;
                    rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                    rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                    rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                    rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                    rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                    rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                    rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec;
                    rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                    rlyd,ilyd : in Standard_Floating_VecVecs.VecVec;
                    rhwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    rlwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    ihwrk : in Standard_Floating_Vectors.Link_to_Vector;
                    ilwrk : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Evaluation and differentiation of the sum of products,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   idx          indexed representation of a sum of products of variables;
  --   rhcff        real high parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   ihcff        imaginary high parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   rlcff        real low parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   ilcff        imaginary low parts of coefficients of the products,
  --                if omitted, then all coefficients are considered as one;
  --   rhx          real high parts of coefficients of series of same degree;
  --   ihx          imaginary high parts of coefficients of series;
  --   rlx          real low parts of coefficients of series of same degree;
  --   ilx          imaginary low parts of coefficients of series;
  --   rhfwd        work space allocated for rhx'last-1 coefficient vectors
  --                of the same fixed degree as the series in rhx;
  --   ihfwd        work space allocated for ihx'last-1 coefficient vectors
  --                of the same fixed degree as the series in ihx;
  --   rlfwd        work space allocated for rlx'last-1 coefficient vectors
  --                of the same fixed degree as the series in rlx;
  --   ilfwd        work space allocated for ilx'last-1 coefficient vectors
  --                of the same fixed degree as the series in ilx;
  --   rhbck        work space allocated for rhx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rhx;
  --   ihbck        work space allocated for ihx'last-2 coefficient vectors
  --                of the same fixed degree as the series in ihx;
  --   rlbck        work space allocated for rlx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rlx;
  --   ilbck        work space allocated for ilx'last-2 coefficient vectors
  --                of the same fixed degree as the series in ilx;
  --   rhcrs        work space allocated for rhx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rhx;
  --   ihcrs        work space allocated for ihx'last-2 coefficient vectors
  --                of the same fixed degree as the series in ihx;
  --   rlcrs        work space allocated for rlx'last-2 coefficient vectors
  --                of the same fixed degree as the series in rlx;
  --   ilcrs        work space allocated for ilx'last-2 coefficient vectors
  --                of the same fixed degree as the series in ilx;
  --   rhyd         vector of range 0..rhx'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   ihyd         vector of range 0..ihx'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rlyd         vector of range 0..rlx'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   ilyd         vector of range 0..ilx'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rhwrk        work space for the coefficients of the same fixed degree;
  --   ihwrk        work space for the coefficients of the same fixed degree;
  --   rlwrk        work space for the coefficients of the same fixed degree;
  --   ilwrk        work space for the coefficients of the same fixed degree.

  -- ON RETURN :
  --   rhyd         rhyd(rhx'last+1) contains the real high parts of
  --                coefficients of the value of the sum of products at x;
  --                rhyd(k) is the real high part of the k-th partial
  --                derivative at x;
  --   ihyd         ihyd(ihx'last+1) contains the imaginary high parts of
  --                coefficients of the value of the sum of products at x,
  --                ihyd(k) is the imaginary high part of the k-th partial
  --                derivative at x.
  --   rlyd         rlyd(rlx'last+1) contains the real low parts of
  --                coefficients of the value of the sum of products at x;
  --                rlyd(k) is the real low part of the k-th partial
  --                derivative at x;
  --   ilyd         ilyd(ilx'last+1) contains the imaginary low parts of
  --                coefficients of the value of the sum of products at x,
  --                ilyd(k) is the imaginary low part of the k-th partial
  --                derivative at x.

  procedure Multiply_Factor
              ( xpk,facidx : in Standard_Integer_Vectors.Link_to_Vector;
                rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhcff,ihcff : in Standard_Floating_Vectors.Link_to_Vector;
                rlcff,ilcff : in Standard_Floating_Vectors.Link_to_Vector;
                rhwrk,ihwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rlwrk,ilwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rhacc,ihacc : in Standard_Floating_Vectors.Link_to_Vector;
                rlacc,ilacc : in Standard_Floating_Vectors.Link_to_Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Multiplies the coefficient with the common factor.

  -- REQUIRED : facidx /= null.

  -- ON ENTRY :
  --   xpk          k-th exponent vector;
  --   facidx       factor index of k-th exponents in xpk;
  --   rhx          real high parts of the values for the variables;
  --   ihx          imaginary high parts of the values for the variables;
  --   rlx          real low parts of the values for the variables;
  --   ilx          imaginary low parts of the values for the variables;
  --   rhcff        real high parts of the coefficient of the monomial;
  --   ihcff        imaginary high parts coefficient of the monomial;
  --   rlcff        real low parts of the coefficient of the monomial;
  --   ilcff        imaginary low parts coefficient of the monomial;
  --   rhwrk        allocated space for the real high parts
  --                of the coefficients of series of same degree;
  --   ihwrk        allocated space for the imaginary high parts
  --                of the coefficients of series of same degree;
  --   rlwrk        allocated space for the real low parts
  --                of the coefficients of series of same degree;
  --   ilwrk        allocated space for the imaginary low parts
  --                of the coefficients of series of same degree;
  --   rhacc        allocated space for the real high parts 
  --                of coefficients of series of same degree;
  --   ihacc        allocated space for the imaginary high parts 
  --                of coefficients of series of same degree;
  --   rlacc        allocated space for the real low parts 
  --                of coefficients of series of same degree;
  --   ilacc        allocated space for the imaginary low parts 
  --                of coefficients of series of same degree;
  --   rhpwt        power table of the real high parts for the values in x;
  --   ihpwt        power table of the imaginary high parts;
  --   rlpwt        power table of the real lowparts for the values in x;
  --   ilpwt        power table of the imaginary low parts.

  -- ON RETURN :
  --   rhacc        accumulates the product of the real high parts
  --                of the coefficients with the evaluated powers of x
  --                as defined by xpk and the factor index in facidx.
  --   ihacc        accumulates the product of the imaginary high parts
  --                of the coefficients with the evaluated powers of x
  --                as defined by xpk and the factor index in facidx;
  --   rlacc        accumulates the product of the real low parts
  --                of the coefficients with the evaluated powers of x
  --                as defined by xpk and the factor index in facidx.
  --   ilacc        accumulates the product of the imaginary low parts
  --                of the coefficients with the evaluated powers of x
  --                as defined by xpk and the factor index in facidx.

  procedure Multiply_Power
              ( multiplier : in integer32;
                rhcff : in Standard_Floating_Vectors.Link_to_Vector; 
                ihcff : in Standard_Floating_Vectors.Link_to_Vector; 
                rlcff : in Standard_Floating_Vectors.Link_to_Vector; 
                ilcff : in Standard_Floating_Vectors.Link_to_Vector ); 

  -- DESCRIPTION :
  --   Multiplies the coefficients of the power series with multiplier.

  -- ON ENTRY:
  --   multiplier   is the multiplier exponent;
  --   rhcff        real high parts of coefficients of a power series;
  --   ihcff        imaginary high parts of coefficients of a power series;
  --   rlcff        real low parts of coefficients of a power series;
  --   ilcff        imaginary low parts of coefficients of a power series.

  -- ON RETURN :
  --   rhcff        real high parts of coefficients multiplied;
  --   ihcff        imaginary high parts of coefficients multiplied;
  --   rhcff        real low parts of coefficients multiplied;
  --   ihcff        imaginary low parts of coefficients multiplied.

  procedure Speel
              ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                rhcff,ihcff : in Standard_Floating_VecVecs.VecVec;
                rlcff,ilcff : in Standard_Floating_VecVecs.VecVec;
                rhx,ihx,rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhfwd,ihfwd : in Standard_Floating_VecVecs.VecVec;
                rlfwd,ilfwd : in Standard_Floating_VecVecs.VecVec;
                rhbck,ihbck : in Standard_Floating_VecVecs.VecVec;
                rlbck,ilbck : in Standard_Floating_VecVecs.VecVec;
                rhcrs,ihcrs : in Standard_Floating_VecVecs.VecVec;
                rlcrs,ilcrs : in Standard_Floating_VecVecs.VecVec;
                rhyd,ihyd,rlyd,ilyd : in Standard_Floating_VecVecs.VecVec;
                rhwrk,ihwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rlwrk,ilwrk : in Standard_Floating_Vectors.Link_to_Vector;
                rhacc,ihacc : in Standard_Floating_Vectors.Link_to_Vector;
                rlacc,ilacc : in Standard_Floating_Vectors.Link_to_Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Evaluation and differentiation of a polynomial,
  --   given in indexed format at a power series.

  -- ON ENTRY :
  --   xps          exponent vector of the monomials;
  --   idx          indexed representation of the variables in exponents;
  --   fac          factor index of the exponents;
  --   rhcff        real high parts of coefficients of the products;
  --   ihcff        imaginary high parts of coefficients of the products;
  --   rlcff        real low parts of coefficients of the products;
  --   ilcff        imaginary low parts of coefficients of the products;
  --   rhx          real high parts of coefficients of series of same degree;
  --   ihx          imaginary high parts of coefficients of series;
  --   rlx          real low parts of coefficients of series of same degree;
  --   ilx          imaginary low parts of coefficients of series;
  --   rhfwd,ihfwd  work space allocated for xr'last-1 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   rlfwd,ilfwd  work space allocated for xi'last-1 coefficient vectors
  --                of the same fixed degree as the series in xi;
  --   rhbck,ihbck  work space allocated for xr'last-2 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   rlbck,ilbck  work space allocated for xi'last-2 coefficient vectors
  --                of the same fixed degree as the series in xi;
  --   rhcrs,ihcrs  work space allocated for xr'last-2 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   ilcrs,ilcrs  work space allocated for xi'last-2 coefficient vectors
  --                of the same fixed degree as the series in xr;
  --   rhyd,ihyd    vector of range 0..xr'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rlyd,ilyd    vector of range 0..xi'last with space allocated for the
  --                coefficients of power series of the same fixed degree;
  --   rhwrk,rlwrk  work space for the coefficients of the same fixed degree;
  --   ilwrk,ilwrk  work space for the coefficients of the same fixed degree;
  --   rhacc,rlacc  work space for the coefficients of the same fixed degree;
  --   ihacc,ilacc  work space for the coefficients of the same fixed degree;
  --   rhpwt,ihpwt  power table of the real parts for the values in x;
  --   rlpwt,ilpwt  power table of the real imaginary for the values in x.

  -- ON RETURN :
  --   rhyd         rhyd(rhx'last+1) contains the real high parts of
  --                the coefficients of the sum of products evaluated at x,
  --                rhyd(k) is the real part of the k-th partial derivative;
  --   ihyd         ihyd(ihx'last+1) contains the imaginary high parts of
  --                the coefficients of the sum of products evaluated at x,
  --                ihyd(k) is the imag part of the k-th partial derivative;
  --   rlyd         rlyd(rlx'last+1) contains the real low parts of
  --                the coefficients of the sum of products evaluated at x,
  --                rlyd(k) is the real part of the k-th partial derivative;
  --   ilyd         ilyd(ilx'last+1) contains the imaginary high parts of
  --                the coefficients of the sum of products evaluated at x,
  --                ilyd(k) is the imag part of the k-th partial derivative.

-- EVALUATION AND DIFFERENTIATION ON CIRCUITS :

  procedure EvalDiff
              ( c : in Circuit;
                rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                rlyd,ilyd : in Standard_Floating_VecVecs.VecVec );


  -- DESCRIPTION :
  --   Wraps the Speel procedure for the convolution circuit c, to
  --   evaluate at the series with 4-vector representations,
  --   with the aid of the power table.

  -- ON ENTRY :
  --   c        a circuit properly defined and allocated;
  --   rhx      real high parts of coefficients of series of same degree;
  --   ihx      imaginary high parts of coefficients of series;
  --   rlx      real low parts of coefficients of series of same degree;
  --   ilx      imaginary low parts of coefficients of series;
  --   ryd      vector of range 0..rx'last+1 with space allocated for the
  --            coefficients of power series of the same fixed degree;
  --   iyd      vector of range 0..ix'last+1 with space allocated for the
  --            coefficients of power series of the same fixed degree;
  --   rhpwt    power table of the real high parts for the values in x;
  --   ihpwt    power table of the imaginary high for the values in x;
  --   rlpwt    power table of the real low parts for the values in x;
  --   ilpwt    power table of the imaginary low for the values in x.

  -- ON RETURN :
  --   rhyd     rhyd(rhx'last+1) contains the real high parts of the
  --            coefficient vector of the value of the sum of products at x,
  --            rhyd(k) is the real high part of the k-th partial derivative;
  --   ihyd     rhyd(ihx'last+1) contains the imaginary high parts of the
  --            coefficient vector of the value of the sum of products at x,
  --            ihyd(k) is the imaginary high part of the k-th partial
  --            derivative;
  --   rlyd     rlyd(rlx'last+1) contains the real high parts of the coefficient
  --            coefficient vector of the value of the sum of products at x,
  --            rlyd(k) is the real high part of the k-th partial derivative;
  --   ilyd     rlyd(ilx'last+1) contains the imaginary high parts of the
  --            coefficient vector of the value of the sum of products at x,
  --            ihyd(k) is the imaginary high part of the k-th partial
  --            derivative.

  procedure EvalDiff
              ( c : in Circuits;
                rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.VecVec;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rhyd,ihyd : in Standard_Floating_VecVecs.VecVec;
                rlyd,ilyd : in Standard_Floating_VecVecs.VecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat );

  -- DESCRIPTION :
  --   Evaluates and differentiates the convolution circuits in c
  --   at the series x.

  -- ON ENTRY :
  --   c        an array of convolution circuits;
  --   rhx      real high parts of coefficients of series of same degree;
  --   ihx      imaginary high parts of coefficients of series;
  --   rlx      real low parts of coefficients of series of same degree;
  --   ilx      imaginary low parts of coefficients of series;
  --   ryd      vector of range 0..rx'last+1 with space allocated for the
  --            coefficients of power series of the same fixed degree;
  --   iyd      vector of range 0..ix'last+1 with space allocated for the
  --            coefficients of power series of the same fixed degree;
  --   rhpwt    power table of the real high parts for the values in x;
  --   ihpwt    power table of the imaginary high for the values in x;
  --   rlpwt    power table of the real low parts for the values in x;
  --   ilpwt    power table of the imaginary low for the values in x.
  --   rhyd     work space of range 1..rx'last+1 to contain the real high
  --            parts of the gradient and the value of the circuits in c;
  --   ihyd     work space of range 1..ix'last+1 to contain the imaginary
  --            high parts of the gradient and the value at x;
  --   rlyd     work space of range 1..rx'last+1 to contain the real low
  --            parts of the gradient and the value of the circuits in c;
  --   ilyd     work space of range 1..ix'last+1 to contain the imaginary
  --            low parts of the gradient and the value at x;
  --   vy       allocated space for the values of the circuits at x,
  --            done by the above procedure Linearized_Allocation,
  --   vm       space allocated for a series of some fixed degree
  --            with matrix coefficients.

  -- ON RETURN :
  --   vy       values of the circuits at x, in linearized form;
  --   vm       the evaluated circuits at x as a series 
  --            of some fixe degree with matrix coefficients.

  procedure Delinearize ( vy,yv : in DoblDobl_Complex_VecVecs.VecVec );

  --  DESCRIPTION :
  --    Converts the linearized representation in vy into the vector yv
  --    of coefficient vectors of the series.
  --    This conversion is convenient for the difference computation
  --    and needed in the application of Newton's method.

  -- REQUIRED :
  --   if vy'range = 0..degree and vy(k)'range = 1..dimension,
  --   then yv'range = 1..dimension and yv(k)'range = 0..degree.

  procedure EvalDiff ( s : in System;
                       rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                       rlx,ilx : in Standard_Floating_VecVecs.VecVec );
  procedure EvalDiff ( s : in Link_to_System;
                       rhx,ihx : in Standard_Floating_VecVecs.VecVec;
                       rlx,ilx : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Wraps the EvalDiff on the convolution circuits in s.crc,
  --   at the power series with coefficients in x,
  --   represented by the vectors rhx, ihx, rlx, and ilx.

  -- REQUIRED :
  --   All data in s are allocated properly with respect to dimension
  --   and degree, the power table is up to data with the given parts
  --   in the four vectors.

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

end DoblDobl_Coefficient_Convolutions;
