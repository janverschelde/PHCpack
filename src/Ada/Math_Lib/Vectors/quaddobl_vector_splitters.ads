with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;

package QuadDobl_Vector_Splitters is

-- DESCRIPTION :
--   An 8-vector representation of a quad double complex vector consists
--   of 8 double vectors, with the highest, second highest, second lowest,
--   and lowest doubles of the real and imaginary parts.
--   A 2-vector representation of a quad double complex vector consists
--   of 2 double vectors, with real and imaginary parts stored in
--   sequences of four doubles, respectively the highest, second highest,
--   second lowest, and lowest double of each quad double number.
--   Convolutions run faster on the 2-vector representation.

-- 8-VECTOR REPRESENTATIONS :

  procedure Split ( x : in Complex_Number;
                    rehihi,imhihi,relohi,imlohi : out double_float;
                    rehilo,imhilo,relolo,imlolo : out double_float );

  -- DESCRIPTION :
  --   Splits the complex number in quad double precision into 8 parts.

  -- ON ENTRY :
  --   x        a quad double complex number.

  -- ON RETURN :
  --   rehihi   highest double of the real part of x,
  --            or high_part(high_part(real_part(x)));
  --   imhihi   highest double of the imaginary part of x;
  --            or high_part(high_part(imag_part(x)));
  --   relohi   second highest double of the real part of x;
  --            or low_part(high_part(real_part(x)));
  --   imlohi   second highest double of the imaginary part of x;
  --            or low_part(high_part(imag_part(x)));
  --   rehilo   second lowest double of the real part of x;
  --            or high_part(low_part(real_part(x)));
  --   imhilo   second lowest double of the imaginary part of x;
  --            or high_part(low_part(imag_part(x)));
  --   relolo   lowest double of the real part of x;
  --            or low_part(low_part(imag_part(x)));
  --   imlolo   lowest double of the imaginary part of x;
  --            or low_part(low_part(real_part(x))).

  procedure Merge ( x : out Complex_Number;
                    rehihi,imhihi,relohi,imlohi : in double_float;
                    rehilo,imhilo,relolo,imlolo : in double_float );

  -- DESCRIPTION :
  --   Merges the 8 doubles into a quad double precision complex number.

  -- ON ENTRY :
  --   rehihi   highest double of the real part for x;
  --   imhihi   highest double of the imaginary part for x;
  --   relohi   second highest double of the real part for x;
  --   imlohi   second highest double of the imaginary part for x;
  --   rehilo   second lowest double of the real part for x;
  --   imhilo   second lowest double of the imaginary part for x;
  --   relolo   lowest double of the real part for x;
  --   imlolo   lowest double of the imaginary part for x.

  -- ON RETURN :
  --   x        a quad double complex number, with
  --            rehihi = high_part(high_part(real_part(x))),
  --            imhihi = high_part(high_part(imag_part(x))),
  --            relohi = low_part(high_part(real_part(x))),
  --            imlohi = low_part(high_part(imag_part(x))),
  --            rehilo = high_part(low_part(real_part(x))),
  --            imhilo = high_part(low_part(imag_part(x))),
  --            relolo = low_part(low_part(real_part(x))),
  --            imlolo = low_part(low_part(imag_part(x))).

  procedure Split ( x : in QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : out Standard_Floating_Vectors.Vector;
                    xrlh,xilh : out Standard_Floating_Vectors.Vector;
                    xrhl,xihl : out Standard_Floating_Vectors.Vector;
                    xrll,xill : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the vector x of quad double precision complex numbers
  --   into 8 vectors with real and imaginary parts,
  --   and for each part into four doubles.

  -- REQUIRED :
  --   x'range = xrhh'range = xihh'range = xrlh'range = xilh'range
  --           = xrhl'range = xihl'range = xrll'range = xill'range.

  -- ON ENTRY :
  --   x        a vector of quad double precision complex numbers.

  -- ON RETURN :
  --   xrhh     highest doubles of the real parts of x;
  --   xihh     highest doubles of the imaginary parts of x;
  --   xrlh     second highest doubles of the real parts of x;
  --   xilh     second highest doubles of the imaginary parts of x;
  --   xrhl     second lowest doubles of the real parts of x;
  --   xihl     second lowest doubles of the imaginary parts of x;
  --   xrll     lowest doubles of the real parts of x;
  --   xill     lowest doubles of the imaginary parts of x.

  procedure Merge ( x : out QuadDobl_Complex_Vectors.Vector;
                    xrhh,xihh : in Standard_Floating_Vectors.Vector;
                    xrlh,xilh : in Standard_Floating_Vectors.Vector;
                    xrhl,xihl : in Standard_Floating_Vectors.Vector;
                    xrll,xill : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Merges the 8 vectors of doubles into one vector of
  --   quad double precision complex numbers.

  -- REQUIRED :
  --   x'range = xrhh'range = xihh'range = xrlh'range = xilh'range
  --           = xrhl'range = xihl'range = xrll'range = xill'range.

  -- ON ENTRY :
  --   xrhh     highest doubles of the real parts for x;
  --   xihh     highest doubles of the imaginary parts for x;
  --   xrlh     second highest doubles of the real parts for x;
  --   xilh     second highest doubles of the imaginary parts for x;
  --   xrhl     second lowest doubles of the real parts for x;
  --   xihl     second lowest doubles of the imaginary parts for x;
  --   xrll     lowest doubles of the real parts for x;
  --   xill     lowest doubles of the imaginary parts for x.

  -- ON RETURN :
  --   x        a vector of quad double precision complex numbers, with
  --            xrhh = high_part(high_part(real_part(x))),
  --            xihh = high_part(high_part(imag_part(x))),
  --            xrlh = low_part(high_part(real_part(x))),
  --            xilh = low_part(high_part(imag_part(x))),
  --            xrhl = high_part(low_part(real_part(x))),
  --            xihl = high_part(low_part(imag_part(x))),
  --            xrll = low_part(low_part(real_part(x))),
  --            xill = low_part(low_part(imag_part(x))).

-- 2-VECTOR REPRESENTATIONS :
--   The 2-vector representations are simpler to work with
--   and seem just as efficient as the 8-vector representations.

  procedure Two_Split
              ( x : in QuadDobl_Complex_Vectors.Vector;
                xr,xi : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Splits the vector of quad double complex numbers into
  --   two vectors of doubles, one with the real and the other
  --   with the imaginary parts.

  -- REQUIRED : x'first = xr'first = xi'first
  --   and xr'last = xi'last = 4*x'last.

  -- ON ENTRY :
  --   x        a vector of quad double complex numbers.

  -- ON RETURN :
  --   xr       the real parts of the complex numbers in x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the real part of x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the real part of x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the real part of x(k),
  --            x(4*(k-1) + 4) = lowest double of the real part of x(k);
  --   xi       the imaginary parts of the complex numbers in x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the imag part of x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the imag part of x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the imag part of x(k),
  --            x(4*(k-1) + 4) = lowest double of the imag part of x(k).

  procedure Two_Merge
              ( x : out QuadDobl_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Merges the 2 vectors of doubles into one vector of
  --   quad double precision complex numbers.

  -- REQUIRED : x'first = xr'first = xi'first
  --   and xr'last = xi'last = 4*x'last.

  -- ON ENTRY :
  --   xr       the real parts of the complex numbers for x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the real part for x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the real part for x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the real part for x(k),
  --            x(4*(k-1) + 4) = lowest double of the real part for x(k);
  --   xi       the imaginary parts of the complex numbers for x,
  --            the k-th number x(k) in x is stored as follows:
  --            x(4*(k-1) + 1) = highest double of the imag part for x(k),
  --            x(4*(k-1) + 2) = 2nd highest double of the imag part for x(k),
  --            x(4*(k-1) + 3) = 2nd lowest double of the imag part for x(k),
  --            x(4*(k-1) + 4) = lowest double of the imag part for x(k).

  -- ON RETURN :
  --   x        a vector of quad double complex numbers,
  --            with real parts taken xr and imaginary parts from xi.

-- SPLITTERS AND MERGERS WITH MEMORY ALLOCATIONS :

  procedure Split_Complex
              ( x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                xr,xi : out Standard_Floating_Vectors.Link_to_Vector );
  procedure Split_Complex
              ( x : in QuadDobl_Complex_VecVecs.VecVec;
                xr,xi : out Standard_Floating_VecVecs.VecVec );
  procedure Split_Complex
              ( x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                xr,xi : out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Splits the complex vector (of vectors) x into two real vectors,
  --   with its real and imaginary parts of the complex numbers in x.
  --   Memory is allocated for the resulting vectors.

-- MEMORY ALLOCATORS :

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return QuadDobl_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns allocated space for the complex coefficients
  --   of a series truncated to degree deg.

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec;
  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return QuadDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of a vector 
  --   of series, all truncated to degree deg.
  --   The vector on return has range 1..dim.

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return QuadDobl_Complex_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return QuadDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns an array of range neqstart..neq,
  --   with allocated vectors of range dimstart..dim.

-- SPLITTERS AND MERGERS ON ALLOCATED VECTORS :

  procedure Complex_Parts
              ( x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( x : in QuadDobl_Complex_VecVecs.VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Complex_Parts
              ( x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Parts the complex vector (of vectors) x into two vectors,
  --   with its real parts in xr and imaginary parts in xr,
  --   with highest, second highest, second lowest, and lowest doubles
  --   stored sequentially in xr and xi.

  -- REQUIRED :
  --   The vectors xr and xi are allocated to hold four times the
  --   numbers as in x, and x'first = xr'first = xi'first,
  --   for both zero or one.

  -- ON ENTRY :
  --   x        a complex vector of quad double precision.

  -- ON RETURN :
  --   xr       real parts of the quad double complex numbers in x;
  --   xi       imaginary parts of the quad double complex numbers in x.

  procedure Complex_Parts
              ( deg : in integer32;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Parts the complex vector (of vectors) x into two vectors,
  --   with its real parts in xr and imaginary parts in xr,
  --   with highest, second highest, second lowest, and lowest doubles
  --   stored sequentially in xr and xi.
  --   The coefficient vectors of the series are parted to degree deg.

  -- REQUIRED :
  --   The vectors xr and xi are allocated to hold four times the
  --   numbers as in x, and x'first = xr'first = xi'first,
  --   for both zero or one, at least up to degree deg.

  -- ON ENTRY :
  --   deg      degree of the series coefficients;
  --   x        a complex vector of quad double precision.

  -- ON RETURN :
  --   xr       real parts of the quad double complex numbers in x;
  --   xi       imaginary parts of the quad double complex numbers in x.

  procedure Complex_Merge
              ( xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure Complex_Merge
              ( xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Complex_Merge
              ( xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Merges the real and imaginary parts,
  --   into the vector (of vectors) of complex numbers.

  -- REQUIRED :
  --   The vector cvx is completely allocated
  --   and has compatible ragnes with xr and xi.

  -- ON ENTRY :
  --   xr      real parts of the numbers of a complex vector;
  --   xi      imaginary parts of the numbers of a complex vector;
  --   cvx     allocated with vectors with compatible ranges.

  -- ON RETURN :
  --   cvx     a complex vector of quad double precision, 
  --           with real parts from the doubles in xr,
  --           and imaginary parts from the doubles in xi.

  procedure Complex_Merge
              ( deg : in integer32;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure Complex_Merge
              ( deg : in integer32;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Complex_Merge
              ( deg : in integer32;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Merges the real and imaginary parts,
  --   into the vector (of vectors) of complex numbers.
  --   The coefficient vectors of the series are merged to degree deg.

  -- REQUIRED :
  --   The vector cvx is allocated, at least till the last index deg,
  --   and has compatible ranges with xr and xi.

  -- ON ENTRY :
  --   deg     degree of the series coefficient vectors;
  --   xr      real parts of the numbers of a complex vector;
  --   xi      imaginary parts of the numbers of a complex vector;
  --   cvx     allocated with vectors with compatible ranges.

  -- ON RETURN :
  --   cvx     a complex vector of quad double precision, 
  --           with real parts from the doubles in xr,
  --           and imaginary parts from the doubles in xi.

-- VECTOR COMPUTATIONS :

  procedure Add ( x,y : in Standard_Floating_Vectors.Vector;
                  z : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Adds two quad double numbers, x and y,
  --   and assigns the sum to the quad double number z.
  --   The vector representation for quad doubles is assumed.

  -- REQUIRED : x'range = y'range = z'range = 0..3.

  -- ON ENTRY :
  --   x        a quad double number, with
  --            x(0) the highest double of x,
  --            x(1) the second highest double of x,
  --            x(2) the second lowest double of x,
  --            x(3) the lowest double of x;
  --   y        a quad double number, with
  --            y(0) the highest double of y,
  --            y(1) the second highest double of y,
  --            y(2) the second lowest double of y,
  --            y(3) the lowest double of y.

  -- ON RETURN :
  --   z        the sum of x + y, with
  --            z(0) the highest double of x + y,
  --            z(1) the second highest double of x + y,
  --            z(2) the second lowest double of x + y,
  --            z(3) the lowest double of x + y.

  procedure Add ( x,y,z : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Add ( offset : in integer32;
                  x,y,z : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds two quad double numbers, x and y,
  --   and assigns the sum to the quad double number z.
  --   The vector representation for quad doubles is assumed.

  -- REQUIRED : x'range = y'range = z'range = 0..3.

  -- ON ENTRY :
  --   x        a quad double number, with
  --            x(offset) the highest double of x,
  --            x(offset+1) the second highest double of x,
  --            x(offset+2) the second lowest double of x,
  --            x(offset+3) the lowest double of x;
  --   y        a quad double number, with
  --            y(offset) the highest double of y,
  --            y(offset+1) the second highest double of y,
  --            y(offset+2) the second lowest double of y,
  --            y(offset+3) the lowest double of y.

  -- ON RETURN :
  --   z        the sum of x + y, with
  --            z(offset) the highest double of x + y,
  --            z(offset+1) the second highest double of x + y,
  --            z(offset+2) the second lowest double of x + y,
  --            z(offset+3) the lowest double of x + y.

  procedure Add ( zrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  zihh : in Standard_Floating_Vectors.Link_to_Vector;
                  zrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  zilh : in Standard_Floating_Vectors.Link_to_Vector;
                  zrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  zihl : in Standard_Floating_Vectors.Link_to_Vector;
                  zrll : in Standard_Floating_Vectors.Link_to_Vector;
                  zill : in Standard_Floating_Vectors.Link_to_Vector;
                  xrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  xihh : in Standard_Floating_Vectors.Link_to_Vector;
                  xrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  xilh : in Standard_Floating_Vectors.Link_to_Vector;
                  xrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  xihl : in Standard_Floating_Vectors.Link_to_Vector;
                  xrll : in Standard_Floating_Vectors.Link_to_Vector;
                  xill : in Standard_Floating_Vectors.Link_to_Vector;
                  yrhh : in Standard_Floating_Vectors.Link_to_Vector;
                  yihh : in Standard_Floating_Vectors.Link_to_Vector;
                  yrlh : in Standard_Floating_Vectors.Link_to_Vector;
                  yilh : in Standard_Floating_Vectors.Link_to_Vector;
                  yrhl : in Standard_Floating_Vectors.Link_to_Vector;
                  yihl : in Standard_Floating_Vectors.Link_to_Vector;
                  yrll : in Standard_Floating_Vectors.Link_to_Vector;
                  yill : in Standard_Floating_Vectors.Link_to_Vector;
                  x,y,z : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds two vectors of quad double complex numbers,
  --   in the 8-vector representation of doubles.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrhh     highest double of the real part of x;
  --   xihh     highest double of the imaginary part of x;
  --   xrlh     second highest double of the real part of x;
  --   xilh     second highest double of the imaginary part of x;
  --   xrhl     second lowest double of the real part of x;
  --   xihl     second lowest double of the imaginary part of x;
  --   xrll     lowest double of the real part of x;
  --   xill     lowest double of the imaginary part of x;
  --   yrhh     highest double of the real part of y;
  --   yihh     highest double of the imaginary part of y;
  --   yrlh     second highest double of the real part of y;
  --   yilh     second highest double of the imaginary part of y;
  --   yrhl     second lowest double of the real part of y;
  --   yihl     second lowest double of the imaginary part of y;
  --   yrll     lowest double of the real part of y;
  --   yill     lowest double of the imaginary part of y;
  --   x        a 4-vector work space of range 0..3;
  --   y        a 4-vector work space of range 0..3;
  --   z        a 4-vector work space of range 0..3.

  -- ON RETURN :
  --   zrhh     highest double of the real part of the sum;
  --   zihh     highest double of the imaginary part of the sum;
  --   zrlh     second highest double of the real part of the sum;
  --   zilh     second highest double of the imaginary part of the sum;
  --   zrhl     second lowest double of the real part of the sum;
  --   zihl     second lowest double of the imaginary part of the sum;
  --   zrll     lowest double of the real part of the sum;
  --   zill     lowest double of the imaginary part of the sum.

  procedure Two_Add 
              ( zr : in Standard_Floating_Vectors.Link_to_Vector;
                zi : in Standard_Floating_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                yr : in Standard_Floating_Vectors.Link_to_Vector;
                yi : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds two vectors of quad double complex numbers,
  --   in the 2-vector representation of doubles.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xr       real parts of the vector x;
  --   xi       imaginary parts of the vector x;
  --   yr       real parts of the vector y;
  --   yi       imaginary parts of the vector y.

  -- ON RETURN :
  --   zr       real parts of the sum;
  --   zi       imaginary parts of the sum.

  procedure Update
              ( offset : in integer32;
                z,y,x : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Updates the quad double starting at z(offset)
  --   with the quad double starting at y(offset),
  --   using the 4-vector x as work space.

  procedure Update
              ( zr : in Standard_Floating_Vectors.Link_to_Vector;
                zi : in Standard_Floating_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                wrk : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Updates z with x, in 2-vector representation.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   zr       current real parts of the vector z;
  --   zi       current imaginary parts of the vector z;
  --   xr       real parts of the vector x;
  --   xi       imaginary parts of the vector x;
  --   wrk      vector of range 0..3 for work space.

  -- ON RETURN :
  --   zr       updated real parts of the sum of x and z;
  --   zi       updated imaginary parts of the sum of x and z.

  procedure Update_Product
              ( zrhh,zihh,zrlh,zilh : in out double_float;
                zrhl,zihl,zrll,zill : in out double_float;
                xrhh,xihh,xrlh,xilh : in double_float;
                xrhl,xihl,xrll,xill : in double_float;
                yrhh,yihh,yrlh,yilh : in double_float;
                yrhl,yihl,yrll,yill : in double_float;
                x,y,z : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Adds the product of two quad double complex numbers to z,
  --   represented by 8 doubles.

  -- ON ENTRY :
  --   zrhh     current highest double in the real part of the inner product;
  --   zihh     current highest double in the imaginary part of the product;
  --   zrlh     current second highest double in the real part of the product;
  --   zilh     current second highest double in the imaginary part;
  --   zrlh     current second lowest double in the real part of the product;
  --   zilh     current second lowest double in the imaginary part;
  --   zrll     current lowest double in the real part of the inner product;
  --   zill     current lowest double in the imaginary part of the product;
  --   xrhh     highest double in the real part of the 1st number
  --   xihh     highest double in the imaginary part of the 1st number
  --   xrlh     second highest double in the real part of the 1st number
  --   xilh     second highest double in the imaginary part of the 1st number;
  --   xrlh     second lowest double in the real part of the 1st number
  --   xilh     second lowest double in the imaginary part of the 1st number;
  --   xrll     lowest double in the real part of the 1st number;
  --   xill     lowest double in the imaginary part of the 1st number;
  --   yrhh     highest double in the real part of the 2nd number;
  --   yihh     highest double in the imaginary part of the 2nd number;
  --   yrlh     second highest double in the real part of the 2nd number;
  --   yilh     second highest double in the imaginary part of the 2nd number;
  --   yrlh     second lowest double in the real part of the 2nd number
  --   yilh     second lowest double in the imaginary part of the 2nd number;
  --   yrll     lowest double in the real part of the 2nd number;
  --   yill     lowest double in the imaginary part of the 2nd number;
  --   x,y,z    vectors of range 0..3 as work space.

  -- ON RETURN :
  --   zrhh     updated highest double in the real part of the inner product;
  --   zihh     updated highest double in the imaginary part of the product;
  --   zrlh     updated second highest double in the real part of the product;
  --   zilh     updated second highest double in the imaginary part;
  --   zrlh     updated second lowest double in the real part of the product;
  --   zilh     updated second lowest double in the imaginary part;
  --   zrll     updated lowest double in the real part of the inner product;
  --   zill     updated lowest double in the imaginary part of the product.

  procedure Inner_Product
              ( zrhh,zihh,zrlh,zilh : out double_float;
                zrhl,zihl,zrll,zill : out double_float;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                yr,yi : in Standard_Floating_Vectors.Link_to_Vector;
                x,y,z : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes the inner product of two quad double complex vectors,
  --   stores as sequences of doubles.

  -- REQUIRED : xr'range = xi'range = yr'range = yi'range.

  -- ON ENTRY :
  --   xr       real parts of the first vector in the inner product;
  --   xi       imaginary parts of the first vector in the inner product;
  --   yr       real parts of the second vector in the inner product;
  --   yi       imaginary parts of the second vector in the inner product;
  --   x,y,z    work space of range 0..3.

  -- ON RETURN :
  --   zrhh     highest double in the real part of the inner product;
  --   zihh     highest double in the imaginary part of the product;
  --   zrlh     second highest double in the real part of the inner product;
  --   zilh     second highest double in the imaginary part of the product;
  --   zrlh     second lowest double in the real part of the inner product;
  --   zilh     second lowest double in the imaginary part of the product;
  --   zrll     lowest double in the real part of the inner product;
  --   zill     lowest double in the imaginary part of the product.

  procedure Multiply
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                yr : in Standard_Floating_Vectors.Link_to_Vector;
                yi : in Standard_Floating_Vectors.Link_to_Vector;
                zr : in Standard_Floating_Vectors.Link_to_Vector;
                zi : in Standard_Floating_Vectors.Link_to_Vector;
                x,y,z : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Multiplies the coefficients (xr, xi) of the first series
  --   with the coefficients (yr, yi) of the second series
  --   and stores the results in (zr, zi).

  -- REQUIRED : xr'first = 0 and
  --   xr'range = xi'range = yr'range = yi'range = zr'range = zi'range,
  --   and x'range = y'range = z'range = 0..3.

  -- ON ENTRY :
  --   xr       real parts of the first vector in the convolution;
  --   xi       imaginary parts of the first vector in the convolution;
  --   yr       real parts of the second vector in the convolution;
  --   yi       imaginary parts of the second vector in the convolution;
  --   x,y,z    work space vectors for quad double addition.

  -- ON RETURN :
  --   zr       real parts of the result of the convolution;
  --   zi       imaginary parts of the results of the convolution.

end QuadDobl_Vector_Splitters;
