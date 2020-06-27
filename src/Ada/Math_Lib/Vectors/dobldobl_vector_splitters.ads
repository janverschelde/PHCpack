with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;

package DoblDobl_Vector_Splitters is

-- DESCRIPTION :
--   A 4-vector representation of a double double complex vector consists
--   of 4 double float vectors, respectively with the high parts of the
--   real and imaginary parts, and the low parts of the real and imaginary
--   parts of the complex numbers in the double double complex vector.
--   Several computations run faster on a 4-vector representation than
--   on a vector of double double complex numbers.

  procedure Split ( x : in Complex_Number;
                    rehi,imhi,relo,imlo : out double_float );

  -- DESCRIPTION :
  --   Splits the complex number in real and imaginary parts,
  --   in high and low parts.

  -- ON ENTRY :
  --   x        a double double complex number.

  -- ON RETURN :
  --   rehi     high part of the real part of x;
  --   imhi     high part of the imaginary part of x;
  --   relo     low part of the real part of x;
  --   imlo     low part of the imaginary part of x.

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    rehi : out Standard_Floating_Vectors.Vector;
                    imhi : out Standard_Floating_Vectors.Vector;
                    relo : out Standard_Floating_Vectors.Vector;
                    imlo : out Standard_Floating_Vectors.Vector );
  
  -- DESCRIPTION :
  --   Splits the vector v into real and imaginary,
  --   into high and low parts.

  -- REQUIRED : v'range = rehi'range = imhi'range = relo'range = imlo'range.

  -- ON ENTRY :
  --   v        a vector of double double complex numbers.

  -- ON RETURN :
  --   rehi     high parts of the real parts of the numbers in v;
  --   imhi     high parts of the imaginary parts of the numbers in v;
  --   relo     low parts of the real parts of the numbers in v;
  --   imlo     low parts of the imaginary parts of the numbers in v.

  procedure Split_Complex
              ( x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                rhpx,ihpx : out Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : out Standard_Floating_Vectors.Link_to_Vector );
  procedure Split_Complex
              ( x : in DoblDobl_Complex_VecVecs.VecVec;
                rhpx,ihpx : out Standard_Floating_VecVecs.VecVec;
                rlpx,ilpx : out Standard_Floating_VecVecs.VecVec );
  procedure Split_Complex
              ( x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                rhpx,ihpx : out Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : out Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Splits the complex vector (of vectors) x into two real vectors,
  --   with its real and imaginary parts of the complex numbers in x.
  --   Memory is allocated for the resulting vectors.

-- MEMORY ALLOCATORS :

  function Allocate_Complex_Coefficients
             ( deg : integer32 )
             return DoblDobl_Complex_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns allocated space for the complex coefficients
  --   of a series truncated to degree deg.

  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate_Complex_Coefficients
             ( dim,deg : integer32 )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns allocated space for the coefficients of a vector 
  --   of series, all truncated to degree deg.
  --   The vector on return has range 1..dim.

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return DoblDobl_Complex_VecVecs.VecVec;
  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return DoblDobl_Complex_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns an array of range neqstart..neq,
  --   with allocated vectors of range dimstart..dim.

-- MERGE PROCEDURES :

  procedure Merge ( x : out Complex_Number;
                    rehi,imhi,relo,imlo : in double_float );

  -- DESCRIPTION :
  --   Merges the high and low parts of real and imaginary parts
  --   into one complex number.

  -- ON ENTRY :
  --   rehi     high part of the real part;
  --   imhi     high part of the imaginary part;
  --   relo     low part of the real part;
  --   imlo     low part of the imaginary part.

  -- ON RETURN :
  --   x        complex number with high part of real part in rehi,
  --            high part of imaginary part in imhi,
  --            low part of real part in imlo,
  --            low part of imaginary part in imlo.

  procedure Merge ( v : out DoblDobl_Complex_Vectors.Vector;
                    rehi : in Standard_Floating_Vectors.Vector;
                    imhi : in Standard_Floating_Vectors.Vector;
                    relo : in Standard_Floating_Vectors.Vector;
                    imlo : in Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --   Merges the high and low parts of real and imaginary parts
  --   into one complex vector v.

  -- REQUIRED : v'range = rehi'range = imhi'range = relo'range = imlo'range.

  -- ON ENTRY :
  --   rehi     high parts of the real parts for the numbers in v;
  --   imhi     high parts of the imaginary parts or the numbers in v;
  --   relo     low parts of the real parts for the numbers in v;
  --   imlo     low parts of the imaginary parts for the numbers in v.

  -- ON RETURN :
  --   v        a vector of double double complex numbers, with
  --            high part of real parts in rehi,
  --            high part of imaginary parts in imhi,
  --            low part of real parts in relo, and
  --            low part of imaginary parts in imlo,

-- PROCEDURES TO PART AND MERGE VECTORS OF VECTORS :

  procedure Complex_Parts
              ( x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( x : in DoblDobl_Complex_VecVecs.VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Complex_Parts
              ( x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Parts the complex vector (of vectors) x into two four vectors,
  --   with its real and imaginary parts of the complex numbers in x,
  --   separated into real high, imaginary high, real low,
  --   and imaginary low.

  -- REQUIRED :
  --   The vectors rpx and ipx are completely allocated
  --   and have the same ranges as x.

  -- ON ENTRY :
  --   x        a complex vector of double double precision.

  -- ON RETURN :
  --   rhpx     high parts of the real numbers in x;
  --   ihpx     high parts of the imaginary numbers in x;
  --   rlpx     low parts of the real numbers in x;
  --   ilpx     low parts of the imaginary numbers in x.

  procedure Complex_Parts
              ( deg : in integer32;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec );
  procedure Complex_Parts
              ( deg : in integer32;
                x : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Parts the complex vector (of vectors) x into two four vectors,
  --   with its real and imaginary parts of the complex numbers in x,
  --   separated into real high, imaginary high, real low,
  --   and imaginary low.
  --   The coefficient vectors of the series are parted only up
  --   to the given degree deg.

  -- REQUIRED :
  --   The vectors rpx and ipx are allocated, at least till degree deg,
  --   and have the same ranges as x, at least up to the index deg.

  -- ON ENTRY :
  --   deg      degree of the series coefficients;
  --   x        a complex vector of double double precision.

  -- ON RETURN :
  --   rhpx     high parts of the real numbers in x;
  --   ihpx     high parts of the imaginary numbers in x;
  --   rlpx     low parts of the real numbers in x;
  --   ilpx     low parts of the imaginary numbers in x.

  procedure Complex_Merge
              ( rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure Complex_Merge
              ( rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Complex_Merge
              ( rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Merges the real and imaginary parts, high and low,
  --   into the vector (of vectors) of complex numbers.

  -- REQUIRED :
  --   The vector cvx is completely allocated
  --   and has the same ranges as rhpx, ihpx, rlpx, and ilpx.

  -- ON ENTRY :
  --   rhpx    high parts of real numbers of a complex vector;
  --   ihpx    high parts of imaginary numbers of a complex vector;
  --   rlpx    low parts of real numbers of a complex vector;
  --   ilpx    low parts of imaginary numbers of a complex vector;
  --   cvx     allocated with vectors of the same ranges.

  -- ON RETURN :
  --   cvx     a complex vector of double double precision, 
  --           with real high from rhpx, imaginary high from ihpx,
  --           real low from rlpx, and imaginary low from ilpx.

  procedure Complex_Merge
              ( deg : in integer32;
                rhpx,ihpx : in Standard_Floating_Vectors.Link_to_Vector;
                rlpx,ilpx : in Standard_Floating_Vectors.Link_to_Vector;
                cvx : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure Complex_Merge
              ( deg : in integer32;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.Link_to_VecVec );
  procedure Complex_Merge
              ( deg : in integer32;
                rhpx,ihpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlpx,ilpx : in Standard_Floating_VecVecs.Link_to_VecVec;
                cvx : in DoblDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Merges the real and imaginary parts, high and low,
  --   into the vector (of vectors) of complex numbers.

  -- REQUIRED :
  --   The vector cvx is allocated, at least up to deg,
  --   and has the same ranges as rhpx, ihpx, rlpx, and ilpx.

  -- ON ENTRY :
  --   deg     degree of the series coefficients;
  --   rhpx    high parts of real numbers of a complex vector;
  --   ihpx    high parts of imaginary numbers of a complex vector;
  --   rlpx    low parts of real numbers of a complex vector;
  --   ilpx    low parts of imaginary numbers of a complex vector;
  --   cvx     allocated with vectors of the same ranges.

  -- ON RETURN :
  --   cvx     a complex vector of double double precision, 
  --           with real high from rhpx, imaginary high from ihpx,
  --           real low from rlpx, and imaginary low from ilpx.

  procedure Add ( zrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  zimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  zrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  zimlo : in Standard_Floating_Vectors.Link_to_Vector;
                  xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                  xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                  yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  yimlo : in Standard_Floating_Vectors.Link_to_Vector);

  -- DESCRIPTION :
  --   Adds two double double complex vectors x and y to form
  --   the result z, in 4-vector representation.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrehi    high parts of the real parts for the numbers in x;
  --   ximhi    high parts of the imaginary parts for the numbers in x;
  --   xrelo    low parts of the real parts for the numbers in x;
  --   ximlo    low parts of the imaginary parts for the numbers in x;
  --   yrehi    high parts of the real parts for the numbers in y;
  --   yimhi    high parts of the imaginary parts or the numbers in y;
  --   yrelo    low parts of the real parts for the numbers in y;
  --   yimlo    low parts of the imaginary parts for the numbers in y.

  -- ON RETURN :
  --   zrehi    high parts of the real parts for the numbers in z;
  --   zimhi    high parts of the imaginary parts for the numbers in z;
  --   zrelo    low parts of the real parts for the numbers in z;
  --   zimlo    low parts of the imaginary parts for the numbers in z.

  procedure Update ( zrehi : in Standard_Floating_Vectors.Link_to_Vector;
                     zimhi : in Standard_Floating_Vectors.Link_to_Vector;
                     zrelo : in Standard_Floating_Vectors.Link_to_Vector;
                     zimlo : in Standard_Floating_Vectors.Link_to_Vector;
                     xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                     ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                     xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                     ximlo : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Updates z with x, in 4-vector representation.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrehi    high parts of the real parts for the numbers in x;
  --   ximhi    high parts of the imaginary parts for the numbers in x;
  --   xrelo    low parts of the real parts for the numbers in x;
  --   ximlo    low parts of the imaginary parts for the numbers in x.

  -- ON RETURN :
  --   zrehi    high parts of the real parts for the numbers in z;
  --   zimhi    high parts of the imaginary parts for the numbers in z;
  --   zrelo    low parts of the real parts for the numbers in z;
  --   zimlo    low parts of the imaginary parts for the numbers in z.

  procedure Update_Product
                ( zrehi,zimhi,zrelo,zimlo : in out double_float;
                  xrehi,ximhi,xrelo,ximlo : in double_float;
                  yrehi,yimhi,yrelo,yimlo : in double_float );

  -- DESCRIPTION :
  --   Updates the double double complex number z,
  --   given in 4-number representation,
  --   with the product of two double double complex numbers x and y,
  --   also given in their 4-number representation.

  -- ON ENTRY :
  --   xrehi    high part of the real part for the number x;
  --   ximhi    high part of the imaginary part for the number x;
  --   xrelo    low part of the real part for the number in x;
  --   ximlo    low part of the imaginary part for the number in x;
  --   yrehi    high part of the real part for the number in y;
  --   yimhi    high part of the imaginary part for the number in y;
  --   yrelo    low part of the real part for the number in y;
  --   yimlo    low part of the imaginary part for the number in y;
  --   zrehi    high part of the real part for the number in z;
  --   zimhi    high part of the imaginary part for the number in z;
  --   zrelo    low part of the real part for the number in z;
  --   zimlo    low part of the imaginary part for the number in z.

  -- ON RETURN :
  --   zrehi    high part of the real part of the number z;
  --   zimhi    high part of the imaginary part of the number z.
  --   zrelo    low part of the real part of the number z;
  --   zimlo    low part of the imaginary part of the number z.

  procedure Inner_Product
                ( zrehi,zimhi,zrelo,zimlo : out double_float;
                  xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                  xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                  yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                  yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                  yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                  yimlo : in Standard_Floating_Vectors.Link_to_Vector);

  -- DESCRIPTION :
  --   Computes the inner product of two double double complex vectors
  --   x and y to form the result z, in 4-vector representation.

  -- REQUIRED : all vectors have the same range.

  -- ON ENTRY :
  --   xrehi    high parts of the real parts for the numbers in x;
  --   ximhi    high parts of the imaginary parts for the numbers in x;
  --   xrelo    low parts of the real parts for the numbers in x;
  --   ximlo    low parts of the imaginary parts for the numbers in x;
  --   yrehi    high parts of the real parts for the numbers in y;
  --   yimhi    high parts of the imaginary parts for the numbers in y;
  --   yrelo    low parts of the real parts for the numbers in y;
  --   yimlo    low parts of the imaginary parts for the numbers in y.

  -- ON RETURN :
  --   zrehi    high parts of the real parts of the inner product;
  --   zimhi    high parts of the imaginary parts of the inner product.
  --   zrelo    low parts of the real parts of the inner product;
  --   zimlo    low parts of the imaginary parts of the inner product.

  procedure Multiply
              ( xrehi : in Standard_Floating_Vectors.Link_to_Vector;
                ximhi : in Standard_Floating_Vectors.Link_to_Vector;
                xrelo : in Standard_Floating_Vectors.Link_to_Vector;
                ximlo : in Standard_Floating_Vectors.Link_to_Vector;
                yrehi : in Standard_Floating_Vectors.Link_to_Vector;
                yimhi : in Standard_Floating_Vectors.Link_to_Vector;
                yrelo : in Standard_Floating_Vectors.Link_to_Vector;
                yimlo : in Standard_Floating_Vectors.Link_to_Vector;
                zrehi : in Standard_Floating_Vectors.Link_to_Vector;
                zimhi : in Standard_Floating_Vectors.Link_to_Vector;
                zrelo : in Standard_Floating_Vectors.Link_to_Vector;
                zimlo : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Convolutes two complex vectors in double double precision,
  --   in the 4-vector representation.

  -- REQUIRED : all vectors have the same range, starting at zero.

  -- ON ENTRY :
  --   xrehi    high parts of the real parts for the numbers in x;
  --   ximhi    high parts of the imaginary parts for the numbers in x;
  --   xrelo    low parts of the real parts for the numbers in x;
  --   ximlo    low parts of the imaginary parts for the numbers in x;
  --   yrehi    high parts of the real parts for the numbers in y;
  --   yimhi    high parts of the imaginary parts for the numbers in y;
  --   yrelo    low parts of the real parts for the numbers in y;
  --   yimlo    low parts of the imaginary parts for the numbers in y.

  -- ON RETURN :
  --   zrehi    high parts of the real parts of the convolution;
  --   zimhi    high parts of the imaginary parts of the convolution;
  --   zrelo    low parts of the real parts of the convolution;
  --   zimlo    low parts of the imaginary parts of the convolution.

end DoblDobl_Vector_Splitters;
