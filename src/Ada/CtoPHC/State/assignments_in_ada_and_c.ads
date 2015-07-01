with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Quad_Double_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with C_Integer_Arrays;                   use C_Integer_Arrays;
with C_Double_Arrays;                    use C_Double_Arrays;

package Assignments_in_Ada_and_C is

-- DESCRIPTION :
--   This package provides assignment operators between Ada and C data.

  procedure Assign ( this : in integer32; to_that : in C_intarrs.Pointer );
  procedure Assign ( this : in C_intarrs.Pointer; to_that : out integer32 );
  procedure Assign ( this : in double_float; to_that : in C_dblarrs.Pointer );
  procedure Assign ( this : in C_dblarrs.Pointer; to_that : out double_float );

  -- DESCRIPTION :
  --   Assigns what is in this to the location pointed by to_that.

  procedure Assign ( this : in double_double;
                     to_that : in C_dblarrs.Pointer );
  procedure Assign ( this : in C_dblarrs.Pointer;
                     to_that : out double_double );
  procedure Assign ( this : in quad_double;
                     to_that : in C_dblarrs.Pointer );
  procedure Assign ( this : in C_dblarrs.Pointer;
                     to_that : out quad_double );

  -- DESCRIPTION :
  --   Assigns the consecutive locations for double doubles and quad doubles
  --   of what is in this to the location pointed by to_that.

  procedure Assign ( ada_cf : in Standard_Complex_Numbers.Complex_Number;
                     c_cf : in C_dblarrs.Pointer );
  procedure Assign ( ada_cf : in DoblDobl_Complex_Numbers.Complex_Number;
                     c_cf : in C_dblarrs.Pointer );
  procedure Assign ( ada_cf : in QuadDobl_Complex_Numbers.Complex_Number;
                     c_cf : in C_dblarrs.Pointer );

  -- DESCRIPTION :
  --   Assigns the complex number in ada_cf to the array in c_cf,
  --   copying real and imaginary part of ada_cf to c_cf
  --   in consecutive locations.

  procedure Assign ( c_cf : in C_dblarrs.Pointer;
                     ada_cf : out Standard_Complex_Numbers.Complex_Number );
  procedure Assign ( c_cf : in C_dblarrs.Pointer;
                     ada_cf : out DoblDobl_Complex_Numbers.Complex_Number );
  procedure Assign ( c_cf : in C_dblarrs.Pointer;
                     ada_cf : out QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Assigns to ada_cf the values of real and imaginary parts in c_cf.
  --   The values for double double and quad doubles are assumed to be
  --   in consecutive locations in memory.

  procedure Assign ( ada_d : in Standard_Natural_Vectors.Vector;
                     c_d : in C_intarrs.Pointer );
  procedure Assign ( ada_d : in Standard_Integer_Vectors.Vector;
                     c_d : in C_intarrs.Pointer );
  procedure Assign ( ada_d : in Standard_Floating_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer );
  procedure Assign ( ada_d : in Double_Double_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer );
  procedure Assign ( ada_d : in Quad_Double_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer );
  procedure Assign ( ada_d : in Standard_Complex_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer );
  procedure Assign ( ada_d : in DoblDobl_Complex_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer );
  procedure Assign ( ada_d : in QuadDobl_Complex_Vectors.Vector;
                     c_d : in C_dblarrs.Pointer );

  -- DESCRIPTION :
  --   Assigns the content of the Ada vector to the C vector.
  --   The complex numbers in ada_d are written to c_d as a sequence
  --   of real and imaginary double floats.
  --   Double double and quad double numbers are written as sequences
  --   of respectively two and four doubles, in order of significance.

  -- REQUIRED : ada_d'range = 1..n, for some n.

  procedure Assign ( v_n : in natural32; c_d : in C_intarrs.Pointer; 
                     ada_d : out Standard_Natural_Vectors.Vector );
  procedure Assign ( v_n : in natural32; c_d : in C_intarrs.Pointer; 
                     ada_d : out Standard_Integer_Vectors.Vector );
  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out Standard_Floating_Vectors.Vector );
  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out Standard_Complex_Vectors.Vector );
  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out DoblDobl_Complex_Vectors.Vector );
  procedure Assign ( v_n : in natural32; c_d : in C_dblarrs.Pointer;
                     ada_d : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Assigns the content of the C vector to the Ada vector.
  --   For complex vectors, it copies the consecutive real and
  --   and imaginary parts in c_d to ada_d.
  --   When writing to a vector of double double or quad double
  --   complex numbers, the doubles are interpreted as parts of
  --   a double double or quad doubles, with their most significant
  --   parts first.

  -- REQUIRED : ada_d'range = 1..N, for some N;
  --    for complex vectors: ada_d'range = 1..v_n/2, for some v_n;
  --    for double double vectors: ada_d'range = 1..v_n/2; and
  --    for quad double vectors: ada_d'range = 1..v_n/4; and
  --    for double double complex vectors: ada_d'range = 1..v_n/4; and
  --    for quad double complex vectors: ada_d'range = 1..v_n/8.

  function C_Integer_Array_to_String
             ( n : natural32; v : C_Integer_Array ) return String;

  -- DESCRIPTION :
  --   Converts an array of n integers into a string of n characters.

  function Pad_with_Spaces ( n : natural32; s : string ) return string;

  -- DESCRIPTION :
  --   If s'last >= n, then the string s is returned, otherwise,
  --   if s'last < n, then the string on return will be of length n,
  --   contain s padded with spaces.

  function String_to_C_Integer_Array
             ( n : natural32; s : string ) return C_Integer_Array;
  function String_to_Integer_Vector
             ( s : string ) return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts a string into an array of integer C type numbers.

end Assignments_in_Ada_and_C;
