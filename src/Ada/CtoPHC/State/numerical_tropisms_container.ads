with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_Vectors;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Quad_Double_Vectors;
with Quad_Double_VecVecs;

package Numerical_Tropisms_Container is

-- DESCRIPTION :
--   A numerical tropisms is a numerical approximation for a tropism,
--   computed typically via a polyhedral end game.
--   The data to represent a numerical tropism is a triplet:
--   1) an estimated winding number, stored as an integer;
--   2) real floating-point coordinates for the exponents of
--   the leading powers of the Puiseux series of the path; and
--   3) one float: an estimate for the error on those exponents.
--   This package provides a data structure to manage numerical tropisms
--   in standard double, double double, and quad double precision.

-- CONSTRUCTOES :

  procedure Standard_Initialize 
              ( w : in Standard_Integer_Vectors.Vector;
                v : in Standard_Floating_VecVecs.VecVec;
                e : in Standard_Floating_Vectors.Vector );
  procedure DoblDobl_Initialize 
              ( w : in Standard_Integer_Vectors.Vector;
                v : in Double_Double_VecVecs.VecVec;
                e : in Double_Double_Vectors.Vector );
  procedure QuadDobl_Initialize 
              ( w : in Standard_Integer_Vectors.Vector;
                v : in Quad_Double_VecVecs.VecVec;
                e : in Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Initializes the data structures to store numerical tropisms
  --   in standard double, double double, or quad double precision.

  -- REQUIRED : w'range = v'range = e'range = 1..#solutions.

  procedure Store_Standard_Tropism
              ( k : integer32; w : in integer32;
                v : in Standard_Floating_Vectors.Vector;
                e : in double_float ); 

  -- DESCRIPTION :
  --   Stores the triplet (w, v, e) as the values as the k-th tropism,
  --   computed in standard double precision.

  -- REQUIRED : Standard_Initialize was executed and Standard_Size >= k.
  --   Moreover, the range of v fits the allocated data structures
  --   for the directions.

  procedure Store_DoblDobl_Tropism
              ( k : integer32; w : in integer32;
                v : in Double_Double_Vectors.Vector;
                e : in double_double ); 

  -- DESCRIPTION :
  --   Stores the triplet (w, v, e) as the values as the k-th tropism,
  --   computed in double double precision.

  -- REQUIRED : DoblDobl_Initialize was executed and DoblDobl_Size >= k.
  --   Moreover, the range of v fits the allocated data structures
  --   for the directions.

  procedure Store_QuadDobl_Tropism
              ( k : integer32; w : in integer32;
                v : in Quad_Double_Vectors.Vector;
                e : in quad_double ); 

  -- DESCRIPTION :
  --   Stores the triplet (w, v, e) as the values as the k-th tropism,
  --   computed in quad double precision.

  -- REQUIRED : QuadDobl_Initialize was executed and QuadDobl_Size >= k.
  --   Moreover, the range of v fits the allocated data structures
  --   for the directions.

-- SELECTORS :

  procedure Standard_Retrieve
              ( w : out Standard_Integer_Vectors.Link_to_Vector;
                v : out Standard_Floating_VecVecs.Link_to_VecVec;
                e : out Standard_Floating_Vectors.Link_to_Vector );
  procedure DoblDobl_Retrieve
              ( w : out Standard_Integer_Vectors.Link_to_Vector;
                v : out Double_Double_VecVecs.Link_to_VecVec;
                e : out Double_Double_Vectors.Link_to_Vector );
  procedure QuadDobl_Retrieve
              ( w : out Standard_Integer_Vectors.Link_to_Vector;
                v : out Quad_Double_VecVecs.Link_to_VecVec;
                e : out Quad_Double_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Returns the entire data structure to store numerical tropisms
  --   in standard double, double double, or quad double precision.

  function Standard_Size return integer32;
  function DoblDobl_Size return integer32;
  function QuadDobl_Size return integer32;

  -- DESCRIPTION :
  --   Returns the number of numerical tropisms in standard double,
  --   double double, or quad double precision.
  --   Returns 0 if the data structures are cleared or not initialized.

  function Standard_Dimension return integer32;
  function DoblDobl_Dimension return integer32;
  function QuadDobl_Dimension return integer32;

  -- DESCRIPTION :
  --   Returns the length of the numerical tropisms in standard double,
  --   double double, or quad double precision.
  --   Returns 0 if the data structures are cleared or not initialized.

  procedure Standard_Retrieve_Tropism
              ( k : in integer32; w : out integer32;
                v : out Standard_Floating_Vectors.Vector;
                e : out double_float );

  -- DESCRIPTION :
  --   Retrieves the triplet (w, v, e) as the values as the k-th tropism,
  --   computed in standard double precision.

  -- REQUIRED : Standard_Initialize was executed and Standard_Size >= k.
  --   Moreover, the range of v fits the allocated data structures
  --   for the directions.

  procedure DoblDobl_Retrieve_Tropism
              ( k : in integer32; w : out integer32;
                v : out Double_Double_Vectors.Vector;
                e : out double_double );

  -- DESCRIPTION :
  --   Retrieves the triplet (w, v, e) as the values as the k-th tropism,
  --   computed in double double precision.

  -- REQUIRED : DoblDobl_Initialize was executed and DoblDobl_Size >= k.
  --   Moreover, the range of v fits the allocated data structures
  --   for the directions.

  procedure QuadDobl_Retrieve_Tropism
              ( k : in integer32; w : out integer32;
                v : out Quad_Double_Vectors.Vector;
                e : out quad_double );

  -- DESCRIPTION :
  --   Retrieves the triplet (w, v, e) as the values as the k-th tropism,
  --   computed in quad double precision.

  -- REQUIRED : QuadDobl_Initialize was executed and QuadDobl_Size >= k.
  --   Moreover, the range of v fits the allocated data structures
  --   for the directions.

-- DESTRUCTORS :

  procedure Standard_Clear;
  procedure DoblDobl_Clear;
  procedure QuadDobl_Clear;

  -- DESCRIPTION :
  --   Deallocates the data structures used to store numerical tropisms
  --   in standard double, double double, or quad double precision.

end Numerical_Tropisms_Container;
