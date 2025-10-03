with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package DEMiCs_Translated is

-- DESCRIPTION :
--   This package provides an interface to the translated DEMiCs,
--   to compute all mixed cells by dynamic enumeration.
--   DEMiCs was developed by Tomohiko Mizutani, Akiko Takeda, and
--   Masakazu Kojima and licensed under GNU GPL Version 2 or higher.
--   The translated DEMiCs is a literal translation of the C++ into Ada.

  function Mixed_Volume
             ( p : Poly_Sys; seednbr : integer32 := 0;
               stablemv : boolean := false;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32;
  function Mixed_Volume
             ( p : Laur_Sys; seednbr : integer32 := 0;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume of the polynomials in p.
  --   This function is equivalent to the main function of DEMiCs.
  --   If the mixed cells are required, then the functions
  --   Mixed_Labels and Mixed_Cells need to be called.

  -- ON ENTRY :
  --   p       a system of polynomial equations;
  --   seednbr is the number to seed the random number generator,
  --           if nonzero;
  --   stablemv flags if the stable mixed volume needs to be computed;
  --   userlifting indicates that the user will be prompted to provide
  --           lifting values interactively for each point in the supports,
  --           otherwise the lifting is randomly generated.

  function Mixed_Labels
             ( p : Poly_Sys; monitor : boolean := true;
               seednbr : integer32 := 0;
               stablemv : boolean := false;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32;
  function Mixed_Labels
             ( p : Laur_Sys; monitor : boolean := true;
               seednbr : integer32 := 0;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32;

  -- DESCRIPTION :
  --   Returns the mixed volume and the labels to the mixed cell indices.
  --   If monitor is true, then the cell is written to screen each time
  --   the cell is added to the data stored in demics_output_cells.

  -- SIDE EFFECT :
  --   The package DEMiCs_Output_Cells contains the labels.

  function Mixed_Cells ( vrblvl : integer32 := 0 ) return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Returns the mixed-cell configuration computed by Mixed_Labels.

  -- REQUIRED :
  --   Mixed_Labels has been executed successfully.

-- INTERFACE :
--   The procedures below have the same specifications as the
--   original interface to DEMiCs.


  procedure Extract_Supports 
              ( p : in Poly_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );
  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the supports and computes the type of mixture.

  procedure Extract_Supports 
              ( p : in Poly_Sys;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );
  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Extracts the supports and computes the type of mixture.
  --   In addition, returns the permutation to put the same supports
  --   in consecutive order.

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_VecVecs.Link_to_VecVec;

  -- DESCRIPTION :
  --   Returns random lifting values for each of the support sets,
  --   as arranged in the order of their mixture, defined by mix.
 
  function Flatten ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                   return Standard_Floating_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   Returns a vector with all values in v,
  --   used when the user lifting is provided as a vector of vectors.

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float := 0.0;
                lft : in Standard_Floating_Vectors.Link_to_Vector := null;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Calls DEMiCs on the given supports in sup, of type of mixture mix.

  -- ON ENTRY :
  --   mix      type of mixture of the supports.
  --   sup      supports of a polynomial system,
  --            sorted along the type of mixture;
  --   stlb     stable lifting bound:
  --            if zero, then no stable mixed volume will be computed,
  --            otherwise, stlb is the lifting for the artificial origins;
  --   lft      null if random numbers need to be generated,
  --            otherwise, as many values as points as in sup,
  --            with sufficiently random values for the lifting,
  --            note that is stlb /= 0, then lft is treated as null;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   supp     artificial origins added if stlb /= 0.0.

  procedure Show_Output;

  -- DESCRIPITON :
  --   Shows the output stored in DEMiCs_Output_Cells.

  procedure Process_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : out Mixed_Subdivision; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Converts the output of DEMiCs to lifted supports
  --   and a mixed cell configuration.

  -- REQUIRED :
  --   Mixed_Labels has been executed successfully.

  -- ON ENTRY :
  --   dim      dimension of the points before lifting;
  --   mix      type of mixture;
  --   sup      supports of a polynomial system;
  --   verbose  flag to indicate if extra output is wanted.

  -- ON RETURN :
  --   lif      lifted support sets;
  --   mcc      a regular mixed-cell configuration.

-- DESTRUCTOR :

  procedure Clear ( vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Deallocates all allocated memory.

end DEMiCs_Translated;
