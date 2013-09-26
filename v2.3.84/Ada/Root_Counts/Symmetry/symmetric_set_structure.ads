with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Symmetry_Group;                     use Symmetry_Group;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Symmetric_Set_Structure is

-- DESCRIPTION :
--   The aim of this package is the construction of a symmetric
--   start system, given a symmetric set structure.

-- CONSTRUCTORS :

  procedure Equivariant_Start_System
                ( n : in natural32; g : in List_of_Permutations;
                  fail : out boolean );

  -- DESCRIPTION :
  --   Constructs an equivariant linear-product start system.
  --   When not fail on return, the package Random_Product_System
  --   contains the data for generating the polynomial system.

  -- REQUIRED :
  --   The data of the package Set_Structure may not be empty!

  -- ON ENTRY :
  --   n          the dimension of the problem;
  --   g          the list of generating permutations.

  -- ON RETURN :
  --   fail       if true, then the set structure was not equivariant.

  procedure Symmetric_Start_System
                ( n,bb : in natural32; lp : in List;
                  v,w : in List_of_Permutations;
                  notsymmetric,degenerate : out boolean );

  -- DESCRIPTION :
  --   After calling this routine, the package Random_Product_System
  --   contains the data for a symmetric random product system,
  --   when notsymmetric and degenerate are false on return.

  -- REQUIRED :
  --   The data of the package Set_Structure may not be empty!

  -- ON ENTRY :
  --   n          the dimension of the problem;
  --   bb         the Bezout number based on the set structure;
  --   lp         list of positions indicating the acceptable
  --              classes in the set structure;
  --   v,w        representations of the symmetry group.

  -- ON RETURN :
  --   notsymmetric  is true if the set structure is not symmetric;
  --   degenerate    is true if the set structure is degenerate.

-- SELECTORS 

  procedure Write_Covering;
  procedure Write_Covering ( file : in file_type );

  procedure Write_Templates ( n : in natural32 );
  procedure Write_Templates ( file : in file_type; n : in natural32 );

  -- DESCRIPTION :
  --   These procedure write an intermediate data structures
  --   in the construction of a symmetric start system.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   All allocated memory space will be freed.

end Symmetric_Set_Structure;
