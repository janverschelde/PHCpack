with text_io;                            use text_io;
with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;

package Drivers_for_Homotopy_Creation is

-- DESCRIPTION :
--   Menu-driven creation of an artificial-parameter polynomial homotopy.

  procedure Default_Homotopy_Settings
                ( d,k : out natural32; a,t : out Complex_Number );
  procedure Default_Homotopy_Settings
                ( d,k : out natural32; a,t : out Complex_Number;
                  proj : out boolean );

  -- DESCRIPTION :
  --   Sets the default settings of the homotopy parameters.

  -- ON RETURN :
  --   d          number of decimal places, by default equals 16;
  --   k          relaxation constant, by default equals 2;
  --   a          randomly chosen complex number of radius 1;
  --   t          1 is default target value of the continuation parameter;
  --   proj       no projective transformation, proj = false.

  procedure Menu_for_Homotopy_Settings
                ( file : in file_type; d,k : in out natural32;
                  a,t : in out Complex_Number; proj : in out boolean );
  procedure Menu_for_Homotopy_Settings
                ( file : in file_type; d,k : in out natural32;
                  a,t : in out Complex_Number );
  procedure Menu_for_Homotopy_Settings
                ( file : in file_type; qd : in boolean;
                  k : in out natural32; a,t : in out Complex_Number );

  -- DESCRIPTION :
  --   Interactive setting of the homotopy parameters.
  --   If proj is left out, then no projective transformation is set.
  --   If qd, then quad doubles, else double doubles will be used.

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 p,q : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 target : out Complex_Number; verbose : in integer32 := 0 );
  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 p,q : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 target : out Complex_Number; verbose : in integer32 := 0 );
  procedure Driver_for_Homotopy_Construction
               ( file : in file_type; dp : in natural32;
                 p,q : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 target : out Complex_Number; verbose : in integer32 := 0 );
  procedure Driver_for_Homotopy_Construction
               ( file : in file_type; p,q : in out Laur_Sys;
                 target : out Complex_Number; deci : in out natural32;
                 verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Drivers to construct polynomial homotopies in double, double double,
  --   quad double, and arbitrary multiprecision.

  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 ls : in String_Splitters.Link_to_Array_of_Strings;
                 p,q : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols : in out Solution_List; target : out Complex_Number;
                 deci : in out natural32; verbose : in integer32 := 0 );
  procedure Driver_for_Homotopy_Construction
               ( file : in file_type;
                 ls : in String_Splitters.Link_to_Array_of_Strings;
                 p,q : in out Standard_Complex_Laur_Systems.Laur_Sys;
                 qsols : in out Solution_List; target : out Complex_Number;
                 deci : in out natural32; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   This is an interactive driver for the construction of an artificial
  --   parameter homotopy.  The user can ask to homogenize the polynomial
  --   system and can determine the homotopy parameters and the target
  --   value for the  continuation parameter.
  --   Note that the starting value for the continuation parameter is stored
  --   with the solutions.

  -- ON ENTRY :
  --   file      for output of the settings of the homotopy;
  --   p         target system;
  --   q         start system;
  --   qsols     solutions of the start system q;
  --   deci      may have a preset value (if different from zero);
  --   verbose   the verbose level.

  -- ON RETURN :
  --   p         target system, eventually with projective transformation;
  --   q         start system, eventually with projective transformation;
  --   qsols     start solutions, eventually with projective transformation;
  --   target    target value of the continuation parameter;
  --   deci      number of decimal places.

end Drivers_for_Homotopy_Creation;
