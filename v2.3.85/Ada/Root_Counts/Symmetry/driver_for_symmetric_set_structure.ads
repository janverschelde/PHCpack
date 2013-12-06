with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems; 
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Driver_for_Symmetric_Set_Structure is

  procedure Symmetric_Set_Structure_Info;

  -- DESCRIPTION :
  --   Displays information on symmetric set structures on screen.

  procedure Driver_for_Symmetric_Random_Product_Systems
                  ( file : in file_type; p : in Poly_Sys; q : out Poly_Sys;
                    qsols : out Solution_List; bs : in out natural32;
                    lpos : in out List );

  -- DESCRIPTION :
  --   Interactive driver for the construction of a
  --   (G,V,W)-symmetric random product start system.

  -- ON ENTRY :
  --   file         output file to write diagnostics on;
  --   p            a polynomial system;
  --   bs           Bezout number based on the set structure;
  --   lpos         list of positions.

  -- ON RETURN :
  --   q            symmetric random linear-product start system;
  --   qsols        solutions of q;
  --   bs           Bezout number based on a symmetric set structure;
  --   lpos         list of positions for the new set structure.

end Driver_for_Symmetric_Set_Structure;
