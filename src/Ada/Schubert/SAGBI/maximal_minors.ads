with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;

procedure Maximal_Minors ( file : in file_type;
                           n,d : in natural32; mat : in Matrix;
                           min,max : out double_float );

-- DESCRIPTION :
--   Computes all maximal minors of a (nxd)-matrix mat, d < n,
--   and computes the minimal and maximal determinant in absolute value.
