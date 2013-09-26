with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

procedure Driver_to_Factor_Components
            ( file : in file_type;
              ep : in Standard_Complex_Poly_Systems.Poly_Sys;
              sols : in Standard_Complex_Solutions.Solution_List;
              dim : in natural32 );

-- DESCRIPTION :
--   Offers a menu of different factorization methods for a pure
--   dimensional component, given as embedding with generic points.

-- ON ENTRY :
--   file     to write diagnostics and results;
--   ep       embedding of polynomial system;
--   sols     list of generic points;
--   dim      dimension of the solution component.
