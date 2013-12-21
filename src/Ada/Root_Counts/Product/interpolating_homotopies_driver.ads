with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;

procedure Interpolating_Homotopies_Driver
              ( file : in file_type; p : in Poly_Sys; z : in Partition;
                b : in out natural32; q : out Poly_Sys;
                qsols : in out Solution_List );

-- DESCRIPTION :
--   This is an interactive driver for the construction of an interpolating
--   homotopy based on an m-homogeneous Bezout number.

-- ON ENTRY :
--   file       to write diagnostics on;
--   p          the polynomial system;
--   z          partition of the set of unknowns of p;
--   b          an m-homogeneous Bezout number.

-- ON RETURN :
--   b          number of interpolating roots;
--   q          an m-homogeneous start system;
--   qsols      solutions of q.
