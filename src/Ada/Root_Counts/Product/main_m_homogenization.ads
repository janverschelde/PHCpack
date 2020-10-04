with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Standard_Complex_Prod_Systems;      use Standard_Complex_Prod_Systems;

package Main_m_Homogenization is

-- DESCRIPTION :
--   The main interactive procedure to compute m-homogeneous Bezout numbers
--   are defined by this package.

  procedure m_Homogenization_Info;

  -- DESCRIPTION :
  --   Displays information about m-homogenization on screen.

  procedure Menu_Prompt
              ( choice : out character; bz : in natural64;
                np,nz : in natural32; z : in Partition );

  -- DESCRIPTION :
  --   Shows a menu with the current m-homogeneous Bezout number
  --   and prompts either to exit or to compute more Bezout numbers.

  -- ON ENTRY :
  --   bz       the current m-homogeneous Bezout number;
  --   np       total number of partitions;
  --   nz       number of sets in the current partition;
  --   z        current partition of the set of unknowns.

  -- ON RETURN :
  --   choice   choice made by the user.

  procedure Menu_Handler
              ( file : in file_type; choice : in character;
                p : in Poly_Sys;
                bz : in out natural64; np,n : in natural32;
                nz : in out natural32; z : in out Partition );

  -- DESCRIPTION :
  --   Handles the choice of the user.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   choice   choice of method, is the outcome of Menu_Prompt;
  --   p        polynomial system for which bz is the Bezout number;
  --   bz       current Bezout number with partition z;
  --   np       total number of partitions;
  --   n        number of variables in the system p;
  --   nz       number of sets in current partition z;
  --   z        current partition is z(1..nz).

  procedure Define_Partition
              ( file : in file_type; bz : in out natural64;
                p : in Poly_Sys; np,n : in natural32;
                nz : in out natural32; z : in out Partition );

  -- DESCRIPTION :
  --   Interactive procedure to define the partition of the set
  --   of variables of a polynomial system.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   bz       current m-homogeneous Bezout number for the partition z;
  --   p        a polynomial system;
  --   np       total number of partitions;
  --   n        number of variables;
  --   nz       number of sets in the partition z;
  --   z        a partition of the set of variables.

  -- ON RETURN :
  --   bz       updated m-homogeneous Bezout number;
  --   nz       updated number of sets in z;
  --   z        updated partition of the set of variables.

  procedure Construct_Start_System
              ( file : in file_type; p : in Poly_Sys;
                bz : in natural64; z : in Partition;
                q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Given a partition of the set of variables,
  --   constructs a start system with the m-homogeneous structure.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   bz       m-homogeneous Bezout number for the partition z;
  --   p        a polynomial system;
  --   nz       number of sets in the partition z;
  --   z        a partition of the set of variables.

  -- ON RETURN :
  --   q        a start system with same m-homogeneous degree
  --            structure as the given system p;
  --   qsols    solutions of q.

  procedure Write_Results ( file : in file_type; bz : in natural64;
                            nz : in natural32; z : in partition  );

  -- DESCRIPTION :
  --   Writes the m-homogeneous Bezout number bz,
  --   with the partition z and number of sets in nz
  --   to a file opened for output.

  procedure Save_Results ( rq : in Prod_Sys; m : in natural32 );

  -- DESCRIPTION :
  --   Prompts for the name of a file to save a linear-product
  --   start system rq for an m-homogeneous Bezout number.

  procedure Save_Results ( qq : in Poly_Sys; qqsols : in Solution_List );

  -- DESCRIPTION :
  --   Prompts for the name of a file to save a start system qq and
  --   all its solutions qqsols to file.

  procedure Main ( file : in file_type; p : in Poly_Sys;
                   b : in out natural64;
                   q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Computation of an m-homogeneous Bezout number with the option
  --   of constructing an m-homogeneous start system.

  -- ON ENTRY :
  --   file        output file;
  --   p           a polynomial system.

  -- ON RETURN :
  --   b           m-homogeneous Bezout number;
  --   q           m-homogeneous start system;
  --   qsols       solutions of q.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for a system and an output file and then computes.

end Main_m_Homogenization;
