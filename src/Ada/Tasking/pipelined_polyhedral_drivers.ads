with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;

package Pipelined_Polyhedral_Drivers is

-- DESCRIPTION :
--   The procedures in this package define an interface to the
--   pipelined multitasked polyhedral path trackers,
--   in standard double, double double, and quad double precision.
--   The drivers are to be called in Drivers_for_MixedVol_Algorithm.

-- SILENT VERSIONS, ON ORDINARY POLYNOMIAL SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p,
  --   in standard double, double double, or quad double precision,
  --   with pipelined production of mixed cells, interleaved with
  --   the tracking of the paths defined by a polyhedral homotopy.
  --   These procedures are wrappers to the driver procedures on
  --   Laurent polynomial systems below.

  -- REQUIRED : nt >= 2.

  -- ON ENTRY :
  --   nt       number of tasks, which includes one task to produce
  --            the mixed cells, and the others to track the paths;
  --   p        a polynomial system.

  -- ON RETURN :
  --   mv       mixed volume of the tuple of Newton polytopes
  --            spanned by the supports of p;
  --   q        a random coefficient system,
  --            with as many solutions as the mixed volume;
  --   qsols    solution to a random coefficient system,
  --            if all went well, then Length_Of(qsols) = mv.

-- WITH OUTPUT TO FILES, ON ORDINARY POLYNOMIAL SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p,
  --   in standard double, double double, or quad double precision,
  --   with pipelined production of mixed cells, interleaved with
  --   the tracking of the paths defined by a polyhedral homotopy.
  --   These procedures are wrappers to the driver procedures on
  --   Laurent polynomial systems below.

  -- REQUIRED : nt >= 2.

  -- ON ENTRY :
  --   file     to write timings and other results;
  --   cfile    the file to write a regular mixed-cell configuration on,
  --            will be used only if misufile equals true;
  --   qfile    the file to write a random coefficient start system on;
  --   nt       number of tasks, which includes one task to produce
  --            the mixed cells, and the others to track the paths;
  --   misufile indicates whether to write the mixed-cell configuration
  --            to a separate file;
  --   contrep  indicates whether output is wanted during continuation;
  --   p        a polynomial system.

  -- ON RETURN :
  --   mv       mixed volume of the tuple of Newton polytopes
  --            spanned by the supports of p;
  --   q        a random coefficient system,
  --            with as many solutions as the mixed volume;
  --   qsols    solution to a random coefficient system,
  --            if all went well, then Length_Of(qsols) = mv.

-- WITHOUT OUTPUT, ON LAURENT POLYNOMIAL SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p,
  --   in standard double, double double, or quad double precision,
  --   with pipelined production of mixed cells, interleaved with
  --   the tracking of the paths defined by a polyhedral homotopy.

  -- REQUIRED : nt >= 2.

  -- ON ENTRY :
  --   nt       number of tasks, which includes one task to produce
  --            the mixed cells, and the others to track the paths;
  --   p        a polynomial system.

  -- ON RETURN :
  --   mv       mixed volume of the tuple of Newton polytopes
  --            spanned by the supports of p;
  --   q        a random coefficient system,
  --            with as many solutions as the mixed volume;
  --   qsols    solution to a random coefficient system,
  --            if all went well, then Length_Of(qsols) = mv.

-- WITH OUTPUT TO FILES, ON LAURENT POLYNOMIAL SYSTEMS :

  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( file,cfile,qfile : in file_type; nt : in integer32;
                misufile,contrep : in boolean;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p,
  --   in standard double, double double, or quad double precision,
  --   with pipelined production of mixed cells, interleaved with
  --   the tracking of the paths defined by a polyhedral homotopy.

  -- REQUIRED : nt >= 2.

  -- ON ENTRY :
  --   file     to write timings and other results;
  --   cfile    the file to write a regular mixed-cell configuration on,
  --            will be used only if misufile equals true;
  --   qfile    the file to write a random coefficient start system on;
  --   nt       number of tasks, which includes one task to produce
  --            the mixed cells, and the others to track the paths;
  --   misufile indicates whether to write the mixed-cell configuration
  --            to a separate file;
  --   contrep  indicates whether output is wanted during continuation;
  --   p        a polynomial system.

  -- ON RETURN :
  --   mv       mixed volume of the tuple of Newton polytopes
  --            spanned by the supports of p;
  --   q        a random coefficient system,
  --            with as many solutions as the mixed volume;
  --   qsols    solution to a random coefficient system,
  --            if all went well, then Length_Of(qsols) = mv.

-- WITH LIFTING BOUND FOR THE ARTIFICIAL ORIGIN :

  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                lq : out Standard_Complex_Laur_Systems.Laur_Sys;
                q : out Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : out Standard_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                lq : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure Pipelined_Polyhedral_Homotopies
              ( nt : in integer32;
                stable : in boolean; stlb : in double_float;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                lq : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                qsols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Driver to the polyhedral homotopies to create
  --   a random coefficient start system to solve the system p,
  --   in standard double, double double, or quad double precision,
  --   with pipelined production of mixed cells, interleaved with
  --   the tracking of the paths defined by a polyhedral homotopy.
  --   This version allows to provide a lifting bound for the
  --   artificial origin in case the stable mixed volume is requested.

  -- REQUIRED : nt >= 2.

  -- ON ENTRY :
  --   nt       number of tasks, which includes one task to produce
  --            the mixed cells, and the others to track the paths;
  --   stable   indicates whether the stable mixed volume is wanted;
  --   stlb     stable lifting bound if stable, otherwise 0.0;
  --   p        a polynomial system.

  -- ON RETURN :
  --   r        number of distinct supports;
  --   mtype    type of mixture stores frequency of each support;
  --   perm     permutation of the supports when computing mixture type;
  --   lif      configuration of lifted supports;
  --   mcc      mixed cell configuration, which contains all mixed cells,
  --            those with and without artificial origin, if stable;
  --   mv       mixed volume of the tuple of Newton polytopes
  --            spanned by the supports of p;
  --   lq       the random coefficient system q as Laurent system;
  --   q        a random coefficient system,
  --            with as many solutions as the mixed volume;
  --   qsols    solution to a random coefficient system,
  --            if all went well, then Length_Of(qsols) = mv.

end Pipelined_Polyhedral_Drivers;
