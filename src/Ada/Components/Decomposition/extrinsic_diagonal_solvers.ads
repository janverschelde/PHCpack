with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;                      use Symbol_Table;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Extrinsic_Diagonal_Solvers is

-- DESCRIPTION :
--   This package offers drivers to compute witness sets on the intersection
--   of positive dimensional components by extrinsic diagonal homotopies.
--   Note that the routines in this package do not solve, but provide the
--   homotopies on which standard continuation can be applied directly.

  procedure Save_Target_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Save_Target_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Save_Target_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys );

  -- DESCRIPTION :
  --   The user is prompted for a file name to save the target system,
  --   in standard double, double double, or quad double precision,
  --   in the homotopy to start the diagonal cascade.

  procedure Save_Start_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List );
  procedure Save_Start_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Save_Start_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   The user is prompted for a file name to save the start system,
  --   in standard double, double double, or quad double precision,
  --   in the homotopy to start the diagonal cascade.

  procedure Test_Solutions
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Solution_List );
  procedure Test_Solutions
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Test_Solutions
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Evaluates the solutions in s at the system p.

  procedure Standard_Randomize_System;
  procedure DoblDobl_Randomize_System;
  procedure QuadDobl_Randomize_System;

  -- DESCRIPTION :
  --   Interactive driver to randomize a system so that the number
  --   of equations equals the co-dimension of the solution componen,
  --   in standard double, double double, or quad double precision.

  procedure Randomize_System;

  -- DESCRIPTION :
  --   Prompts the user for the precision (d, dd, or qd), and then
  --   calls the proper randomization procedure so that the system
  --   has the right number of equations with respect to the solutions.

  procedure Build_Cascade_Homotopy
              ( file : in file_type;
                p1e,p2e : in Standard_Complex_Poly_Systems.Poly_Sys;
                dim1,dim2 : in natural32;
                sols1e,sols2e : in Standard_Complex_Solutions.Solution_List;
                s1e,s2e : in Array_of_Symbols );
  procedure Build_Cascade_Homotopy
              ( file : in file_type;
                p1e,p2e : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                dim1,dim2 : in natural32;
                sols1e,sols2e : in DoblDobl_Complex_Solutions.Solution_List;
                s1e,s2e : in Array_of_Symbols );
  procedure Build_Cascade_Homotopy
              ( file : in file_type;
                p1e,p2e : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                dim1,dim2 : in natural32;
                sols1e,sols2e : in QuadDobl_Complex_Solutions.Solution_List;
                s1e,s2e : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Builds the homotopy to start the cascade for all components of
  --   the intersection of two solution components (of dimensions dim1
  --   and dim2) defined by the embedded systems ep1 and ep2,
  --   in standard double, double double, and quad double precision.

  -- ON ENTRY :
  --   p1e      1st embedded polynomial system;
  --   p2e      2nd embedded polynomial system;
  --   dim1     dimension of the component defined by the 1st system;
  --   dim2     dimension of the component defined by the 2nd system;
  --   sols1e   witness points on the 1st solution component;
  --   sols1e   witness points on the 2nd solution component;
  --   s1e      symbols used in the 1st polynomial system;
  --   s2e      symbols used in the 2nd polynomial system.

  procedure Standard_Diagonal_Cascade;
  procedure DoblDobl_Diagonal_Cascade;
  procedure QuadDobl_Diagonal_Cascade;

  -- DESCRIPTION :
  --   Interactive driver to build a homotopy to start a cascade
  --   of homotopies to compute all positive dimensional components
  --   of the intersection of two components,
  --   in standard double, double double, or quad double precision.

  procedure Build_Diagonal_Cascade;

  -- DESCRIPTION :
  --   Prompts the user for the precision and then calls then
  --   Standard_Diagonal_Cascade, DoblDobl_Diagonal_Cascade,
  --   or QuadDobl_Diagonal_Cascade, depending on th precision.

  procedure Collapse_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                dim,add2dim : in natural32;
                r : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Collapse_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                dim,add2dim : in natural32;
                r : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure Collapse_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                dim,add2dim : in natural32;
                r : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Removes the duplicate variables from p and sols,
  --   in standard double, double double, or quad double precision,
  --   embedding the result as a component of dimension dim+addtodim.

  -- ON ENTRY :
  --   p        an embedded polynomial system, as many slacks as dim;
  --   sols     solutions to p, with the proper embedding;
  --   dim      number of slack variables currently in the embedding;
  --   add2dim  number of slack variables that must be added.

  -- ON RETURN :
  --   sols     solutions with dim + add2dim slack variables;
  --   r        system with the diagonal removed, embedded with
  --            as many as dim + add2dim slack variables.

  procedure Standard_Collapse_Diagonal_System;
  procedure DoblDobl_Collapse_Diagonal_System;
  procedure QuadDobl_Collapse_Diagonal_System;

  -- DESCRIPTION :
  --   Once the target system in the diagonal homotopy has been solved,
  --   this procedure transforms the problem to the original coordinates,
  --   in standard double, double double, or quad double precision.

  procedure Collapse_Diagonal_System;

  -- DESCRIPTION :
  --   Prompts the user for the precision and then calls the proper
  --   procedure to collapse the problem onto the original coordinates.

end Extrinsic_Diagonal_Solvers;
