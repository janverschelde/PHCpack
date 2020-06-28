with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Series_Vectors;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;

package Homotopy_Series_Readers is

-- DESCRIPTION :
--   Provides interactive procedures to setup of homotopies of series,
--   in double, double double, and quad double precision.
--   The homotopy is an artificial parameter homotopy
--   or a natural parameter homotopy.

  procedure Standard_Projective_Transformation
              ( target,start
                 : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure DoblDobl_Projective_Transformation
              ( target,start
                 : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys );
  procedure QuadDobl_Projective_Transformation
              ( target,start
                 : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys );

  -- DESCRIPTION :
  --   Transforms the target and start system into homogeneous coordinates,
  --   adding one random linear equation to the target system and Z0 = 1 
  --   to the start system, adding 1 to every start solution,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   target   target system in an artificial-parameter homotopy;
  --   start    start system in an artificial-parameter homotopy.

  -- ON RETURN :
  --   target   target system in homogeneous coordinates
  --            and with one random linear equation added;
  --   start    start system in homogeneous coordinates
  --            and with the equation Z0 = 1 added;

  procedure Standard_Projective_Transformation
              ( target : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Projective_Transformation
              ( target : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Projective_Transformation
              ( target : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Transforms the target, start system, and its start solutions
  --   into homogeneous coordinates, adding one random linear equation to
  --   the target system and Z0 = 1 to the start system, adding 1 to
  --   every start solution.

  -- ON ENTRY :
  --   target   target system in an artificial-parameter homotopy;
  --   start    start system in an artificial-parameter homotopy;
  --   sols     start solutions.

  -- ON RETURN :
  --   target   target system in homogeneous coordinates
  --            and with one random linear equation added;
  --   start    start system in homogeneous coordinates
  --            and with the equation Z0 = 1 added;
  --   sols     solutions extended with 1 as last coordinate.

  procedure Standard_Multi_Projective_Transformation
              ( target,start
                  : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition );
  procedure DoblDobl_Multi_Projective_Transformation
              ( target,start
                  : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition );
  procedure QuadDobl_Multi_Projective_Transformation
              ( target,start
                  : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                m : in natural32; z : in Partition );

  -- DESCRIPTION :
  --   Transforms the target and start system into m-homogeneous coordinates,
  --   adding m random linear equations to the target system and Zi = 1,
  --   for i in 1..m, to the start system, adding 1 to every start solution,
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   target   target system in an artificial-parameter homotopy;
  --   start    start system in an artificial-parameter homotopy;
  --   m        number of sets in the partition;
  --   z        defines the multi-homogenization.

  -- ON RETURN :
  --   target   target system in homogeneous coordinates
  --            and with one random linear equation added;
  --   start    start system in homogeneous coordinates
  --            and with the equation Z0 = 1 added;

  procedure Standard_Multi_Projective_Transformation
              ( target : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out Standard_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition );
  procedure DoblDobl_Multi_Projective_Transformation
              ( target : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition );
  procedure QuadDobl_Multi_Projective_Transformation
              ( target : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                start : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                m : in natural32; z : in Partition );

  -- DESCRIPTION :
  --   Transforms the target, start system, and its start solutions
  --   into homogeneous coordinates, adding one random linear equation to
  --   the target system and Z0 = 1 to the start system, adding 1 to
  --   every start solution.

  -- ON ENTRY :
  --   target   target system in an artificial-parameter homotopy;
  --   start    start system in an artificial-parameter homotopy;
  --   sols     start solutions;
  --   m        number of sets in the partition;
  --   z        defines the multi-homogenization.

  -- ON RETURN :
  --   target   target system in homogeneous coordinates
  --            and with one random linear equation added;
  --   start    start system in homogeneous coordinates
  --            and with the equation Z0 = 1 added;
  --   sols     solutions extended with 1 as last coordinate.

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in Standard_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false );
  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in DoblDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false );
  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32;
                gamma : in QuadDobl_Complex_Numbers.Complex_Number;
                homcrd,rabin : in boolean := false );

  -- DESCRIPTION :
  --   Prompts for a target system, a start system with start solutions.
  --   The target and start system are stored in the Homotopy package.
  --   The given power tpow of the continuation parameter and 
  --   the given accessibility constant gamma will be used.

  -- ON ENTRY :
  --   tpow     power of the continuation parameter
  --            in the artificial parameter homotopy;
  --   gamma    value for the accessibility constant;
  --   homcrd   if homogeneous coordinates need to be used;
  --   rabin    if the user should be prompted for the Rabinowitsch trick.

  -- ON RETURN :
  --   nbequ    number of equations in the systems in the homotopy;
  --   sols     start solutions in the homotopy.

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false );
  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false );
  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
               -- tpow : in natural32 := 2;
                homcrd,rabin : in boolean := false );

  -- DESCRIPTION :
  --   Prompts for a target system, a start system with start solutions.
  --   The target and start system are stored in the Homotopy package.
  --   A random gamma constant will be generated by default.

  -- ON ENTRY :
  --   tpow     power of the continuation parameter
  --            in the artificial parameter homotopy.
  --   homcrd   if homogeneous coordinates need to be used;
  --   rabin    if the user should be prompted for the Rabinowitsch trick.

  -- ON RETURN :
  --   nbequ    number of equations in the systems in the homotopy;
  --   sols     start solutions in the homotopy.

  procedure Standard_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Parameter_Reader
              ( nbequ,nbvar,idxpar : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a natural parameter homotopy and start solutions
  --   in double, double double, or quad double precision.
  --   Assumes there is only one natural parameter.

  -- ON RETURN :
  --   nbequ    number of equations;
  --   nbvar    number of variables;
  --   idxpar   index of the variable which is the natural parameter;
  --   sols     start solutions in the homotopy.

  procedure Standard_Series_Newton
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out Standard_Complex_Series_Vectors.Vector );
  procedure Standard_Series_Newton
              ( sol : in Standard_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out Standard_Complex_Series_Vectors.Vector );
  procedure DoblDobl_Series_Newton
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out DoblDobl_Complex_Series_Vectors.Vector );
  procedure DoblDobl_Series_Newton
              ( sol : in DoblDobl_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out DoblDobl_Complex_Series_Vectors.Vector );
  procedure QuadDobl_Series_Newton
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out QuadDobl_Complex_Series_Vectors.Vector );
  procedure QuadDobl_Series_Newton
              ( sol : in QuadDobl_Complex_Solutions.Solution;
                idx,nbequ : in integer32; nbterms,nbiters : in natural32;
                srv,eva : out QuadDobl_Complex_Series_Vectors.Vector );

  -- DESCRIPTION :
  --   Given a start solution in the stored homotopy,
  --   runs Newton's method to compute a series solution,
  --   in double, double double, or quad double precision.

  -- ON ENTRY :
  --   idx      index of the parameter in the series, relative to
  --            the variables in the homotopy, in an arficial parameter
  --            homotopy for square systems, idx = nbequ + 1;
  --   nbequ    number of equations in the homotopy;
  --   nbterms  number of the terms in the series;
  --   nbiters  maximum number of iterations in Newton's method.

  -- ON RETURN :
  --   srv      series approximation for the solution;
  --   eva      evaluated series.

end Homotopy_Series_Readers;
