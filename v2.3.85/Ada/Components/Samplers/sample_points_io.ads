with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Sample_Points;                     use Sample_Points;

package Sample_Points_io is

-- DESCRIPTION :
--   This package provides input-output operations on points sampled
--   from a component of solutions to a polynomial system.
--   The expected format on input is the same as the format on output.

  procedure get ( n,k : in integer32; sample : out Standard_Sample );
  procedure get ( file : in file_type; 
                  n,k : in integer32; sample : out Standard_Sample );
  procedure get ( n,k : in integer32; sample : out Multprec_Sample );
  procedure get ( file : in file_type; 
                  n,k : in integer32; sample : out Multprec_Sample );

  -- DESCRIPTION :
  --   Reads the solution and the hyperplane sections from file or
  --   standard input, where n is the length of the solution vector
  --   and k the number of hyperplane sections.

  procedure put ( sample : in Standard_Sample );
  procedure put ( file : in file_type; sample : in Standard_Sample );
  procedure put ( sample : in Multprec_Sample );
  procedure put ( file : in file_type; sample : in Multprec_Sample );

  -- DESCRIPTION :
  --   Writes the sample on file or on standard input.

end Sample_Points_io;
