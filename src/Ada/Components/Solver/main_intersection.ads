with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Main_Intersection is

-- DESCRIPTION :
--   Constructs witness set representations for the intersection of two
--   algebraic sets given in their witness set representations.

  procedure Read_Witness_Set
              ( w : in string; k : in natural32; p : out Link_to_Poly_Sys;
                sols : out Solution_List; dim : out natural32;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Reads witness set k from file with name in the string w.

  -- ON ENTRY :
  --   w        name of file, if empty, the user is asked for a name;
  --   k        number of witness set, used for prompting;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   p        polynomial equations defining the witness set;
  --   sols     points in the witness set;

  procedure Intersect_Witness_Sets
              ( file : in file_type; filename : in string;
                p1,p2 : in Poly_Sys; w1,w2 : in Solution_List;
                d1,d2 : in natural32; vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Intersects two witness sets of dimensions d1 and d2, d1 >= d2.

  -- ON ENTRY :
  --   file     log file to be opened for output;
  --   filename is the name of the output file, will be used to create
  --            files with witness sets of the intersection;
  --   p1,p2    polynomials defining the two solution components;
  --   w1,w2    witness points of the two solution components;
  --   d1,d2    dimensions of the two witness sets;
  --   vrblvl   is the verbose level.

  procedure Main ( witset_one,witset_two,logfile : in string;
                   vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies diagonal homotopies to intersect to witness sets.
  --   Defines phc -w.

  -- REQUIRED : the algebraic sets live in the same ambient space,
  --   if the file names are empty strings, then the user will be
  --   prompted to provide the file names.

  -- ON ENTRY :
  --   witset_one   witness set for the first algebraic set;
  --   witset_two   witness set for the second algebraic set;
  --   logfile      file name to write diagnostics on;
  --   vrblvl       is the verbose level.

  -- ON RETURN :
  --   logfile_w0, logfile_w1, etc. will contain the super witness sets
  --   for components of dimension 0, 1, etc. of the intersection of the
  --   two witness sets.

end Main_Intersection;
