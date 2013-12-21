with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Symbol_Table;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Monomial_Maps;            use Standard_Monomial_Maps;

package Standard_Monomial_Maps_io is

-- DESCRIPTION :
--   This package provides input and output routines for monomial maps.
--   To enable the reading of a monomial map as a binomial system,
--   the d parameters are always written first, as a header.
--   This makes the input of a monomial map straightforward
--   via the reading of a binomial system.

-- MANAGEMENT OF SYMBOL TABLE :

  procedure Insert_Parameter_Symbols
               ( d : in natural32;
                 symbols : in Symbol_Table.Array_of_Symbols );

  -- DESCRIPTION :
  --   Clears the symbol table and inserts symbols for d parameters
  --   before the symbols in the given array of symbols.

  procedure Insert_Parameter_Symbols ( d : in natural32 );

  -- DESCRIPTION :
  --   Inserts d parameters into the symbol table.

  procedure Check_Parameter_Symbols ( n,d : in natural32 );

  -- DESCRIPTION :
  --   Checks whether the symbol table contains the right number of symbols
  --   for a d-dimensional monomial map in n-space.

  -- REQUIRED :
  --   The symbol table contains at least n symbols.

  procedure Remove_Parameter_Symbols ( n : in natural32 );

  -- DESCRIPTION :
  --   Removes the first d parameter symbols, d = Symbol_Table.number - n,
  --   and restores the symbol table to the first original n symbols.
  --   This procedure is called by the Show_Ideal procedure.

-- INPUT ROUTINES :

  procedure get ( map : out Link_to_Monomial_Map );
  procedure get ( file : in file_type; map : out Link_to_Monomial_Map );

  -- DESCRIPTION :
  --   Reads a monomial map from standard input or from file.
  --   The symbol table gets cleared before and after the reading.

  procedure get ( maps : out Monomial_Map_Array );
  procedure get ( file : in file_type; maps : out Monomial_Map_Array );

  -- DESCRIPTION :
  --   Reads the maps from standard input or from file.

  procedure get ( maps : out Link_to_Monomial_Map_Array );
  procedure get ( file : in file_type; maps : out Link_to_Monomial_Map_Array );

  -- DESCRIPTION :
  --   Reads first from standard input or from file the number of maps,
  --   followed by exactly as many maps as that number.

  procedure get ( maps : out Monomial_Map_List );
  procedure get ( file : in file_type; maps : out Monomial_Map_List );

  -- DESCRIPTION :
  --   Reads first from standard input or from file the number of maps,
  --   followed by exactly as many maps as that number.

-- OUTPUT ROUTINES :

  procedure put ( map : in Monomial_Map );
  procedure put ( file : in file_type; map : in Monomial_Map );

  -- DESCRIPTION :
  --   Writes the map to standard output or to file.

  -- REQUIRED : Insert_Parameter_Symbols(map.d) has been done.

  procedure put ( map : in Link_to_Monomial_Map );
  procedure put ( file : in file_type; map : in Link_to_Monomial_Map );

  -- DESCRIPTION :
  --   Writes the map to standard output or to file.

  procedure put ( maps : in Monomial_Map_Array );
  procedure put ( file : in file_type; maps : in Monomial_Map_Array );

  -- DESCRIPTION :
  --   Writes the monomial maps to standard output or to file.

  -- CAUTION : for maps of different dimension!

  procedure put ( maps : in Monomial_Map_List );
  procedure put ( file : in file_type; maps : in Monomial_Map_List );

  -- DESCRIPTION :
  --   Writes the monomial maps to standard output or to file.

  procedure put ( maps : in Array_of_Monomial_Map_Lists );
  procedure put ( file : in file_type;
                  maps : in Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Writes the monomial maps to standard output or to file.

-- IDEAL REPRESENTATIONS :

  procedure Show_Ideal ( p : in Laur_Sys; map : in Monomial_Map );
  procedure Show_Ideal ( file : in file_type;
                         p : in Laur_Sys; map : in Monomial_Map );

  -- DESCRIPTION :
  --   Shows the ideal defined by the map.

  procedure Show_Ideals ( p : in Laur_Sys; maps : in Monomial_Map_List );
  procedure Show_Ideals ( file : in file_type;
                          p : in Laur_Sys; maps : in Monomial_Map_List );

  -- DESCRIPTION :
  --   Shows the ideals defined by the maps in the list.

  procedure Show_Ideals ( p : in Laur_Sys;
                          maps : in Array_of_Monomial_Map_Lists );
  procedure Show_Ideals ( file : in file_type;
                          p : in Laur_Sys;
                          maps : in Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Shows the ideals defined by the maps in the array of lists.

-- APPEND MAPS TO A FILE :

  procedure Append ( name : in string; maps : in Monomial_Map_List );
  procedure Append ( name : in string;
                     maps : in Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Appends the solution maps to the file with the given name.

-- STATISTICS :

  procedure Write_Lengths ( c : in Array_of_Monomial_Map_Lists );
  procedure Write_Lengths 
              ( file : in file_type; c : in Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Writes the number of all components over all dimensions.

  procedure Show_Degrees ( maps : in Monomial_Map_List );
  procedure Show_Degrees ( maps : in Monomial_Map_List; sum : out natural32 );
  procedure Show_Degrees ( file : in file_type;
                           maps : in Monomial_Map_List );
  procedure Show_Degrees ( file : in file_type;
                           maps : in Monomial_Map_List; sum : out natural32 );
  procedure Show_Degrees ( p : in Laur_Sys;
                           maps : in Monomial_Map_List );
  procedure Show_Degrees ( p : in Laur_Sys;
                           maps : in Monomial_Map_List; sum : out natural32 );
  procedure Show_Degrees ( file : in file_type; p : in Laur_Sys;
                           maps : in Monomial_Map_List );
  procedure Show_Degrees ( file : in file_type; p : in Laur_Sys;
                           maps : in Monomial_Map_List; sum : out natural32 );

  procedure Show_Degrees ( maps : in Array_of_Monomial_Map_Lists );
  procedure Show_Degrees ( file : in file_type;
                           maps : in Array_of_Monomial_Map_Lists );
  procedure Show_Degrees ( p : in Laur_Sys;
                           maps : in Array_of_Monomial_Map_Lists );
  procedure Show_Degrees ( file : in file_type; p : in Laur_Sys;
                           maps : in Array_of_Monomial_Map_Lists );

  -- DESCRIPTION :
  --   Shows the degrees of the maps and if optionally, if p is provided,
  --   also shows the ideal representations of the maps.
  --   The sum on return is the sum of the degrees.

-- READ SYSTEM AND MAPS :

  procedure Read_System_and_Maps
               ( p : out Link_to_Laur_Sys;
                 maps : out Monomial_Map_List );
  procedure Read_System_and_Maps
               ( file : in file_type; p : out Link_to_Laur_Sys;
                 maps : out Monomial_Map_List );

  -- DESCRIPTION :
  --   Prompts a user for a file when the file is not provided
  --   and then reads a system and its solution maps from file.

end Standard_Monomial_Maps_io;
