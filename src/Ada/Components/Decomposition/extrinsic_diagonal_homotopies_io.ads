with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Permutations;                      use Permutations;
with Symbol_Table;                      use Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Extrinsic_Diagonal_Homotopies_io is

-- DESCRIPTION :
--   This package provides operations to double the symbol table,
---  when diagonal homotopies are used.

  procedure Write ( sa : in Array_of_Symbols );
  procedure Write ( file : in file_type; sa : in Array_of_Symbols ); 

  -- DESCRIPTION :
  --   Writes all symbols in the array to screen or to the file.

  function Add_Suffix ( sb : Symbol; c : character ) return Symbol;
  function Add_Suffix ( sa : Array_of_Symbols; c : character )
                      return Array_of_Symbols;

  -- DESCRIPTION :
  --   Adds the character c to the symbol sb or to the symbols in sa.

  function Suffix ( sb : Symbol ) return character;

  -- DESCRIPTION :
  --   Returns the last character of the symbol.

  function Remove_Suffix ( sb : Symbol ) return Symbol;

  -- DESCRIPTION :
  --   Removes the last character from the symbol.

  function Get_Symbols return Array_of_Symbols;
  function Get_Link_to_Symbols return Link_to_Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns the array with all symbols from the symbol table.

  function Suffixed_Symbols ( c : character ) return Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns an array with all symbols from the symbol table,
  --   appended with the suffix c.

  function Is_Embed_Symbol ( sb : Symbol ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the symbol is a symbol used in the embedding,
  --   i.e.: whether the first two characters are z.

  function Retrieve_Suffixed_Symbols
             ( n : in integer32; c : character ) return Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns an array of range 1..n with all symbols that end with c
  --   and for which Is_Embed_Symbol returns false.

  procedure Retrieve_Suffixed_Symbols
              ( n : in integer32; c : character;
                s : out Array_of_Symbols; p : out Permutation );

  -- DESCRIPTION :
  --   In addition to returning in s the symbols whose last character is c,
  --   this procedure also returns in p the corresponding indices where 
  --   the symbols were found in the symbol table.

  -- REQUIRED : s'range = p'range = 1..n.

  function Remove_Embed_Symbols
              ( sa : Array_of_Symbols ) return Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns the array with embedding symbols removed from sa.

  function Equals_mod_Suffix ( sb1,sb2 : Symbol ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the symbols are the same, disregarding the
  --   last character.

  function Look_for_Position 
              ( sa : Array_of_Symbols; sb : Symbol ) return integer32;

  -- DESCRIPTION :
  --   Returns 0 if the symbol sb does not occur in the array sa,
  --   otherwise the position of the symbol in the array is returned.

  function Search_Position 
              ( sa : Array_of_Symbols; sb : Symbol ) return integer32;

  -- DESCRIPTION :
  --   Returns 0 if the symbol sb does not occur in the array sa,
  --   otherwise the position of the symbol in the array is returned,
  --   suffixes are ignored.

  function Match_Symbols
              ( s1,s2 : in Array_of_Symbols ) return Permutation;

  -- DESCRIPTION :
  --   Matches the symbols in second array with those in the first.
  --   Creates a permutation, for every symbol in s2, the i-th position
  --   in the permutation equals the position of the same symbol in s1.

  function Matching_Permutation
              ( s1,s2 : Array_of_Symbols ) return Permutation;

  -- DESCRIPTION :
  --   Matches the symbols in second array with those in the first.
  --   Creates a permutation, for every symbol in s2, the i-th position
  --   in the permutation equals the position of the same symbol in s1.
  --   The suffixes (i.e.: last character of every symbol) are ignored.

  -- REQUIRED : s1'range = s2'range.

  function Combine_Permutations
              ( n,d : integer32; p1,p2,prm : Permutation ) return Permutation;

  -- DESCRIPTION :
  --   Returns one permutation of range 1..2*n+d, combining the information
  --   from the three permutations p1, p2, and prm.

  -- ON ENTRY :
  --   n        dimensions before the doubling and embedding;
  --   d        number of embed variables;
  --   p1       indices to the first group of symbols in the table;
  --   p1       indices to the second group of symbols in the table;
  --   prm      prm matches the two groups of symbols.

  procedure Permute ( prm : in Permutation; s : in out Array_of_Symbols );

  -- DESCRIPTION :
  --   Permutes the symbols in the array according the permutation.

  procedure Write_Symbol_Table;

  -- DESCRIPTION :
  --   Writes the current contents of the symbol table to screen.

  procedure Assign_Symbol_Table ( sa : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Assigns the symbol table with the symbols in the symbol array.

  procedure Assign_Symbol_Table ( s1,s2 : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Assigns the symbol table with the symbols in the symbol arrays.

  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List );
  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_Witness_Set
              ( file : in file_type; name : in string; d : in natural32;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the witness set of dimenension d defined by the system p
  --   and the solutions in sols to the file with the name filename,
  --   followed by the suffix _swd, where d is the value of d.
  --   Writes one line to the file, given as first argument.

end Extrinsic_Diagonal_Homotopies_io;
