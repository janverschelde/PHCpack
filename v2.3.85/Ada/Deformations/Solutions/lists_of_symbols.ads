with generic_lists;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Symbol_Table;                       use Symbol_Table;

package Lists_of_Symbols is

-- DESCRIPTION :
--   This packages defines a list of symbols for use in the packages
--   to convert solution lists to Maple or Python dictionary formats.

  package Symbols_Lists is new Generic_Lists(Symbol);
  type Symbol_List is new Symbols_Lists.List;

  procedure Create_Symbol_Table ( L : in Symbol_List );

  -- DESCRIPTION :
  --   Creates the symbol table using the symbols in L.

  -- REQUIRED : Symbol_Table.Number = 0.

  function Equal ( sb : Symbol; s : string ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the symbol in sb matches the string s.

  function Classify_Symbol ( sb : Symbol ) return natural32;

  -- DESCRIPTION :
  --   Returns a number between 0 and 5 with the following meaning:
  --     0 if the symbol sb is not special
  --     1 if sb = "time", 2 if sb = "multiplicity"
  --     3 if sb = "err", 4 if sb = "rco", 5 if sb = "res".

  procedure Clear ( L : in out Symbol_List );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the list L.

end Lists_of_Symbols;
