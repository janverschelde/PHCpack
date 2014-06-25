with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package Symbol_Table is

-- DESCRIPTION :
--   This package provides management of a table of symbols, useful
--   in the input/output of polynomials in several variables.

  delimiter : constant character := ';';  -- separates polynomials
  type Power is ('*','^');

-- DATA STRUCTURES :

  subtype Symbol is string(1..80);
  type Array_of_Symbols is array ( integer32 range <> ) of Symbol;
  type Link_to_Array_of_Symbols is access Array_of_Symbols;

-- EXCEPTIONS :
    
    OVERFLOW_IN_THE_SYMBOL_TABLE : exception;
     -- occurs when a new symbol is added to a full symbol table

    INDEX_OUT_OF_RANGE : exception;
     -- occurs when a symbol is asked that is not in the range of the table

    INVALID_SYMBOL : exception;
     -- occurs when an invalid symbol is read, see "Is_Valid" below

-- CREATORS :

  procedure Init ( max : in natural32 );

  -- DESCRIPTION :
  --   A new symbol table is created with place for max symbols.

  function Standard_Symbols ( n : integer32 ) return Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns the array of symbols x1, x2, .., xn.

  procedure Init ( s : in Array_of_Symbols );

  -- DESCRIPTION :
  --   Initializes the symbol table with the symbols in s.
 
  procedure Enlarge ( max : in natural32 );

  -- DESCRIPTION :
  --   Enlarges the symbol table so that it can contain as many symbols
  --   as the previous maximum plus the new max.

  procedure Downsize ( n : in natural32 );

  -- DESCRIPTION :
  --   Reduces the maximal size of the symbol table 
  --   to hold n symbols less.  This is useful after
  --   removing variables from the symbol table.

  procedure Replace ( i : in natural32; sb : in Symbol );

  -- DESCRIPTION :
  --   Replaces the ith symbol by the given sb.

  function Create ( s : string ) return Symbol;

  -- DESCRIPTION :
  --   Returns the symbol created from the string.

-- CONSTRUCTORS :

  procedure Add_String ( s : in string );
  procedure Add_String ( s : in string; pos : out natural32 );

  -- DESCRIPTION :
  --   Adds the string s as a new symbol to the table.
  --   The optional "pos" on return equals the position of s.

  procedure Add ( sb : in Symbol );
  procedure Add ( sb : in Symbol; pos : out natural32 );

  -- DESCRIPTION :
  --   A new symbol is added to the symbol table;
  --   pos is the entry of the added symbol in the table.

  procedure Remove ( sb : in Symbol );
  procedure Remove ( i : in natural32 );

  -- DESCRIPTION :
  --   Removes the ith symbol sb from the symbol table.

-- SELECTORS :

  function Length ( s : Symbol ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of non blanks in the symbol.

  function Equal ( s1,s2 : Symbol ) return boolean;

  -- DESCRIPTION :
  --   Returns s1 = s2, ignoring trailing spaces.

  function Is_Valid ( sb : Symbol ) return boolean;

  -- DESCRIPTION :
  --   A symbol is valid if 
  --   (1) it does not start with a digit, and
  --   (2) does not contain any of the arithmetical operation symbols, and
  --   (3) does not contain the semi colon;
  --   (4) does not contain and round brackets.

  function Check_Symbol ( n : natural32; sb : Symbol ) return natural32;

  -- DESCRIPTION :
  --   Checks whether the symbol is a valid name for a variable.
  --   The exception "INVALID_SYMBOL" is raised if sb is invalid,
  --   otherwise, the "OVERFLOW_OF_UNKNOWNS" might be raised if
  --   it cannot be added to the symbol table (<= n elements).
  --   In case no exceptions occur, the number of the variable symbol
  --   in the table is returned.

  function "<" ( s1,s2 : Symbol ) return boolean;
  function ">" ( s1,s2 : Symbol ) return boolean;

  -- DESCRIPTION :
  --   The order relation on the symbols is defined by the order
  --   in which the symbols were added to the symbol table.

  function Maximal_Size return natural32;

  -- DESCRIPTION :
  --   Returns the maximal number of symbols the table can contain.

  function Number return natural32;
    
  -- DESCRIPTION :
  --   Returns the number of current symbols in the table.

  function Empty return boolean;

  -- DESCRIPTION :
  --   Returns true if the symbol table has not been initialized yet,
  --   or if a Clear has been done.

  function Get ( sb : Symbol ) return natural32;

  -- DESCRIPTION :
  --  The entry of the symbol in the table is returned.
  --  If the symbol does not occur in the table, then 0 is returned.

  function Get ( i : natural32 ) return Symbol;

  -- DESCRIPTION :
  --   The symbol corresponding with the ith unknown is returned.

  function Content return Array_of_Symbols;

  -- DESCRIPTION :
  --   Returns the content of the symbol table in the array of symbols.
  -- REQUIRED : Number > 0.

-- DESTRUCTORS :

  procedure Clear ( ls : in out Link_to_Array_of_Symbols );

  -- DESCRIPTION :
  --   Releases the memory occupied by the array of symbols.

  procedure Clear;
  
  -- DESCRIPTION :
  --   The allocated memory space is freed.

end Symbol_Table;
