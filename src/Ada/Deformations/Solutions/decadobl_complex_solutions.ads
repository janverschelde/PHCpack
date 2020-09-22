with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Solutions;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Generic_Lists;
with DecaDobl_Complex_Numbers;           use DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Vectors;           use DecaDobl_Complex_Vectors;
with Multprec_Complex_Solutions;

package DecaDobl_Complex_Solutions is

-- DESCRIPTION :
--   This package provides an abstraction of a list and an array
--   of solutions as vectors of deca double complex numbers.

-- DATA STRUCTURES :
 
  type Solution ( n : integer32 ) is record
    t : Complex_Number;         -- continuation parameter t
    m : integer32;              -- multiplicity of the solution; or flag
    v : Vector(1..n);           -- the solution
    err : deca_double;          -- error = |correction| from Newton
    rco : deca_double;          -- inverse of condition number
    res : deca_double;          -- norm of residual vector
  end record;

  type Link_to_Solution is access Solution;

  package List_of_Solutions is new Generic_Lists(Link_to_Solution);
  type Solution_List is new List_of_Solutions.List;

  type Solution_Array is array ( integer32 range <> ) of Link_to_Solution;

  type Array_of_Solution_Lists is 
    array ( integer32 range <> ) of Solution_List;
  type Link_to_Array_of_Solution_Lists is access Array_of_Solution_Lists;

-- CREATORS :

  function Create ( sl : Solution_List ) return Solution_Array;
  function Create ( sa : Solution_Array ) return Solution_List;

  -- DESCRIPTION :
  --   Allows the transition from a list to an array and vice versa.

  function Create ( s : Standard_Complex_Solutions.Solution )
                  return DecaDobl_Complex_Solutions.Solution;
  function Create ( s : Multprec_Complex_Solutions.Solution )
                  return DecaDobl_Complex_Solutions.Solution;

  -- DESCRIPTION :
  --   Converts the solution as a standard/multprec complex vector
  --   into a representation that uses deca double numbers.

  function Create ( L : Standard_Complex_Solutions.Solution_List )
                  return DecaDobl_Complex_Solutions.Solution_List;
  function Create ( L : Multprec_Complex_Solutions.Solution_List )
                  return DecaDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   Converts the list of solutions as standard/multprec vectors
  --   into a representation that uses deca double numbers.

-- COMPARISON and COPYING :

  function Equal ( s1,s2 : Solution; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if for each component 
  --   |s1.v(i)-s2.v(i)|/|s1.v(i)| < tol, for i in s1.v'range.

  function Equal ( s1,s2 : Solution_List; tol : double_float )
                 return boolean;
  function Equal ( s1,s2 : Solution_Array; tol : double_float )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true if both lists of arrays are equal to each other, upon
  --   the given tolerance for the relative error.

  procedure Equals ( sols : in out Solution_List; flag : in integer32;
                     tol : in double_float; same : out boolean );
  -- DESCRIPTION :
  --   The solutions that are equal to each other are marked with a flag.

  procedure Equals ( sa : in Solution_Array; x : in Vector; i : in integer32;
                     tol : in double_float; j : in out integer32 );
  -- DESCRIPTION :
  --   Compares the first i-1 vectors in sa with x.

  -- ON ENTRY :
  --   sa        a solution array, containing at least i-1 elements;
  --   x         a vector;
  --   i         an index, normally the entry of x in sa;
  --   tol       tolerance for relative error on two vectors;
  --   j         must be equal to sa'first.

  -- ON RETURN :
  --   j         the entry for which sa(j) equals x.

  procedure Copy ( s1 : in Solution; s2 : in out Solution );

  -- DESCRIPTION :
  --   Makes a deep copy of the solution s1 into the solution s2.

  procedure Copy ( s1 : in Solution_List; s2 : in out Solution_List );
  procedure Copy ( s1 : in Solution_Array; s2 : in out Solution_Array );

  -- DESCRIPTION :
  --   Makes a deep copy of the list or the array of solutions.

-- SELECTORS :

  function Number ( sols : Solution_List; flag : integer32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of solutions in the list with a multiplicity
  --   equal to flag.

  function Is_In ( sols : Solution_List; s : Solution; tol : double_float ) 
		 return boolean;
  function Is_In ( sa : Solution_Array; s : Solution; tol : double_float ) 
		 return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution s belongs to the list or to the array.

  function Retrieve ( sols : Solution_List; p : natural32 )
                    return Link_to_Solution;

  -- DESCRIPTION :
  --   Returns the solution at position p in the list sols,
  --   or otherwise a null pointer if p is out of range.

-- CONSTRUCTORS :

  procedure Push ( s0 : in Solution_List; s1 : in out Solution_List );

  -- DESCRIPTION :
  --   Pushes all solutions of s0 in front of the list s1.
  --   Note: only the pointers of solutions in s0 are pushed.

  procedure Append ( first,last : in out Solution_List;
                     ls : in Link_to_Solution );

  -- DESCRIPTION :
  --   The link to the solution is appended to the list first,
  --   where last is a pointer to the last element of the list.

  procedure Append ( first,last : in out Solution_List; s : in Solution );

  -- DESCRIPTION :
  --   The solution sol is appended to the list first, where
  --   last is a pointer to the last element of the list.

  procedure Add ( sols : in out Solution_List; s : in Solution );

  -- DESCRIPTION :
  --   The solution sol is appended to the list sols.

  procedure Add ( sols : in out Solution_List; s : in Solution;
                  tol : in double_float; other : out natural32 );

  -- DESCRIPTION :
  --   Append the solution to the list, if it does not already belong to it.

-- MODIFIERS :

  procedure Change ( sols : in out Solution_List; pos : in natural32;
                     s : in Solution; tol : in double_float; 
                     other : out natural32 );

  -- DESCRIPTION :
  --   Changes the solution at the given position into s, if the solution
  --   does not already occur.

  -- REQUIRED : pos <= Length_Of(sols).

  procedure Set_Continuation_Parameter
                ( sols : in out Solution_List; t : in Complex_Number );

  -- DESCRIPTION :
  --   All solutions in the list will be given the continuation parameter t.

  procedure Change_Multiplicity
                ( sols : in out Solution_List; pos : in natural32;
                  m : in integer32 );

  -- DESCRIPTION :
  --   Changes the multiplicity of the solution with the given position
  --   into m.

  -- REQUIRED : pos <= Length_Of(sols).

  procedure Remove ( sols : in out Solution_List; pos : in natural32 );

  -- DESCRIPTION :
  --   Removes the solution with given position from the list.

  -- REQUIRED : pos <= Length_Of(sols).

  generic
    with function To_Be_Removed ( flag : in integer32 ) return boolean;
  procedure Delete ( sols : in out Solution_List );

  -- DESCRIPTION :
  --   Removes all solutions in s for which To_Be_Removed(s.m) holds.

  procedure Remove_All ( sols : in out Solution_List; flag : in integer32 );

  -- DESCRIPTION :
  --   All solutions with a multiplicity equal to flag are removed.

-- DESTRUCTORS :

  procedure Clear ( sa : in out Solution_Array );
  procedure Clear ( s : in out Solution );
  procedure Clear ( ls : in out Link_to_Solution );
  procedure Shallow_Clear ( sl : in out Solution_List );
  procedure    Deep_Clear ( sl : in out Solution_List );

  -- DESCRIPTION :
  --   Deallocation of the occupied memory space.
  --   A shallow clear only deallocates the pointers,
  --   so that the data may still be accessible by sharing,
  --   whereas a deep clear also makes the data inaccessible.

  procedure Clear ( s : in out Array_of_Solution_Lists );
  procedure Clear ( s : in out Link_to_Array_of_Solution_Lists );

  -- DESCRIPTION :
  --   Deallocation of all (deep clear) the space occupied by s.

end DecaDobl_Complex_Solutions;
