with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

generic

  type Item is private;

package Generic_Lists is

-- DESCRIPTION :
--   This generic package allows to implement lists of items.

  type List is private;

  Null_List : constant List;

  Overflow     : exception;
  List_Is_Null : exception;

-- CONSTRUCTORS :

  function Create ( i : Item ) return List;

  -- DESCRIPTION :
  --   Returns a list with i as its only item.

  procedure Construct ( i : in Item; L : in out List );

  -- DESCRIPTION :
  --   Adds the item i to the front of the list l.

  procedure Set_Head ( L : in out List; i : in Item );

  -- DESCRIPTION :
  --   Sets the first element in the list to item i.

  -- REQUIRED : not Is_Null(l).

  procedure Swap_Tail ( L1,L2 : in out List );

  -- DESCRIPTION :
  --   Swaps the tail of list l1 with the list l2.

  procedure Append ( first,last : in out List; i : in Item );

  -- DESCRIPTION :
  --   Appends the item i to the list, where first points to the first
  --   element and last to its last element.

  procedure Concat ( first,last : in out List; L : in List );

  -- DESCRIPTION :
  --   Concatenates the list l to the list first, where last points to
  --   the last element of the list.

  procedure Copy ( L1 : in List; L2 : in out List );

  -- DESCRIPTION :
  --   Makes a copy from the list L1 to the list L2.

-- SELECTORS :

  function Is_Equal ( L1,L2 : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both lists are equal.

  function Length_Of ( L : List ) return natural32;

  -- DESCRIPTION :
  --   Returns the length of the list, i.e.: the number of elements.

  function Is_Null ( L : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the list is empty, false otherwise.

  function Head_Of ( L : List) return Item;

  -- DESCRIPTION :
  --   Returns the first element in the list.

  -- REQUIRED : not Is_Null(l).

  function Tail_Of ( L : List ) return List;

  -- DESCRIPTION :
  --   Returns the tail of the list L.

  -- REQUIRED : not Is_Null(L).

-- DESTRUCTOR :

  procedure Clear ( L : in out List );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the elements in the list.
    
private

  type Node;
  type List is access Node;
  Null_List : constant List := null;

end Generic_Lists;
