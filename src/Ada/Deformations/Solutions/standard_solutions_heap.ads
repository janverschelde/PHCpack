with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with generic_lists;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Solutions_Heap is

-- DESCRIPTION :
--   To efficiently determine whether there are clusters in long
--   lists of solutions, a heap is constructed.
--   Items on the heap are the following four items:
--   (1) two weighted values computed with random
--       weight vectors on the solution vectors, 
--   (2) the index number of the solution, and
--       the link to the solution.
--   Only in case the two weight values are the same, 
--   are the coordinates of two solutions compared.
--   The heap is implemented via a vector,
--   but because of the long length of solution lists,
--   this vector is implemented as a bucketed vector.

-- THE HEAP DATA STRUCTURE :

  type Heap_Item is record
    first,second : double_float;
    idx : integer32;
    ls : Link_to_Solution;
  end record;

  type Heap_Items is array ( integer32 range <> ) of Heap_Item;
  type Link_to_Heap_Items is access Heap_Items;

  package Buckets is new Generic_Lists(Link_to_Heap_Items);
  type Bucket_Vector is new Buckets.List;

  bucketsize : constant integer32 := 1024;

 -- The original Heap type is kept for reference.
 -- type Heap ( nbr : integer32 ) is record
 --   bottom : integer32 := -1; -- next free spot on heap
 --   values : Heap_Items(0..nbr);
 -- end record;

  type Heap is record
    bottom : integer32 := -1; -- next free spot on heap
    values,values_last : Bucket_Vector;
  end record;

-- WEIGHING SOLUTION VECTORS :

  function Random_Weight_Vector
             ( nbr : in integer32 )
             return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of nbr random floating-point numbers,
  --   in the interval [-1, +1].

  function Weight ( v : Standard_Complex_Vectors.Vector;
                    w : Standard_Floating_Vectors.Vector )
                  return double_float;

  -- DESCRIPTION :
  --   Returns the value of the vector v for the weights w,
  --   obtained by multiplying the real and imaginary parts
  --   of the components of v with the weights in w.

-- OPERATIONS ON BUCKET VECTORS :

  function Size ( b : Bucket_Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of elements in the bucket vector b.

  function New_Heap_Items return Link_to_Heap_Items;

  -- DESCRIPTION :
  --   Returns an allocated vector of heap items,
  --   with everything initialized to zero.

  procedure Assign ( b,b_last : in out Bucket_Vector;
                     i : in integer32; v : in Heap_Item );

  -- DESCRIPTION :
  --   Assigns to b(i) the value v.
  --   The b_last points to the last array in b.
  --   Allocates as many new buckets as needed.

  function Retrieve ( b : Bucket_Vector; i : integer32 ) return Heap_Item;

  -- DESCRIPTION :
  --   Returns the item at position i in the vector b.
  --   If i exceeds the size of b, then the idx of the returned item is -1.

  procedure Clear ( h : in out Link_to_Heap_Items );

  -- DESCRIPTION :
  --   Deallocates the space for the heap items in h.
 
  procedure Clear ( b : in out Bucket_Vector );

  -- DESCRIPTION :
  --   Deallocates all space allocated for the bucket vector b.

-- HEAP OPERATIONS :

  procedure Swap_from_Bottom ( h : in out Heap; p : in integer32 );

  -- DESCRIPTION :
  --   Swaps items of the heap h, starting at position p,
  --   as long as the item at position p has a larger value
  --   than the value of its parent.

  -- REQUIRED : p is a valid index for h.values,
  --   initially the position of an item pushed into the heap h.

  procedure Push ( h : in out Heap; i : in Heap_Item );

  -- DESCRIPTION :
  --   Pushes the item i into the heap h.

  function Max_Child ( h : Heap; p : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns -1 or the index of the largest child of the node
  --   at position p in the heap h.
  --   The value of h.bottom is the position of the empty spot
  --   after the last valid item on the heap.

  procedure Swap_from_Top ( h : in out Heap; p : in integer32 );

  -- DESCRIPTION :
  --   Swaps items of the heap h, starting at position p,
  --   as long as the item at position p has a smaller value
  --   than the value of any of its children.

  procedure Pop ( h : in out Heap );

  -- DESCRIPTION :
  --   Removes the top item at h.values(0).

  procedure Push ( h : in out Heap;
                   w1,w2 : in double_float; i : in integer32;
                   ls : in Link_to_Solution );

  -- DESCRIPTION :
  --   Makes an item with weight w and index i
  --   and pushes that item into the heap h.

  procedure Clear ( h : in out Heap );

  -- DESCRIPTION :
  --   Clears all space allocated to the heap h.

-- CLUSTER REPORT :

  procedure Count_Clusters
              ( h : in out Heap; tol : in double_float;
                cnt : out natural32; verbose : in boolean := true );

  -- DESCRIPTION :
  --   Computes the count of clustered solutions,
  --   if verbose, then writes extra output to screen for each solution.

  -- ON ENTRY :
  --   h        a heap constructed for a solution list;
  --   tol      tolerance to decide whether two doubles are close enough;
  --   verbose  if true, then extra output is written to screen.

  -- ON RETURN :
  --   cnt      number of clustered solutions.

end Standard_Solutions_Heap;
