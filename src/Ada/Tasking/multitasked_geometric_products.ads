with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;              use Quad_Double_Numbers;

package Multitasked_Geometric_Products is

-- DESCRIPTION :
--   A geometric inner product is the sum of the componentwise
--   product of two geometric sequences.

  function Inner_Product ( dim : integer32; rtx,rty : double_float ) 
                         return double_float;
  function Inner_Product ( dim : integer32; rtx,rty : double_double ) 
                         return double_double;
  function Inner_Product ( dim : integer32; rtx,rty : quad_double ) 
                         return quad_double;

  -- DESCRIPTION :
  --   Returns the inner product of two geometric sequences of length dim,
  --   where the ratios of the sequences are rtx and rty.

  function Inner_Product ( idx_start,idx_end : integer32;
                           rtx,rty : double_float ) 
                         return double_float;
  function Inner_Product ( idx_start,idx_end : integer32;
                           rtx,rty : double_double ) 
                         return double_double;
  function Inner_Product ( idx_start,idx_end : integer32;
                           rtx,rty : quad_double ) 
                         return quad_double;

  -- DESCRIPTION :
  --   Returns the inner product of two geometric sequences,
  --   where the ratios of the sequences are rtx and rty,
  --   starting at the index idx_start and ending at idx_end.

  procedure Double_Sequential_Run;

  -- DESCRIPTION :
  --   Runs an instance which takes about one second on one thread,
  --   in double precision.

  procedure Double_Double_Sequential_Run;

  -- DESCRIPTION :
  --   Runs an instance in double double precision.

  procedure Quad_Double_Sequential_Run;

  -- DESCRIPTION :
  --   Runs an instance in quad double precision.

  procedure Double_Parallel_Run ( p : in integer32 );

  -- DESCRIPTION :
  --   Runs with p threads, in double precision.

  procedure Double_Double_Parallel_Run ( p : in integer32 );

  -- DESCRIPTION :
  --   Runs with p threads, in double double precision.

  procedure Quad_Double_Parallel_Run ( p : in integer32 );

  -- DESCRIPTION :
  --   Runs with p threads, in quad double precision.

end Multitasked_Geometric_Products;
