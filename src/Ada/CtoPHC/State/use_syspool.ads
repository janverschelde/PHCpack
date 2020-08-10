with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;

function use_syspool ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32;

-- DESCRIPTION :
--   Provides a gateway to the systems pool.

-- ON ENTRY :
--   job     =  0 : initialize standard system pool with n = a[0];
--           =  1 : get the size of the system pool returned in a[0];
--           =  2 : reads a system p iand with k = a[0], the data of p
--                  for the k-th system in the pool is created;
--           =  3 : writes the k-th system, with k = a[0];
--           =  4 : takes system of the standard container to make k-th 
--                  system in the standard system pool, k = a[0];
--           =  5 : refines a solution using the k-th system in the container,
--                  k = a[0], n = a[1], m = b[0] and c contains the floating
--                  part of the solution;
--           =  6 : copies the k-th system in the pool to the standard system
--                  container, where k is given as a[0];
--           =  7 : copies the k-th system in the pool to the dobldobl system
--                  container, where k is given as a[0];
--           =  8 : copies the k-th system in the pool to the quaddobl system
--                  container, where k is given as a[0];
--           =  9 : the size of the dobldobl systems pool is returned in a[0];
--           = 10 : the size of the dobldobl systems pool is returned in a[0];
--           = 11 : initialize dobldobl system pool with n = a[0];
--           = 12 : initialize quaddobl system pool with n = a[0];
--           = 13 : clears the standard system pool;
--           = 14 : clears the dobldobl system pool;
--           = 15 : clears the quaddobl system pool;
--           = 16 : takes system of the dobldobl container to make k-th 
--                  system in the dobldobl system pool, k = a[0];
--           = 17 : takes system of the quaddobl container to make k-th 
--                  system in the quaddobl system pool, k = a[0].
-- 
--   a        memory allocated a natural number, either the size
--            of the systems pool or the index of a system in the pool;
--   b        used to pass and return multiplicity of a solution;
--   c        for the floating-point numbers in a solution.

-- ON RETURN :
--   0 if the operation was successful, otherwise something went wrong,
--   or job not in the right range.
