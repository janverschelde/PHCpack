with Interfaces.C;
with Interfaces.C.Pointers;

package C_Integer_Arrays is

-- DESCRIPTION :
--   This package defines the type "C_Integer_Array" to work with arrays
--   of C integers and instantiates the C pointers packages.

  type C_Integer_Array is 
    array ( Interfaces.C.size_T range <> ) of aliased Interfaces.C.int;

  package C_intarrs is
    new Interfaces.C.Pointers(Interfaces.C.size_T,
                              Interfaces.C.int,
                              C_Integer_Array,0);

end C_Integer_Arrays;
