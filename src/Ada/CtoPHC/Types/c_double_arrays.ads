with Interfaces.C;
with Interfaces.C.Pointers;

package C_Double_Arrays is

-- DESCRIPTION :
--   This package defines the type "C_Double_Array" to work with arrays
--   of C doubles and instantiates the C pointers packages.
--   The "Pointer_to_C_Double_Array" is needed to make matrices as
--   arrays of arrays of doubles.

  type C_Double_Array is
    array ( Interfaces.C.size_t range <> ) of aliased Interfaces.C.double;

  type Pointer_to_C_Double_Array is access C_Double_Array;

  package C_DblArrs is
    new Interfaces.C.Pointers(Interfaces.C.size_t,
                              Interfaces.C.double,
                              C_Double_Array,0.0);

  function Concat ( a,b : C_Double_Array ) return C_Double_Array;

  -- DESCRIPTION :
  --   Concatenates the array b after the array a.

end C_Double_Arrays;
