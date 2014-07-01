with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Mathematical_Functions;     use DoblDobl_Mathematical_Functions;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Vector_Norms;       use DoblDobl_Complex_Vector_Norms;
with Interfaces.C;
with C_Double_Arrays;                     use C_Double_Arrays;
with C_to_Ada_Arrays;

procedure ts_cpu2norm_dd is

-- DESCRIPTION :
--   Tests the computation of the 2-norm of a vector by a C function,
--   in complex double double arithmetic.
--   The C function is called "ts_cpu2norm_in_c" and must be seen as 
--   a preparation to calling the GPU.

  procedure DoblDobl_Complex_Vector_Norm_in_Ada ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   All computations are defined by Ada code.

    v : DoblDobl_Complex_Vectors.Vector(1..dim);
    d_dim : double_float := double_float(dim);
    dd_dim : double_double := Double_Double_Numbers.Create(d_dim);
    nrm : double_double;

  begin
    v := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    nrm := Norm2(v);
    put("The 2-norm of v : "); put(nrm); new_line;
    put("SQRT(dimension) : ");
    put(SQRT(dd_dim)); new_line;
  end DoblDobl_Complex_Vector_Norm_in_Ada;

  procedure DoblDobl_Complex_Vector_Norm_to_C ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   The vector is passed to C for the norm computation.

    v : DoblDobl_Complex_Vectors.Vector(1..dim);
    return_of_call : integer32;

    function normC ( n : integer32;        -- n is the dimension
                     x : C_Double_Array;   -- contains 2*n doubles
                     y : C_Double_Array )  -- on return is y(0) 
                   return integer32;
    pragma import(C, normC, "cpu2norm_dd_in_c");

  begin
    v := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    declare
      v_a : constant C_Double_Array := C_to_Ada_Arrays.Convert(v);
      res : C_Double_Array(0..1);
      lo,hi : double_float;
      result : double_double;
    begin
      res(0) := 0.0;
      res(1) := 0.0;
      return_of_call := normC(dim,v_a,res);
      hi := double_float(res(0));
      lo := double_float(res(1));
      result := Double_Double_Numbers.Create(hi,lo);
      put("The result of C : "); put(result); new_line;
      put("The sqrt of it  : "); put(SQRT(result)); new_line;
    end;
  end DoblDobl_Complex_Vector_Norm_to_C;

  procedure Main is

    n : integer32 := 0;

  begin
    new_line;
    put_line("Accelerating the 2-norm computation of a complex vector ...");
    new_line;
    put("Give the dimension : "); get(n);
    put("-> the dimension is "); put(n,1); put_line(".");
    new_line;
    DoblDobl_Complex_Vector_Norm_in_Ada(n);
    DoblDobl_Complex_Vector_Norm_to_C(n);
  end Main;

begin
  Main;
end ts_cpu2norm_dd;
