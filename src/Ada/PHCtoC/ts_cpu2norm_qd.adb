with text_io;                             use text_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Mathematical_Functions;     use QuadDobl_Mathematical_Functions;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;       use QuadDobl_Complex_Vector_Norms;
with Interfaces.C;
with C_Double_Arrays;                     use C_Double_Arrays;
with C_to_Ada_Arrays;

procedure ts_cpu2norm_qd is

-- DESCRIPTION :
--   Tests the computation of the 2-norm of a vector by a C function,
--   in complex quad double arithmetic.
--   The C function is called "ts_cpu2norm_in_c" and must be seen as 
--   a preparation to calling the GPU.

  procedure QuadDobl_Complex_Vector_Norm_in_Ada ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   All computations are defined by Ada code.

    v : QuadDobl_Complex_Vectors.Vector(1..dim);
    d_dim : double_float := double_float(dim);
    dd_dim : double_double := Double_Double_Numbers.Create(d_dim);
    qd_dim : quad_double := Quad_Double_Numbers.Create(dd_dim);
    nrm : quad_double;

  begin
    v := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    nrm := Norm2(v);
    put("The 2-norm of v : "); put(nrm); new_line;
    put("SQRT(dimension) : ");
    put(SQRT(qd_dim)); new_line;
  end QuadDobl_Complex_Vector_Norm_in_Ada;

  procedure QuadDobl_Complex_Vector_Norm_to_C ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   The vector is passed to C for the norm computation.

    v : QuadDobl_Complex_Vectors.Vector(1..dim);
    return_of_call : integer32;

    function normC ( n : integer32;        -- n is the dimension
                     x : C_Double_Array;   -- contains 2*n doubles
                     y : C_Double_Array )  -- on return is y(0) 
                   return integer32;
    pragma import(C, normC, "cpu2norm_qd_in_c");

  begin
    v := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    declare
      v_a : constant C_Double_Array := C_to_Ada_Arrays.Convert(v);
      res : C_Double_Array(0..3);
      hihi,lohi,hilo,lolo : double_float;
      result : quad_double;
    begin
      res(0) := 0.0;
      res(1) := 0.0;
      res(2) := 0.0;
      res(3) := 0.0;
      return_of_call := normC(dim,v_a,res);
      hihi := double_float(res(0));
      lohi := double_float(res(1));
      hilo := double_float(res(2));
      lolo := double_float(res(3));
      result := Quad_Double_Numbers.Create(hihi,lohi,hilo,lolo);
      put("The result of C : "); put(result); new_line;
      put("The sqrt of it  : "); put(SQRT(result)); new_line;
    end;
  end QuadDobl_Complex_Vector_Norm_to_C;

  procedure Main is

    n : integer32 := 0;

  begin
    new_line;
    put_line("Accelerating the 2-norm computation of a complex vector ...");
    new_line;
    put("Give the dimension : "); get(n);
    put("-> the dimension is "); put(n,1); put_line(".");
    new_line;
    QuadDobl_Complex_Vector_Norm_in_Ada(n);
    QuadDobl_Complex_Vector_Norm_to_C(n);
  end Main;

begin
  Main;
end ts_cpu2norm_qd;
