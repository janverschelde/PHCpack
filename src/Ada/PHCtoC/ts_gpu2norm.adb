with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with DoblDobl_Mathematical_Functions;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with QuadDobl_Mathematical_Functions;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Standard_Complex_Vector_Norms;       use Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;       use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;       use QuadDobl_Complex_Vector_Norms;
with C_Double_Arrays;                     use C_Double_Arrays;
with C_to_Ada_Arrays;

procedure ts_gpu2norm is

-- DESCRIPTION :
--   Tests the computation of the 2-norm of a vector by the GPU,
--   via the C functions "gpu2norm_d", "gpu2norm_dd", and "gpu2norm_qd",
--   for double, double double, and quad double precision.

  procedure Standard_Complex_Vector_Norm_in_Ada ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   All computations are defined by Ada code.

    use Standard_Mathematical_Functions;

    v : Standard_Complex_Vectors.Vector(1..dim);
    nrm : double_float;

  begin
    v := Standard_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    nrm := Norm2(v);
    put("The 2-norm of v :"); put(nrm); new_line;
    put("SQRT(dimension) :");
    put(SQRT(double_float(dim))); new_line;
  end Standard_Complex_Vector_Norm_in_Ada;

  procedure DoblDobl_Complex_Vector_Norm_in_Ada ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   All computations are defined by Ada code.

    use DoblDobl_Mathematical_Functions;

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

  procedure QuadDobl_Complex_Vector_Norm_in_Ada ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   All computations are defined by Ada code.

    use QuadDobl_Mathematical_Functions;

    v : QuadDobl_Complex_Vectors.Vector(1..dim);
    d_dim : double_float := double_float(dim);
    dd_dim : quad_double := Quad_Double_Numbers.Create(d_dim);
    nrm : quad_double;

  begin
    v := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    nrm := Norm2(v);
    put("The 2-norm of v : "); put(nrm); new_line;
    put("SQRT(dimension) : ");
    put(SQRT(dd_dim)); new_line;
  end QuadDobl_Complex_Vector_Norm_in_Ada;

  procedure Standard_Complex_Vector_Norm_to_C ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   The vector is passed to C for the norm computation.

    v : Standard_Complex_Vectors.Vector(1..dim);
    return_of_call : integer32;

    function normC ( n : integer32;        -- n is the dimension
                     x : C_Double_Array;   -- contains 2*n doubles
                     y : C_Double_Array )  -- on return is y(0) 
                   return integer32;
    pragma import(C, normC, "gpu2norm_d");

  begin
    v := Standard_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    declare
      v_a : constant C_Double_Array := C_to_Ada_Arrays.Convert(v);
      res : C_Double_Array(0..0);
    begin
      res(res'first) := 0.0;
      return_of_call := normC(dim,v_a,res);
      put("The result of C : "); put(double_float(res(0))); new_line;
    end;
  end Standard_Complex_Vector_Norm_to_C;

  procedure DoblDobl_Complex_Vector_Norm_to_C ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   The vector is passed to C for the norm computation.

    v : DoblDobl_Complex_Vectors.Vector(1..dim);
    return_of_call : integer32;

    function normC ( n : integer32;        -- n is the dimension
                     x : C_Double_Array;   -- contains 4*n doubles
                     y : C_Double_Array )  -- on return is y(0) 
                   return integer32;
    pragma import(C, normC, "gpu2norm_dd");

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
      hi := double_float(v_a(0)); -- res(0));
      lo := double_float(v_a(1)); -- res(1));
      put("hi = "); put(hi); new_line;
      put("lo = "); put(lo); new_line;
      result := Double_Double_Numbers.Create(hi,lo);
      put("The result of C : "); put(result); new_line;
     -- put("The sqrt of it  : "); put(SQRT(result)); new_line;
    end;
  end DoblDobl_Complex_Vector_Norm_to_C;

  procedure QuadDobl_Complex_Vector_Norm_to_C ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex vector of range 1..dim
  --   and computes its 2-norm, which should equal sqrt(dim)
  --   for a vector where all numbers are on the complex unit circle.
  --   The vector is passed to C for the norm computation.

    v : QuadDobl_Complex_Vectors.Vector(1..dim);
    return_of_call : integer32;

    function normC ( n : integer32;        -- n is the dimension
                     x : C_Double_Array;   -- contains 8*n doubles
                     y : C_Double_Array )  -- on return is y(0) 
                   return integer32;
    pragma import(C, normC, "gpu2norm_qd");

  begin
   -- Standard_Random_Numbers.Set_Seed(1234);
    v := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    put_line("A random complex vector v : "); put_line(v);
    declare
      v_a : constant C_Double_Array := C_to_Ada_Arrays.Convert(v);
      res : C_Double_Array(0..3);
      hihi,hilo,lohi,lolo : double_float;
     -- re0 : quad_double := QuadDobl_Complex_Numbers.REAL_PART(v(1));
      result : quad_double;
    begin
      res(0) := 0.0;
      res(1) := 0.0;
      res(2) := 0.0;
      res(3) := 0.0;
      -- put("v[0] = "); put(hihi_part(re0)); new_line;
      -- put("v[1] = "); put(lohi_part(re0)); new_line;
      -- put("v[2] = "); put(hilo_part(re0)); new_line;
      -- put("v[3] = "); put(lolo_part(re0)); new_line;
      -- put("v_a[0] = "); put(double_float(v_a(0))); new_line;
      -- put("v_a[1] = "); put(double_float(v_a(1))); new_line;
      -- put("v_a[2] = "); put(double_float(v_a(2))); new_line;
      -- put("v_a[3] = "); put(double_float(v_a(3))); new_line;
      return_of_call := normC(dim,v_a,res);
      hihi := double_float(v_a(0)); -- res(0));
      lohi := double_float(v_a(1)); -- res(1));
      hilo := double_float(v_a(2)); -- res(2));
      lolo := double_float(v_a(3)); -- res(3));
      -- put("hihi = "); put(hihi); new_line;
      -- put("lohi = "); put(lohi); new_line;
      -- put("hilo = "); put(hilo); new_line;
      -- put("lolo = "); put(lolo); new_line;
      result := Quad_Double_Numbers.Create(hihi,lohi,hilo,lolo);
      put("The result of C : "); put(result); new_line;
     -- put("The sqrt of it  : "); put(SQRT(result)); new_line;
    end;
  end QuadDobl_Complex_Vector_Norm_to_C;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension, generates a random vector
  --   and then computes the norm in Ada and C in double precision.

    n : integer32 := 0;

  begin
    new_line;
    put_line("Accelerating in hardware double precision ...");
    new_line;
    put("Give the dimension : "); get(n);
    put("-> the dimension is "); put(n,1); put_line(".");
    new_line;
    Standard_Complex_Vector_Norm_in_Ada(n);
    Standard_Complex_Vector_Norm_to_C(n);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension, generates a random vector
  --   and then computes the norm in Ada and C in double double precision.

    n : integer32 := 0;

  begin
    new_line;
    put_line("Accelerating in double double precision ...");
    new_line;
    put("Give the dimension : "); get(n);
    put("-> the dimension is "); put(n,1); put_line(".");
    new_line;
    DoblDobl_Complex_Vector_Norm_in_Ada(n);
    DoblDobl_Complex_Vector_Norm_to_C(n);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension, generates a random vector
  --   and then computes the norm in Ada and C in quad double precision.

    n : integer32 := 0;

  begin
    new_line;
    put_line("Accelerating in quad double precision ...");
    new_line;
    put("Give the dimension : "); get(n);
    put("-> the dimension is "); put(n,1); put_line(".");
    new_line;
    QuadDobl_Complex_Vector_Norm_in_Ada(n);
    QuadDobl_Complex_Vector_Norm_to_C(n);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches
  --   the proper main program.

    ans : character;

  begin
    new_line;
    put_line("Accelerating the 2-norm computation of a complex vector ...");
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_gpu2norm;
