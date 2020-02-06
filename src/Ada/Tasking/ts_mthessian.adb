with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Hessian_Convolution_Circuits;       use Hessian_Convolution_Circuits;
with Multitasked_Hessian_Convolutions;   use Multitasked_Hessian_Convolutions;

procedure ts_mthessian is

-- DESCRIPTION :
--   Tests the development of the Hessian criterion,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : Standard_Complex_Vectors.Link_to_Vector;
    val : double_float;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := Standard_Complex_Numbers.REAL_PART(lnk(i));
        put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in DoblDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    val : double_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := DoblDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Write_Singular_Values 
              ( values : in QuadDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    val : quad_double;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := QuadDobl_Complex_Numbers.REAL_PART(lnk(i));
        put(" "); put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    vx : Standard_Complex_Vectors.Vector(1..dim);
    A,U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : Standard_Complex_VecVecs.VecVec(1..dim);
    values : Standard_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    nbt : integer32 := 0;

  begin
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    for i in 1..dim loop
      declare
        vls : Standard_Complex_Vectors.Vector(1..dim);
      begin
        Singular_Values(s.crc(i),vx,A,U,V,e,vls);
        svl(i) := new Standard_Complex_Vectors.Vector'(vls);
      end;
    end loop;
    put_line("All singular values :");
    Write_Singular_Values(svl);
    new_line;
    put("Give the number of tasks : "); get(nbt);
    Multitasked_Singular_Values(nbt,s,vx,values);
    put_line("All singular values computed by multitasking :");
    Write_Singular_Values(values);
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    vx : DoblDobl_Complex_Vectors.Vector(1..dim);
    A,U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : DoblDobl_Complex_Vectors.Vector(1..dim);
    svl : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    values : DoblDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    nbt : integer32 := 0;

  begin
    DoblDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    for i in 1..dim loop
      declare
        vls : DoblDobl_Complex_Vectors.Vector(1..dim);
      begin
        Singular_Values(s.crc(i),vx,A,U,V,e,vls);
        svl(i) := new DoblDobl_Complex_Vectors.Vector'(vls);
      end;
    end loop;
    put_line("All singular values :");
    Write_Singular_Values(svl);
    new_line;
    put("Give the number of tasks : "); get(nbt);
    Multitasked_Singular_Values(nbt,s,vx,values);
    put_line("All singular values computed by multitasking :");
    Write_Singular_Values(values);
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    vx : QuadDobl_Complex_Vectors.Vector(1..dim);
    A,U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    e : QuadDobl_Complex_Vectors.Vector(1..dim);
    svl : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    values : QuadDobl_Complex_VecVecs.VecVec(1..dim) := Allocate(dim,dim);
    nbt : integer32 := 0;

  begin
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    for i in 1..dim loop
      declare
        vls : QuadDobl_Complex_Vectors.Vector(1..dim);
      begin
        Singular_Values(s.crc(i),vx,A,U,V,e,vls);
        svl(i) := new QuadDobl_Complex_Vectors.Vector'(vls);
      end;
    end loop;
    put_line("All singular values :");
    Write_Singular_Values(svl);
    new_line;
    put("Give the number of tasks : "); get(nbt);
    Multitasked_Singular_Values(nbt,s,vx,values);
    put_line("All singular values computed by multitasking :");
    Write_Singular_Values(values);
  end QuadDobl_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Launches the test after prompting for the parameters
  --   to generate a random problem.

    dim,deg,nbr,pwr : integer32 := 0;
    precision : character;

  begin
    new_line;
    put_line("Testing the Hessian criterion ...");
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the power series : "); get(deg);
    put("Give the number of terms in each circuit : "); get(nbr);
    put("Give the largest power of the variables : "); get(pwr);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    case precision is
      when '0' => Standard_Random_Test(dim,deg,nbr,pwr);
      when '1' => DoblDobl_Random_Test(dim,deg,nbr,pwr);
      when '2' => QuadDobl_Random_Test(dim,deg,nbr,pwr);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mthessian;
