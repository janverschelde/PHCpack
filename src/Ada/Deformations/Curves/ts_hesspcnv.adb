with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
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
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Hessians;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Hessians;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Hessians;
with Singular_Values_of_Hessians;        use Singular_Values_of_Hessians;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Solution_Drops;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Series_and_Homotopies;
with Test_Series_Predictors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Jacobian_Convolution_Circuits;
with Hessian_Convolution_Circuits;

procedure ts_hesspcnv is

-- DESCRIPTION :
--   Development of the Hessian criterion on convolution circuits.

  function AbsSum ( m : Standard_Complex_Matrices.Matrix )
                  return double_float is

  -- DESCRIPTION :
  --   Returns the absolute value of the sum of all elements in m.

    res : double_float := 0.0;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        res := res + Standard_Complex_Numbers.AbsVal(m(i,j));
      end loop;
    end loop;
    return res;
  end AbsSum;

  function AbsSum ( m : DoblDobl_Complex_Matrices.Matrix )
                  return double_double is

  -- DESCRIPTION :
  --   Returns the absolute value of the sum of all elements in m.

    res : double_double := create(0.0);

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        res := res + DoblDobl_Complex_Numbers.AbsVal(m(i,j));
      end loop;
    end loop;
    return res;
  end AbsSum;

  function AbsSum ( m : QuadDobl_Complex_Matrices.Matrix )
                  return quad_double is

  -- DESCRIPTION :
  --   Returns the absolute value of the sum of all elements in m.

    res : quad_double := create(0.0);

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        res := res + QuadDobl_Complex_Numbers.AbsVal(m(i,j));
      end loop;
    end loop;
    return res;
  end AbsSum;

  procedure Standard_Test
              ( polh : in Standard_Complex_Poly_Systems.Poly_Sys;
                cnvh : in Standard_Speelpenning_Convolutions.Link_to_System;
                solv : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given a homotopy with a start solution, evaluates the Jacobian
  --   matrix with the symbolic procedures and compares with the outcome
  --   of the convolution circuits, in double precision.

  -- ON ENTRY :
  --   polh     polynomials in an artificial-parameter homotopy,
  --            the last variable in polh is the continuation parameter t;
  --   cnvh     convolution circuit representation of the homotopy,
  --            the continuation parameter t is the series parameter.

    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    polht0 : constant Standard_Complex_Poly_Systems.Poly_Sys(polh'range)
           := Standard_Complex_Poly_SysFun.Eval(polh,zero,solv'last+1);
    jac : Standard_Complex_Jaco_Matrices.Jaco_Mat(polh'range,solv'range)
        := Standard_Complex_Jaco_Matrices.Create(polht0);
    ejm1 : Standard_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    ejm2 : Standard_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    dff : Standard_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    vm : Standard_Complex_VecMats.VecMat(cnvh.crc'range);
    val : double_float;
    hjm : Standard_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hss : Standard_Complex_Hessians.Link_to_Array_of_Hessians;
    solt0 : Standard_Complex_Vectors.Vector(solv'first..solv'last+1);

    use Standard_Complex_Matrices;

  begin
    ejm1 := Standard_Complex_Jaco_Matrices.Eval(jac,solv);
    ejm2 := Jacobian_Convolution_Circuits.Jacobian(cnvh.crc,solv);
    dff := ejm1 - ejm2;
    val := AbsSum(dff);
    put("Absolute sum of differences of Jacobians :"); put(val,3); new_line;
    Standard_Complex_Jaco_Matrices.Clear(jac);
    Standard_Jacobian_Hessians_of_Homotopy(hjm,hss);
    vm := Hessian_Convolution_Circuits.Hessians(cnvh.crc,solv);
    solt0(solv'range) := solv;
    solt0(solt0'last) := zero;
    for k in vm'range loop
      declare
        H : constant Standard_Complex_Matrices.Matrix(polh'range,polh'range)
          := Standard_Complex_Hessians.eval(hss(k),solt0);
        Hdiff : Standard_Complex_Matrices.Matrix(H'range(1),H'range(2));
      begin
        Hdiff := vm(k).all - H; val := AbsSum(Hdiff);
        put("Hessian "); put(k,1); put(" diffsum :"); put(val,3); new_line;
      end;
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test
              ( polh : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                cnvh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                solv : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given a homotopy with a start solution, evaluates the Jacobian
  --   matrix with the symbolic procedures and compares with the outcome
  --   of the convolution circuits, in double double precision.

  -- ON ENTRY :
  --   polh     polynomials in an artificial-parameter homotopy,
  --            the last variable in polh is the continuation parameter t;
  --   cnvh     convolution circuit representation of the homotopy,
  --            the continuation parameter t is the series parameter;
  --   solv     a start solution for t = 0.

    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    polht0 : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(polh'range)
           := DoblDobl_Complex_Poly_SysFun.Eval(polh,zero,solv'last+1);
    jac : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(polh'range,solv'range)
        := DoblDobl_Complex_Jaco_Matrices.Create(polht0);
    ejm1 : DoblDobl_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    ejm2 : DoblDobl_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    dff : DoblDobl_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    vm : DoblDobl_Complex_VecMats.VecMat(cnvh.crc'range);
    val : double_double;
    hjm : DoblDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hss : DoblDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    solt0 : DoblDobl_Complex_Vectors.Vector(solv'first..solv'last+1);

    use DoblDobl_Complex_Matrices;

  begin
    ejm1 := DoblDobl_Complex_Jaco_Matrices.Eval(jac,solv);
    ejm2 := Jacobian_Convolution_Circuits.Jacobian(cnvh.crc,solv);
    dff := ejm1 - ejm2;
    val := AbsSum(dff);
    put("Absolute sum of differences of Jacobians : "); put(val,3); new_line;
    DoblDobl_Complex_Jaco_Matrices.Clear(jac);
    vm := Hessian_Convolution_Circuits.Hessians(cnvh.crc,solv);
    DoblDobl_Jacobian_Hessians_of_Homotopy(hjm,hss);
    solt0(solv'range) := solv;
    solt0(solt0'last) := zero;
    for k in vm'range loop
      declare
        H : constant DoblDobl_Complex_Matrices.Matrix(polh'range,polh'range)
          := DoblDobl_Complex_Hessians.eval(hss(k),solt0);
        Hdiff : DoblDobl_Complex_Matrices.Matrix(H'range(1),H'range(2));
      begin
        Hdiff := vm(k).all - H; val := AbsSum(Hdiff);
        put("Hessian "); put(k,1); put(" diffsum : "); put(val,3); new_line;
      end;
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( polh : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                cnvh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                solv : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given a homotopy with a start solution, evaluates the Jacobian
  --   matrix with the symbolic procedures and compares with the outcome
  --   of the convolution circuits, in double double precision.

  -- ON ENTRY :
  --   polh     polynomials in an artificial-parameter homotopy,
  --            the last variable in polh is the continuation parameter t;
  --   cnvh     convolution circuit representation of the homotopy,
  --            the continuation parameter t is the series parameter;
  --   solv     a start solution for t = 0.

    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    polht0 : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(polh'range)
           := QuadDobl_Complex_Poly_SysFun.Eval(polh,zero,solv'last+1);
    jac : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(polh'range,solv'range)
        := QuadDobl_Complex_Jaco_Matrices.Create(polht0);
    ejm1 : QuadDobl_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    ejm2 : QuadDobl_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    dff : QuadDobl_Complex_Matrices.Matrix(jac'range(1),jac'range(2));
    val : quad_double;
    vm : QuadDobl_Complex_VecMats.VecMat(cnvh.crc'range);
    hjm : QuadDobl_Complex_Jaco_Matrices.Link_to_Jaco_Mat;
    hss : QuadDobl_Complex_Hessians.Link_to_Array_of_Hessians;
    solt0 : QuadDobl_Complex_Vectors.Vector(solv'first..solv'last+1);

    use QuadDobl_Complex_Matrices;

  begin
    ejm1 := QuadDobl_Complex_Jaco_Matrices.Eval(jac,solv);
    ejm2 := Jacobian_Convolution_Circuits.Jacobian(cnvh.crc,solv);
    dff := ejm1 - ejm2;
    val := AbsSum(dff);
    put("Absolute sum of differences of Jacobians : "); put(val,3); new_line;
    QuadDobl_Complex_Jaco_Matrices.Clear(jac);
    vm := Hessian_Convolution_Circuits.Hessians(cnvh.crc,solv);
    QuadDobl_Jacobian_Hessians_of_Homotopy(hjm,hss);
    solt0(solv'range) := solv;
    solt0(solt0'last) := zero;
    for k in vm'range loop
      declare
        H : constant QuadDobl_Complex_Matrices.Matrix(polh'range,polh'range)
          := QuadDobl_Complex_Hessians.eval(hss(k),solt0);
        Hdiff : QuadDobl_Complex_Matrices.Matrix(H'range(1),H'range(2));
      begin
        Hdiff := vm(k).all - H; val := AbsSum(Hdiff);
        put("Hessian "); put(k,1); put(" diffsum : "); put(val,3); new_line;
      end;
    end loop;
  end QuadDobl_Test;

  procedure Standard_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The Standard_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    use Standard_Complex_Solutions;

    hom : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
        := Standard_Homotopy.Homotopy_System;
    serhom : constant Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idxpar);
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    cnvhom := Make_Convolution_System(serhom,0);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put_line("Checking the start solutions ...");
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        y : constant Standard_Complex_Vectors.Vector
          := Standard_Speelpenning_Convolutions.Eval(cnvhom.crc,ls.v,zero);
      begin
        put("Value at solution "); put(k,1); put_line(" :");
        put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Test(hom,cnvhom,Head_Of(sols).v);
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The DoblDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    use DoblDobl_Complex_Solutions;

    hom : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
        := DoblDobl_Homotopy.Homotopy_System;
    serhom : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idxpar);
    cnvhom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));

  begin
    cnvhom := Make_Convolution_System(serhom,0);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put_line("Checking the start solutions ...");
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        y : constant DoblDobl_Complex_Vectors.Vector
          := DoblDobl_Speelpenning_Convolutions.Eval(cnvhom.crc,ls.v,zero);
      begin
        put("Value at solution "); put(k,1); put_line(" :");
        put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Test(hom,cnvhom,Head_Of(sols).v);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( nq,idxpar : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The QuadDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    use QuadDobl_Complex_Solutions;

    hom : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
        := QuadDobl_Homotopy.Homotopy_System;
    serhom : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idxpar);
    cnvhom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));

  begin
    cnvhom := Make_Convolution_System(serhom,0);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put_line("Checking the start solutions ...");
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        y : constant QuadDobl_Complex_Vectors.Vector
          := QuadDobl_Speelpenning_Convolutions.Eval(cnvhom.crc,ls.v,zero);
      begin
        put("Value at solution "); put(k,1); put_line(" :");
        put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    QuadDobl_Test(hom,cnvhom,Head_Of(sols).v);
  end QuadDobl_Test_Prediction;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq,idxpar : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      Standard_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant Standard_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        Standard_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq,idxpar : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.DoblDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      DoblDobl_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant DoblDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        DoblDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq,idxpar : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.QuadDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      QuadDobl_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant QuadDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        QuadDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    case precision is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hesspcnv;
