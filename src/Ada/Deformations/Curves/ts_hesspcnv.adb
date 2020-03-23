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
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
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
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy_Convolutions_io;
with Jacobian_Convolution_Circuits;
with Hessian_Convolution_Circuits;
with Test_Predictor_Convolutions;

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
              ( cnvh : in Standard_Speelpenning_Convolutions.Link_to_System;
                solv : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given a homotopy with a start solution, evaluates the Jacobian
  --   matrix with the symbolic procedures and compares with the outcome
  --   of the convolution circuits, in double precision.

  -- ON ENTRY :
  --   cnvh     convolution circuit representation of the homotopy,
  --            the continuation parameter t is the series parameter;
  --   solv     a start solution for t = 0.

    polh : constant Standard_Complex_Poly_Systems.Poly_Sys
         := Standard_Homotopy.Homotopy_System;
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
              ( cnvh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                solv : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given a homotopy with a start solution, evaluates the Jacobian
  --   matrix with the symbolic procedures and compares with the outcome
  --   of the convolution circuits, in double double precision.

  -- ON ENTRY :
  --   cnvh     convolution circuit representation of the homotopy,
  --            the continuation parameter t is the series parameter;
  --   solv     a start solution for t = 0.

    polh : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
         := DoblDobl_Homotopy.Homotopy_System;
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
              ( cnvh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                solv : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Given a homotopy with a start solution, evaluates the Jacobian
  --   matrix with the symbolic procedures and compares with the outcome
  --   of the convolution circuits, in double double precision.

  -- ON ENTRY :
  --   cnvh     convolution circuit representation of the homotopy,
  --            the continuation parameter t is the series parameter;
  --   solv     a start solution for t = 0.

    polh : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
         := QuadDobl_Homotopy.Homotopy_System;
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

  procedure Standard_Test_Prediction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;

  begin
    Standard_Homotopy_Convolutions_io.get(0,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.Standard_Check_Solutions(cnvhom,sols);
    end if;
    Standard_Test(cnvhom,Standard_Complex_Solutions.Head_Of(sols).v);
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;

  begin
    DoblDobl_Homotopy_Convolutions_io.get(0,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.DoblDobl_Check_Solutions(cnvhom,sols);
    end if;
    DoblDobl_Test(cnvhom,DoblDobl_Complex_Solutions.Head_Of(sols).v);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar : integer32;
    ans : character;

  begin
    QuadDobl_Homotopy_Convolutions_io.get(0,cnvhom,sols,idxpar);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.QuadDobl_Check_Solutions(cnvhom,sols);
    end if;
    QuadDobl_Test(cnvhom,QuadDobl_Complex_Solutions.Head_Of(sols).v);
  end QuadDobl_Test_Prediction;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Test_Prediction;
      when '1' => DoblDobl_Test_Prediction;
      when '2' => QuadDobl_Test_Prediction;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hesspcnv;
