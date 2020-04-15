with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Hyperplane_Solution_Scaling;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Series_Path_Trackers;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy_Convolutions_io;
with Hyperplane_Convolution_Scaling;

procedure ts_scalecnv is

-- DESCRIPTION :
--   Development of homotopy convolution circuits with scaling of
--   the solutions, after a projective coordinate transformation.

  procedure Standard_Test
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Checks whether the solutions satisfy the last equation
  --   of the homotopy, to test whether the projective transformation
  --   was executed correctly, in double precision.

    use Standard_Complex_Numbers,Standard_Complex_Solutions;
    use Standard_Speelpenning_Convolutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    val : Complex_Number;
    zero : constant Complex_Number := Create(0.0);
    ans : character;

  begin
    new_line;
    put_line("The coefficients of the last circuit : ");
    put(chom.crc(chom.crc'last).cst);
    put(chom.crc(chom.crc'last).cff);
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put("Solution vector "); put(k,1); put_line(" :"); put_line(ls.v);
      val := Eval(chom.crc(chom.crc'last),ls.v,zero);
      put_line("The vector evaluated at the last equation :");
      put(val); new_line;
      Hyperplane_Solution_Scaling.Scale(ls.v);
      put_line("Solution vector after scaling :"); put_line(ls.v);
      Hyperplane_Convolution_Scaling.Adjust_Last_Constant(chom,ls.v);
      val := Eval(chom.crc(chom.crc'last),ls.v,zero);
      put_line("The vector evaluated at the last equation :");
      put(val); new_line;
      tmp := Tail_Of(tmp);
      if not Is_Null(tmp) then
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
    end loop;
  end Standard_Test;

  procedure DoblDobl_Test
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Checks whether the solutions satisfy the last equation
  --   of the homotopy, to test whether the projective transformation
  --   was executed correctly, in double double precision.

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Solutions;
    use DoblDobl_Speelpenning_Convolutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    val : Complex_Number;
    zero : constant Complex_Number := Create(integer(0));
    ans : character;

  begin
    new_line;
    put_line("The coefficients of the last circuit : ");
    put(chom.crc(chom.crc'last).cst);
    put(chom.crc(chom.crc'last).cff);
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put("Solution vector "); put(k,1); put_line(" :"); put_line(ls.v);
      val := Eval(chom.crc(chom.crc'last),ls.v,zero);
      put_line("The vector evaluated at the last equation :");
      put(val); new_line;
      Hyperplane_Solution_Scaling.Scale(ls.v);
      put_line("Solution vector after scaling :"); put_line(ls.v);
      Hyperplane_Convolution_Scaling.Adjust_Last_Constant(chom,ls.v);
      val := Eval(chom.crc(chom.crc'last),ls.v,zero);
      put_line("The vector evaluated at the last equation :");
      put(val); new_line;
      tmp := Tail_Of(tmp);
      if not Is_Null(tmp) then
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Checks whether the solutions satisfy the last equation
  --   of the homotopy, to test whether the projective transformation
  --   was executed correctly, in quad double precision.

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Solutions;
    use QuadDobl_Speelpenning_Convolutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    val : Complex_Number;
    zero : constant Complex_Number := Create(integer(0));
    ans : character;

  begin
    new_line;
    put_line("The coefficients of the last circuit : ");
    put(chom.crc(chom.crc'last).cst);
    put(chom.crc(chom.crc'last).cff);
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put("Solution vector "); put(k,1); put_line(" :"); put_line(ls.v);
      val := Eval(chom.crc(chom.crc'last),ls.v,zero);
      put_line("The vector evaluated at the last equation :");
      put(val); new_line;
      tmp := Tail_Of(tmp);
      Hyperplane_Solution_Scaling.Scale(ls.v);
      put_line("Solution vector after scaling :"); put_line(ls.v);
      Hyperplane_Convolution_Scaling.Adjust_Last_Constant(chom,ls.v);
      val := Eval(chom.crc(chom.crc'last),ls.v,zero);
      put_line("The vector evaluated at the last equation :");
      put(val); new_line;
      tmp := Tail_Of(tmp);
      if not Is_Null(tmp) then
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
    end loop;
  end QuadDobl_Test;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Prompts the user for the settings of the homotopy,
  --   in standard double precision.

    art : constant boolean := Series_Path_Trackers.Prompt_for_Artificial;
    deg,idxpar : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    if not art
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get
      (deg,art,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
    Standard_Test(cnvhom,sols);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for the settings of the homotopy,
  --   in double double precision.

    art : constant boolean := Series_Path_Trackers.Prompt_for_Artificial;
    deg,idxpar : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    cnvhom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    if not art
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    DoblDobl_Homotopy_Convolutions_io.get
      (deg,art,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
    DoblDobl_Test(cnvhom,sols);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Prompts the user for the settings of the homotopy,
  --   in quad double precision.

    art : constant boolean := Series_Path_Trackers.Prompt_for_Artificial;
    deg,idxpar : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    cnvhom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    if not art
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    QuadDobl_Homotopy_Convolutions_io.get
      (deg,art,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
    QuadDobl_Test(cnvhom,sols);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then
  --   launches the corresponding test.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_scalecnv;
