with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;

package body Test_Predictor_Convolutions is

  procedure Standard_Check_Solutions
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    put_line("Checking the start solutions ...");
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        y : constant Standard_Complex_Vectors.Vector
          := Standard_Speelpenning_Convolutions.Eval(chom.crc,ls.v,zero);
      begin
        put("Value at solution "); put(k,1); put_line(" :"); put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Check_Solutions;

  procedure DoblDobl_Check_Solutions
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));

  begin
    put_line("Checking the start solutions ...");
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        y : constant DoblDobl_Complex_Vectors.Vector
          := DoblDobl_Speelpenning_Convolutions.Eval(chom.crc,ls.v,zero);
      begin
        put("Value at solution "); put(k,1); put_line(" :"); put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Check_Solutions;

  procedure QuadDobl_Check_Solutions
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));

  begin
    put_line("Checking the start solutions ...");
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        y : constant QuadDobl_Complex_Vectors.Vector
          := QuadDobl_Speelpenning_Convolutions.Eval(chom.crc,ls.v,zero);
      begin
        put("Value at solution "); put(k,1); put_line(" :"); put_line(y);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end QuadDobl_Check_Solutions;

end Test_Predictor_Convolutions;
