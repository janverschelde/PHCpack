with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;

package body Homotopy_Convolution_Circuits is

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit ) is

    cf : Standard_Complex_Vectors.Link_to_Vector;

    use Standard_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := Standard_Complex_Numbers.Create(1.0);
    end loop;
    if c.cst /= null then
      c.cst(1) := Standard_Complex_Numbers.Create(1.0);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit ) is

    cf : DoblDobl_Complex_Vectors.Link_to_Vector;
    one : constant double_double := create(1.0);

    use DoblDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := DoblDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := DoblDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit ) is

    cf : QuadDobl_Complex_Vectors.Link_to_Vector;
    one : constant quad_double := create(1.0);

    use QuadDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := QuadDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := QuadDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.
                       Convolution_Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Convolution_Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Convolution_Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Parameter_to_Constant
              ( c : in Standard_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(1.0);
    else
      c.cst := new Standard_Complex_Vectors.Vector'(0..deg => Create(0.0));
      c.cst(1) := Create(1.0);
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new DoblDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new QuadDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

end Homotopy_Convolution_Circuits;
