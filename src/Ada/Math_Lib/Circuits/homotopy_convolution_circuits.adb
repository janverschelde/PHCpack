with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;

package body Homotopy_Convolution_Circuits is

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit ) is

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
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

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
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    cf : TripDobl_Complex_Vectors.Link_to_Vector;
    one : constant triple_double := create(1.0);

    use TripDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := TripDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := TripDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

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
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    cf : PentDobl_Complex_Vectors.Link_to_Vector;
    one : constant penta_double := create(1.0);

    use PentDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := PentDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := PentDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    cf : OctoDobl_Complex_Vectors.Link_to_Vector;
    one : constant octo_double := create(1.0);

    use OctoDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := OctoDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := OctoDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    cf : DecaDobl_Complex_Vectors.Link_to_Vector;
    one : constant deca_double := create(1.0);

    use DecaDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := DecaDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := DecaDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    cf : HexaDobl_Complex_Vectors.Link_to_Vector;
    one : constant hexa_double := create(1.0);

    use HexaDobl_Complex_Vectors;

  begin
    for k in c.cff'range loop
      cf := c.cff(k);
      cf(1) := HexaDobl_Complex_Numbers.Create(one);
    end loop;
    if c.cst /= null then
      c.cst(1) := HexaDobl_Complex_Numbers.Create(one);
    end if;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Continuation_Parameter
              ( c : in HexaDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in c'range loop
      Add_Continuation_Parameter(c(k));
    end loop;
  end Add_Continuation_Parameter;

  procedure Add_Parameter_to_Constant
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
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
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
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
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 ) is

    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new TripDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
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
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 ) is

    use PentDobl_Complex_Numbers;
    use PentDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new PentDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 ) is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new OctoDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 ) is

    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new DecaDobl_Complex_Vectors.Vector'
                     (0..deg => Create(integer32(0)));
      c.cst(1) := Create(integer32(1));
    end if;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 ) is

    use HexaDobl_Complex_Numbers;
    use HexaDobl_Complex_Vectors;

  begin
    if c.cst /= null then
      c.cst(1) := Create(integer32(1));
    else
      c.cst := new HexaDobl_Complex_Vectors.Vector'
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
              ( s : in TripDobl_Speelpenning_Convolutions.Link_to_System ) is
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

  procedure Add_Parameter_to_Constant
              ( s : in PentDobl_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in OctoDobl_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in DecaDobl_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Add_Parameter_to_Constant
              ( s : in HexaDobl_Speelpenning_Convolutions.Link_to_System ) is
  begin
    for k in s.crc'range loop
      Add_Parameter_to_Constant(s.crc(k),s.deg);
    end loop;
  end Add_Parameter_to_Constant;

  procedure Set_Solution_Constant
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                z : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new Standard_Complex_Vectors.Vector'(0..deg => Create(0.0));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant double_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new DoblDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in TripDobl_Complex_Vectors.Vector ) is

    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Vectors;
    use TripDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant triple_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new TripDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant quad_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new QuadDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in PentDobl_Complex_Vectors.Vector ) is

    use PentDobl_Complex_Numbers;
    use PentDobl_Complex_Vectors;
    use PentDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant penta_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new PentDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in OctoDobl_Complex_Vectors.Vector ) is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Vectors;
    use OctoDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant octo_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new OctoDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DecaDobl_Complex_Vectors.Vector ) is

    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Vectors;
    use DecaDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant deca_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new DecaDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in HexaDobl_Complex_Vectors.Vector ) is

    use HexaDobl_Complex_Numbers;
    use HexaDobl_Complex_Vectors;
    use HexaDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant hexa_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
      else
        deg := c.cff(1)'last;
        c.cst := new HexaDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
      end if;
    end if;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                z : in Standard_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                z : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits;
                z : in TripDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                z : in QuadDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits;
                z : in PentDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                z : in OctoDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                z : in DecaDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Set_Solution_Constant
              ( c : in HexaDobl_Speelpenning_Convolutions.Circuits;
                z : in HexaDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Set_Solution_Constant(c(k),z);
    end loop;
  end Set_Solution_Constant;

  procedure Newton_Homotopy
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                z : in Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(0.0);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new Standard_Complex_Vectors.Vector'(0..deg => Create(0.0));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant double_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new DoblDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in TripDobl_Complex_Vectors.Vector ) is

    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Vectors;
    use TripDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant triple_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new TripDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy 
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant quad_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new QuadDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy 
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in PentDobl_Complex_Vectors.Vector ) is

    use PentDobl_Complex_Numbers;
    use PentDobl_Complex_Vectors;
    use PentDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant penta_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new PentDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy 
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in OctoDobl_Complex_Vectors.Vector ) is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Vectors;
    use OctoDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant octo_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new OctoDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy 
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DecaDobl_Complex_Vectors.Vector ) is

    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Vectors;
    use DecaDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant deca_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new DecaDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy 
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in HexaDobl_Complex_Vectors.Vector ) is

    use HexaDobl_Complex_Numbers;
    use HexaDobl_Complex_Vectors;
    use HexaDobl_Speelpenning_Convolutions;

    y : Complex_Number;
    deg : integer32;
    zero : constant hexa_double := create(0.0);

  begin
    if c /= null then
      y := Eval(c,z);
      if c.cst /= null then
        c.cst(0) := c.cst(0) - y;
        c.cst(1) := y;
        for k in 2..c.cst'last loop
          c.cst(k) := Create(zero);
        end loop;
      else
        deg := c.cff(1)'last;
        c.cst := new HexaDobl_Complex_Vectors.Vector'(0..deg => Create(zero));
        c.cst(0) := -y;
        c.cst(1) := y;   -- c turned into c - y + t*y
      end if;
    end if;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                z : in Standard_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                z : in DoblDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits;
                z : in TripDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                z : in QuadDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits;
                z : in PentDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                z : in OctoDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                z : in DecaDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

  procedure Newton_Homotopy
              ( c : in HexaDobl_Speelpenning_Convolutions.Circuits;
                z : in HexaDobl_Complex_Vectors.Vector ) is
  begin
    for k in c'range loop
      Newton_Homotopy(c(k),z);
    end loop;
  end Newton_Homotopy;

end Homotopy_Convolution_Circuits;
