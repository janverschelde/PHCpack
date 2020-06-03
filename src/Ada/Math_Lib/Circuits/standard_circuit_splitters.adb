with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Vector_Splitters;

package body Standard_Circuit_Splitters is

  function Split ( c : Standard_Complex_Circuits.Circuit )
                 return Standard_Coefficient_Circuits.Circuit is

    res : Standard_Coefficient_Circuits.Circuit(c.nbr);
    fwd : constant Standard_Floating_Vectors.Vector(1..c.dim-1)
        := (1..c.dim-1 => 0.0);
    bck : constant Standard_Floating_Vectors.Vector(1..c.dim-2)
        := (1..c.dim-2 => 0.0);
    crs : constant Standard_Floating_Vectors.Vector(1..c.dim-2)
        := (1..c.dim-2 => 0.0);
 
  begin
    res.dim := c.dim;
    res.xps := c.xps;
    res.idx := c.idx;
    res.fac := c.fac;
    Standard_Vector_Splitters.Split_Complex(c.cff,res.rcf,res.icf);
    res.rcst := Standard_Complex_Numbers.REAL_PART(c.cst);
    res.icst := Standard_Complex_Numbers.IMAG_PART(c.cst);
    res.rfwd := new Standard_Floating_Vectors.Vector'(fwd);
    res.ifwd := new Standard_Floating_Vectors.Vector'(fwd);
    res.rbck := new Standard_Floating_Vectors.Vector'(bck);
    res.ibck := new Standard_Floating_Vectors.Vector'(bck);
    res.rcrs := new Standard_Floating_Vectors.Vector'(crs);
    res.icrs := new Standard_Floating_Vectors.Vector'(crs);
    return res;
  end Split;

  function Split ( c : Standard_Complex_Circuits.Link_to_Circuit )
                 return Standard_Coefficient_Circuits.Link_to_Circuit is

    res : Standard_Coefficient_Circuits.Link_to_Circuit;

    use Standard_Complex_Circuits;

  begin
    if c /= null then
      declare
        cs : Standard_Coefficient_Circuits.Circuit(c.nbr);
      begin
        cs := Split(c.all);
        res := new Standard_Coefficient_Circuits.Circuit'(cs);
      end;
    end if;
    return res;
  end Split;

  function Split ( c : Standard_Complex_Circuits.Circuits )
                 return Standard_Coefficient_Circuits.Circuits is

    res : Standard_Coefficient_Circuits.Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := Split(c(k));
    end loop;
    return res;
  end Split;

  function Split ( s : Standard_Complex_Circuits.System )
                 return Standard_Coefficient_Circuits.System is

    res : Standard_Coefficient_Circuits.System(s.neq,s.dim);

  begin
    res.crc := Split(s.crc);
    res.mxe := s.mxe;
    Standard_Vector_Splitters.Split_Complex(s.pwt,res.rpwt,res.ipwt);
    Standard_Vector_Splitters.Split_Complex(s.yd,res.ryd,res.iyd);
    Standard_Coefficient_Circuits.Allocate_Hessian_Space
      (s.dim,res.hrp,res.hip);
    return res;
  end Split;

  function Split ( s : Standard_Complex_Circuits.Link_to_System )
                 return Standard_Coefficient_Circuits.Link_to_System is

    res : Standard_Coefficient_Circuits.Link_to_System;

    use Standard_Complex_Circuits;

  begin
    if s /= null then
      declare
        cs : Standard_Coefficient_Circuits.System(s.neq,s.dim);
      begin
        cs := Split(s.all);
        res := new Standard_Coefficient_Circuits.System'(cs);
      end;
    end if;
    return res;
  end Split;

end Standard_Circuit_Splitters;
