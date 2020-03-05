with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors_cv;        use DoblDobl_Complex_Vectors_cv;
with QuadDobl_Complex_Vectors_cv;        use QuadDobl_Complex_Vectors_cv;
with Varbprec_VecVec_Conversions;
with Standard_Complex_Series;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Series;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Series;
with Exponent_Indices;

package body System_Convolution_Circuits is

  function Is_Zero ( d : Standard_Complex_Polynomials.Degrees )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all entries of d are zero,
  --   return false otherwise.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( d : Standard_CSeries_Polynomials.Degrees )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all entries of d are zero,
  --   return false otherwise.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( d : DoblDobl_Complex_Polynomials.Degrees )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all entries of d are zero,
  --   return false otherwise.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( d : DoblDobl_CSeries_Polynomials.Degrees )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all entries of d are zero,
  --   return false otherwise.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( d : QuadDobl_Complex_Polynomials.Degrees )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all entries of d are zero,
  --   return false otherwise.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Zero ( d : QuadDobl_CSeries_Polynomials.Degrees )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all entries of d are zero,
  --   return false otherwise.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Nonzero_Constant
             ( c : Standard_Complex_Numbers.Complex_Number )
             return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if c is nonzero, 0 otherwise.

    use Standard_Complex_Numbers;

    cmplxzero : constant Complex_Number := Create(0.0);

  begin
    if c = cmplxzero
     then return 0;
     else return 1;
    end if;
  end Nonzero_Constant;

  function Nonzero_Constant
             ( c : DoblDobl_Complex_Numbers.Complex_Number )
             return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if c is nonzero, 0 otherwise.

    use DoblDobl_Complex_Numbers;

    cmplxzero : constant Complex_Number := Create(integer32(0));

  begin
    if c = cmplxzero
     then return 0;
     else return 1;
    end if;
  end Nonzero_Constant;

  function Nonzero_Constant
             ( c : QuadDobl_Complex_Numbers.Complex_Number )
             return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if c is nonzero, 0 otherwise.

    use QuadDobl_Complex_Numbers;

    cmplxzero : constant Complex_Number := Create(integer32(0));

  begin
    if c = cmplxzero
     then return 0;
     else return 1;
    end if;
  end Nonzero_Constant;

  function Nonzero_Constant
             ( s : Standard_Complex_Series.Link_to_Series )
             return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if s /= null, 0 otherwise.

    use Standard_Complex_Series;

  begin
    if s = null
     then return 0;
     else return 1;
    end if;
  end Nonzero_Constant;

  function Nonzero_Constant
             ( s : DoblDobl_Complex_Series.Link_to_Series )
             return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if s /= null, 0 otherwise.

    use DoblDobl_Complex_Series;

  begin
    if s = null
     then return 0;
     else return 1;
    end if;
  end Nonzero_Constant;

  function Nonzero_Constant
             ( s : QuadDobl_Complex_Series.Link_to_Series )
             return integer32 is

  -- DESCRIPTION :
  --   Returns 1 if s /= null, 0 otherwise.

    use QuadDobl_Complex_Series;

  begin
    if s = null
     then return 0;
     else return 1;
    end if;
  end Nonzero_Constant;

  function Make_Convolution_Circuit
             ( p : Standard_Complex_Polynomials.Poly;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Circuit is

    use Standard_Speelpenning_Convolutions;

    deg : constant integer32 := integer32(d);
    dim : constant integer32
        := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p));
    zero : Standard_Complex_Polynomials.Degrees
         := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    cmplxzero : constant Standard_Complex_Numbers.Complex_Number
              := Standard_Complex_Numbers.Create(0.0);
    cst : constant Standard_Complex_Numbers.Complex_Number
        := Standard_Complex_Polynomials.Coeff(p,zero);
    nbr : constant integer32
        := integer32(Standard_Complex_Polynomials.Number_of_Terms(p))
           - Nonzero_Constant(cst);
    res : Circuit(nbr,dim,dim-1,dim-2);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           continue : out boolean ) is

      use Standard_Complex_Polynomials;

    begin
      if not Is_Zero(t.dg) then
        idx := idx + 1;
        declare
          cf : Standard_Complex_Vectors.Vector(0..deg) := (0..deg => cmplxzero);
          xp : Standard_Integer_Vectors.Vector(1..dim);
        begin
          cf(0) := t.cf;
          res.cff(idx) := new Standard_Complex_Vectors.Vector'(cf);
          for i in t.dg'range loop
            xp(i) := integer32(t.dg(i));
          end loop;
          res.xps(idx) := new Standard_Integer_Vectors.Vector'(xp);
        end;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cst := new Standard_Complex_Vectors.Vector'(0..deg => cmplxzero);
    res.cst(0) := cst;
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(zero));
    return res;
  end Make_Convolution_Circuit;

  function Make_Convolution_Circuit
             ( p : DoblDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    deg : constant integer32 := integer32(d);
    dim : constant integer32
        := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p));
    zero : DoblDobl_Complex_Polynomials.Degrees
         := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    cmplxzero : constant DoblDobl_Complex_Numbers.Complex_Number
              := DoblDobl_Complex_Numbers.Create(integer32(0));
    cst : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Complex_Polynomials.Coeff(p,zero);
    nbr : constant integer32
        := integer32(DoblDobl_Complex_Polynomials.Number_of_Terms(p))
           - Nonzero_Constant(cst);
    res : Circuit(nbr,dim,dim-1,dim-2);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      use DoblDobl_Complex_Polynomials;

    begin
      if not Is_Zero(t.dg) then
        idx := idx + 1;
        declare
          cf : DoblDobl_Complex_Vectors.Vector(0..deg) := (0..deg => cmplxzero);
          xp : Standard_Integer_Vectors.Vector(1..dim);
        begin
          cf(0) := t.cf;
          res.cff(idx) := new DoblDobl_Complex_Vectors.Vector'(cf);
          for i in t.dg'range loop
            xp(i) := integer32(t.dg(i));
          end loop;
          res.xps(idx) := new Standard_Integer_Vectors.Vector'(xp);
        end;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cst := new DoblDobl_Complex_Vectors.Vector'(0..deg => cmplxzero);
    res.cst(0) := cst;
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(zero));
    return res;
  end Make_Convolution_Circuit;

  function Make_Convolution_Circuit
             ( p : QuadDobl_Complex_Polynomials.Poly;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    deg : constant integer32 := integer32(d);
    dim : constant integer32
        := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p));
    zero : QuadDobl_Complex_Polynomials.Degrees
         := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    cmplxzero : constant QuadDobl_Complex_Numbers.Complex_Number
              := QuadDobl_Complex_Numbers.Create(integer32(0));
    cst : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Complex_Polynomials.Coeff(p,zero);
    nbr : constant integer32
        := integer32(QuadDobl_Complex_Polynomials.Number_of_Terms(p))
           - Nonzero_Constant(cst);
    res : Circuit(nbr,dim,dim-1,dim-2);
    idx : integer32 := 0;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      use QuadDobl_Complex_Polynomials;

    begin
      if not Is_Zero(t.dg) then
        idx := idx + 1;
        declare
          cf : QuadDobl_Complex_Vectors.Vector(0..deg) := (0..deg => cmplxzero);
          xp : Standard_Integer_Vectors.Vector(1..dim);
        begin
          cf(0) := t.cf;
          res.cff(idx) := new QuadDobl_Complex_Vectors.Vector'(cf);
          for i in t.dg'range loop
            xp(i) := integer32(t.dg(i));
          end loop;
          res.xps(idx) := new Standard_Integer_Vectors.Vector'(xp);
        end;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cst := new QuadDobl_Complex_Vectors.Vector'(0..deg => cmplxzero);
    res.cst(0) := cst;
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(zero));
    return res;
  end Make_Convolution_Circuit;

  function Make_Convolution_Circuit
             ( p : Standard_CSeries_Polynomials.Poly )
             return Standard_Speelpenning_Convolutions.Circuit is

    use Standard_Speelpenning_Convolutions;

    dim : constant integer32
        := integer32(Standard_CSeries_Polynomials.Number_of_Unknowns(p));
    zero : Standard_CSeries_Polynomials.Degrees
         := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    cst : constant Standard_Complex_Series.Link_to_Series
        := Standard_CSeries_Polynomials.Coeff(p,zero);
    nbr : constant integer32
        := integer32(Standard_CSeries_Polynomials.Number_of_Terms(p))
           - Nonzero_Constant(cst);
    res : Circuit(nbr,dim,dim-1,dim-2);
    idx,deg : integer32 := 0;

    procedure Visit_Term ( t : in Standard_CSeries_Polynomials.Term;
                           continue : out boolean ) is

      use Standard_CSeries_Polynomials;

    begin
      if not Is_Zero(t.dg) then
        idx := idx + 1;
        res.cff(idx) := new Standard_Complex_Vectors.Vector'(t.cf.cff);
        declare
          xp : Standard_Integer_Vectors.Vector(1..dim);
        begin
          for i in t.dg'range loop
            xp(i) := integer32(t.dg(i));
          end loop;
          res.xps(idx) := new Standard_Integer_Vectors.Vector'(xp);
        end;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cst := new Standard_Complex_Vectors.Vector'(cst.cff);
    deg := res.cst'last;
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(zero));
    return res;
  end Make_Convolution_Circuit;

  function Make_Convolution_Circuit
             ( p : DoblDobl_CSeries_Polynomials.Poly )
             return DoblDobl_Speelpenning_Convolutions.Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    dim : constant integer32
        := integer32(DoblDobl_CSeries_Polynomials.Number_of_Unknowns(p));
    zero : DoblDobl_CSeries_Polynomials.Degrees
         := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    cst : constant DoblDobl_Complex_Series.Link_to_Series
        := DoblDobl_CSeries_Polynomials.Coeff(p,zero);
    nbr : constant integer32
        := integer32(DoblDobl_CSeries_Polynomials.Number_of_Terms(p))
           - Nonzero_Constant(cst);
    res : Circuit(nbr,dim,dim-1,dim-2);
    idx,deg : integer32 := 0;

    procedure Visit_Term ( t : in DoblDobl_CSeries_Polynomials.Term;
                           continue : out boolean ) is

      use DoblDobl_CSeries_Polynomials;

    begin
      if not Is_Zero(t.dg) then
        idx := idx + 1;
        res.cff(idx) := new DoblDobl_Complex_Vectors.Vector'(t.cf.cff);
        declare
          xp : Standard_Integer_Vectors.Vector(1..dim);
        begin
          for i in t.dg'range loop
            xp(i) := integer32(t.dg(i));
          end loop;
          res.xps(idx) := new Standard_Integer_Vectors.Vector'(xp);
        end;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cst := new DoblDobl_Complex_Vectors.Vector'(cst.cff);
    deg := res.cst'last;
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(zero));
    return res;
  end Make_Convolution_Circuit;

  function Make_Convolution_Circuit
             ( p : QuadDobl_CSeries_Polynomials.Poly )
             return QuadDobl_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    dim : constant integer32
        := integer32(QuadDobl_CSeries_Polynomials.Number_of_Unknowns(p));
    zero : QuadDobl_CSeries_Polynomials.Degrees
         := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    cst : constant QuadDobl_Complex_Series.Link_to_Series
        := QuadDobl_CSeries_Polynomials.Coeff(p,zero);
    nbr : constant integer32
        := integer32(QuadDobl_CSeries_Polynomials.Number_of_Terms(p))
           - Nonzero_Constant(cst);
    res : Circuit(nbr,dim,dim-1,dim-2);
    idx,deg : integer32 := 0;

    procedure Visit_Term ( t : in QuadDobl_CSeries_Polynomials.Term;
                           continue : out boolean ) is

      use QuadDobl_CSeries_Polynomials;

    begin
      if not Is_Zero(t.dg) then
        idx := idx + 1;
        res.cff(idx) := new QuadDobl_Complex_Vectors.Vector'(t.cf.cff);
        declare
          xp : Standard_Integer_Vectors.Vector(1..dim);
        begin
          for i in t.dg'range loop
            xp(i) := integer32(t.dg(i));
          end loop;
          res.xps(idx) := new Standard_Integer_Vectors.Vector'(xp);
        end;
      end if;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cst := new QuadDobl_Complex_Vectors.Vector'(cst.cff);
    deg := res.cst'last;
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(zero));
    return res;
  end Make_Convolution_Circuit;

  function Make_Convolution_Circuits
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return Standard_Speelpenning_Convolutions.Circuits is

    use Standard_Speelpenning_Convolutions;

    res : Circuits(p'range);

  begin
    for i in p'range loop
      res(i) := new Circuit'(Make_Convolution_Circuit(p(i),d));
    end loop;
    return res;
  end Make_Convolution_Circuits;

  function Make_Convolution_Circuits
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return DoblDobl_Speelpenning_Convolutions.Circuits is

    use DoblDobl_Speelpenning_Convolutions;

    res : Circuits(p'range);

  begin
    for i in p'range loop
      res(i) := new Circuit'(Make_Convolution_Circuit(p(i),d));
    end loop;
    return res;
  end Make_Convolution_Circuits;

  function Make_Convolution_Circuits
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               d : natural32 )
             return QuadDobl_Speelpenning_Convolutions.Circuits is

    use QuadDobl_Speelpenning_Convolutions;

    res : Circuits(p'range);

  begin
    for i in p'range loop
      res(i) := new Circuit'(Make_Convolution_Circuit(p(i),d));
    end loop;
    return res;
  end Make_Convolution_Circuits;

  function Make_Convolution_Circuits
             ( p : Standard_CSeries_Poly_Systems.Poly_Sys )
             return Standard_Speelpenning_Convolutions.Circuits is

    use Standard_Speelpenning_Convolutions;

    res : Circuits(p'range);

  begin
    for i in p'range loop
      res(i) := new Circuit'(Make_Convolution_Circuit(p(i)));
    end loop;
    return res;
  end Make_Convolution_Circuits;

  function Make_Convolution_Circuits
             ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys )
             return DoblDobl_Speelpenning_Convolutions.Circuits is

    use DoblDobl_Speelpenning_Convolutions;

    res : Circuits(p'range);

  begin
    for i in p'range loop
      res(i) := new Circuit'(Make_Convolution_Circuit(p(i)));
    end loop;
    return res;
  end Make_Convolution_Circuits;

  function Make_Convolution_Circuits
             ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys )
             return QuadDobl_Speelpenning_Convolutions.Circuits is

    use QuadDobl_Speelpenning_Convolutions;

    res : Circuits(p'range);

  begin
    for i in p'range loop
      res(i) := new Circuit'(Make_Convolution_Circuit(p(i)));
    end loop;
    return res;
  end Make_Convolution_Circuits;

  function to_double
	     ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
	     return Standard_Speelpenning_Convolutions.Circuit is

    use DoblDobl_Complex_Vectors;
    use Standard_Speelpenning_Convolutions;

    n : constant integer32 := c.nbr;
    d : constant integer32 := c.dim;
    d1 : constant integer32 := d-1;
    d2 : constant integer32 := d-2;
    deg : constant integer32 := c.cff(1)'last;
    res : Standard_Speelpenning_Convolutions.Circuit(n,d,d1,d2);

  begin
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    res.cff := Varbprec_VecVec_Conversions.dd2d(c.cff);
    if c.cst /= null then
      declare
        cst : constant Standard_Complex_Vectors.Vector(c.cst'range)
            := DoblDobl_Complex_to_Standard(c.cst.all);
      begin
        res.cst := new Standard_Complex_Vectors.Vector'(cst);
      end;
    end if;
    res.forward := Allocate_Coefficients(d1,deg);
    res.backward := Allocate_Coefficients(d2,deg);
    res.cross := Allocate_Coefficients(d2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end to_double;

  function to_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
	     return Standard_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Complex_Vectors;
    use Standard_Speelpenning_Convolutions;

    n : constant integer32 := c.nbr;
    d : constant integer32 := c.dim;
    d1 : constant integer32 := d-1;
    d2 : constant integer32 := d-2;
    deg : constant integer32 := c.cff(1)'last;
    res : Standard_Speelpenning_Convolutions.Circuit(n,d,d1,d2);

  begin
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    res.cff := Varbprec_VecVec_Conversions.qd2d(c.cff);
    if c.cst /= null then
      declare
        cst : constant Standard_Complex_Vectors.Vector(c.cst'range)
            := QuadDobl_Complex_to_Standard(c.cst.all);
      begin
        res.cst := new Standard_Complex_Vectors.Vector'(cst);
      end;
    end if;
    res.forward := Allocate_Coefficients(d1,deg);
    res.backward := Allocate_Coefficients(d2,deg);
    res.cross := Allocate_Coefficients(d2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end to_double;

  function to_double_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
	     return DoblDobl_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Complex_Vectors;
    use DoblDobl_Speelpenning_Convolutions;

    n : constant integer32 := c.nbr;
    d : constant integer32 := c.dim;
    d1 : constant integer32 := d-1;
    d2 : constant integer32 := d-2;
    deg : constant integer32 := c.cff(1)'last;
    res : DoblDobl_Speelpenning_Convolutions.Circuit(n,d,d1,d2);

  begin
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    res.cff := Varbprec_VecVec_Conversions.qd2dd(c.cff);
    if c.cst /= null then
      declare
        cst : constant DoblDobl_Complex_Vectors.Vector(c.cst'range)
            := QuadDobl_Complex_to_DoblDobl(c.cst.all);
      begin
        res.cst := new DoblDobl_Complex_Vectors.Vector'(cst);
      end;
    end if;
    res.forward := Allocate_Coefficients(d1,deg);
    res.backward := Allocate_Coefficients(d2,deg);
    res.cross := Allocate_Coefficients(d2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end to_double_double;

  function to_double
	     ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return Standard_Speelpenning_Convolutions.Link_to_Circuit is

    res : Standard_Speelpenning_Convolutions.Link_to_Circuit;

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      res := new Standard_Speelpenning_Convolutions.Circuit'(to_double(c.all));
    end if;
    return res;
  end to_double;

  function to_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return Standard_Speelpenning_Convolutions.Link_to_Circuit is

    res : Standard_Speelpenning_Convolutions.Link_to_Circuit;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      res := new Standard_Speelpenning_Convolutions.Circuit'(to_double(c.all));
    end if;
    return res;
  end to_double;

  function to_double_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
             return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit is

    res : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      res := new DoblDobl_Speelpenning_Convolutions.Circuit'
                   (to_double_double(c.all));
    end if;
    return res;
  end to_double_double;

  function to_double
	     ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
	     return Standard_Speelpenning_Convolutions.Circuits is

    res : Standard_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := to_double(c(k));
    end loop;
    return res;
  end to_double;

  function to_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
	     return Standard_Speelpenning_Convolutions.Circuits is

    res : Standard_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := to_double(c(k));
    end loop;
    return res;
  end to_double;

  function to_double_double
	     ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
	     return DoblDobl_Speelpenning_Convolutions.Circuits is

    res : DoblDobl_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := to_double_double(c(k));
    end loop;
    return res;
  end to_double_double;

  function to_double
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return Standard_Speelpenning_Convolutions.System is

    res : Standard_Speelpenning_Convolutions.System
            (s.neq,s.neq1,s.dim,s.dim1,s.deg);

    use Standard_Speelpenning_Convolutions;

  begin
    res.crc := to_double(s.crc);
    Standard_Integer_Vectors.Copy(s.mxe,res.mxe);
    res.pwt := Allocate(res.mxe,s.deg);
    res.yd := Allocate_Coefficients(s.dim+1,s.deg);
    res.vy := Linearized_Allocation(s.dim,s.deg);
    res.yv := Allocate_Coefficients(s.dim,s.deg);
    res.vm := Allocate_Coefficients(s.neq,s.dim,s.deg);
    return res;
  end to_double;

  function to_double
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return Standard_Speelpenning_Convolutions.Link_to_System is

    res : Standard_Speelpenning_Convolutions.Link_to_System;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if s /= null then
      declare
        ds : constant Standard_Speelpenning_Convolutions.System
                        (s.neq,s.neq1,s.dim,s.dim+1,s.deg)
           := to_double(s.all);
      begin
        res := new Standard_Speelpenning_Convolutions.System'(ds);
      end;
    end if;
    return res;
  end to_double;

  function to_double_double
             ( s : QuadDobl_Speelpenning_Convolutions.System )
             return DoblDobl_Speelpenning_Convolutions.System is

    res : DoblDobl_Speelpenning_Convolutions.System
            (s.neq,s.neq1,s.dim,s.dim1,s.deg);

    use DoblDobl_Speelpenning_Convolutions;

  begin
    res.crc := to_double_double(s.crc);
    Standard_Integer_Vectors.Copy(s.mxe,res.mxe);
    res.pwt := Allocate(res.mxe,s.deg);
    res.yd := Allocate_Coefficients(s.dim+1,s.deg);
    res.vy := Linearized_Allocation(s.dim,s.deg);
    res.yv := Allocate_Coefficients(s.dim,s.deg);
    res.vm := Allocate_Coefficients(s.neq,s.dim,s.deg);
    return res;
  end to_double_double;

  function to_double_double
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System is

    res : DoblDobl_Speelpenning_Convolutions.Link_to_System;

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if s /= null then
      declare
        ds : constant DoblDobl_Speelpenning_Convolutions.System
                        (s.neq,s.neq1,s.dim,s.dim1,s.deg)
           := to_double_double(s.all);
      begin
        res := new DoblDobl_Speelpenning_Convolutions.System'(ds);
      end;
    end if;
    return res;
  end to_double_double;

end System_Convolution_Circuits;
