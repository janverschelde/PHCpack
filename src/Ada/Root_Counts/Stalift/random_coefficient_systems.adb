with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;

package body Random_Coefficient_Systems is

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return Standard_Complex_Polynomials.Poly is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := Standard_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return DoblDobl_Complex_Polynomials.Poly is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := DoblDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return QuadDobl_Complex_Polynomials.Poly is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := QuadDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return Standard_Complex_Laurentials.Poly is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := Standard_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := integer32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return DoblDobl_Complex_Laurentials.Poly is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use DoblDobl_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := DoblDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := integer32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Integer_Vectors.List )
                  return QuadDobl_Complex_Laurentials.Poly is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := QuadDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := integer32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return Standard_Complex_Polynomials.Poly is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := Standard_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return DoblDobl_Complex_Polynomials.Poly is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := DoblDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return QuadDobl_Complex_Polynomials.Poly is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := QuadDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := natural32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return Standard_Complex_Laurentials.Poly is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := Standard_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := integer32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return DoblDobl_Complex_Laurentials.Poly is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    use DoblDobl_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := DoblDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := integer32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Lists_of_Floating_Vectors.List )
                  return QuadDobl_Complex_Laurentials.Poly is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;
    tmp : List := L;
    l2v : Link_to_Vector;
    t : Term;

  begin
    t.dg := new Standard_Integer_Vectors.Vector(1..integer32(n));
    while not Is_Null(tmp) loop
      l2v := Head_Of(tmp);
      t.cf := QuadDobl_Random_Numbers.Random1;
      for i in 1..integer32(n) loop
        t.dg(i) := integer32(l2v(i));
      end loop;
      Add(res,t);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(t);
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in res'range loop
      res(i) := Create(n,L(i));
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,l(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,l(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,l(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( n : natural32; mix : Standard_Integer_Vectors.Vector;
                    L : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..integer32(n));
    ind : integer32 := 0;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        ind := ind+1;
        res(ind) := Create(n,L(i));
      end loop;
    end loop;
    return res;
  end Create;

  procedure Add_Constant ( p : in out Standard_Complex_Polynomials.Poly ) is

    c : constant Standard_Complex_Numbers.Complex_Number
      := Standard_Random_Numbers.Random1;

  begin
    Add_Constant(p,c);
  end Add_Constant;

  procedure Add_Constant ( p : in out DoblDobl_Complex_Polynomials.Poly ) is

    c : constant DoblDobl_Complex_Numbers.Complex_Number
      := DoblDobl_Random_Numbers.Random1;

  begin
    Add_Constant(p,c);
  end Add_Constant;

  procedure Add_Constant ( p : in out QuadDobl_Complex_Polynomials.Poly ) is

    c : constant QuadDobl_Complex_Numbers.Complex_Number
      := QuadDobl_Random_Numbers.Random1;

  begin
    Add_Constant(p,c);
  end Add_Constant;

  procedure Add_Constant ( p : in out Standard_Complex_Laurentials.Poly ) is

    c : constant Standard_Complex_Numbers.Complex_Number
      := Standard_Random_Numbers.Random1;

  begin
    Add_Constant(p,c);
  end Add_Constant;

  procedure Add_Constant ( p : in out DoblDobl_Complex_Laurentials.Poly ) is

    c : constant DoblDobl_Complex_Numbers.Complex_Number
      := DoblDobl_Random_Numbers.Random1;

  begin
    Add_Constant(p,c);
  end Add_Constant;

  procedure Add_Constant ( p : in out QuadDobl_Complex_Laurentials.Poly ) is

    c : constant QuadDobl_Complex_Numbers.Complex_Number
      := QuadDobl_Random_Numbers.Random1;

  begin
    Add_Constant(p,c);
  end Add_Constant;

  procedure Add_Constant ( p : in out Standard_Complex_Polynomials.Poly;
                           c : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers,Standard_Complex_Polynomials;
    t : Term;
    n : constant natural32 := Number_of_Unknowns(p);
    cp : Complex_Number;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    cp := Coeff(p,t.dg);
    if ((REAL_PART(cp) = 0.0) and (IMAG_PART(cp) = 0.0))
     then t.cf := c; Add(p,t);
    end if;
    Clear(t);
  end Add_Constant;

  procedure Add_Constant ( p : in out DoblDobl_Complex_Polynomials.Poly;
                           c : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Polynomials;
    t : Term;
    n : constant natural32 := Number_of_Unknowns(p);
    cp : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    cp := Coeff(p,t.dg);
    if ((REAL_PART(cp) = zero) and (IMAG_PART(cp) = zero))
     then t.cf := c; Add(p,t);
    end if;
    Clear(t);
  end Add_Constant;

  procedure Add_Constant ( p : in out QuadDobl_Complex_Polynomials.Poly;
                           c : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Polynomials;
    t : Term;
    n : constant natural32 := Number_of_Unknowns(p);
    cp : Complex_Number;
    zero : constant quad_double := create(0.0);

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    cp := Coeff(p,t.dg);
    if ((REAL_PART(cp) = zero) and (IMAG_PART(cp) = zero))
     then t.cf := c; Add(p,t);
    end if;
    Clear(t);
  end Add_Constant;

  procedure Add_Constant ( p : in out Standard_Complex_Laurentials.Poly;
                           c : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_Complex_Numbers,Standard_Complex_Laurentials;
    t : Term;
    n : constant natural32 := Number_of_Unknowns(p);
    cp : Complex_Number;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    cp := Coeff(p,t.dg);
    if ((REAL_PART(cp) = 0.0) and (IMAG_PART(cp) = 0.0))
     then t.cf := c; Add(p,t);
    end if;
    Clear(t);
  end Add_Constant;

  procedure Add_Constant ( p : in out DoblDobl_Complex_Laurentials.Poly;
                           c : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Laurentials;
    t : Term;
    n : constant natural32 := Number_of_Unknowns(p);
    cp : Complex_Number;
    zero : constant double_double := create(0.0);

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    cp := Coeff(p,t.dg);
    if ((REAL_PART(cp) = zero) and (IMAG_PART(cp) = zero))
     then t.cf := c; Add(p,t);
    end if;
    Clear(t);
  end Add_Constant;

  procedure Add_Constant ( p : in out QuadDobl_Complex_Laurentials.Poly;
                           c : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Laurentials;
    t : Term;
    n : constant natural32 := Number_of_Unknowns(p);
    cp : Complex_Number;
    zero : constant quad_double := create(0.0);

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    cp := Coeff(p,t.dg);
    if ((REAL_PART(cp) = zero) and (IMAG_PART(cp) = zero))
     then t.cf := c; Add(p,t);
    end if;
    Clear(t);
  end Add_Constant;

  procedure Add_Constants
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Add_Constant(p(i));
    end loop;
  end Add_Constants;

  procedure Add_Constants
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Add_Constant(p(i));
    end loop;
  end Add_Constants;

  procedure Add_Constants
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Add_Constant(p(i));
    end loop;
  end Add_Constants;

  procedure Add_Constants
              ( p : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in p'range loop
      Add_Constant(p(i));
    end loop;
  end Add_Constants;

  procedure Add_Constants
              ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in p'range loop
      Add_Constant(p(i));
    end loop;
  end Add_Constants;

  procedure Add_Constants
              ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in p'range loop
      Add_Constant(p(i));
    end loop;
  end Add_Constants;

end Random_Coefficient_Systems;
