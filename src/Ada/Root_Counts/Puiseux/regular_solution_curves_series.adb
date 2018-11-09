with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vector_Norms;      use Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with DoblDobl_Complex_Laurentials_io;    use DoblDobl_Complex_Laurentials_io;
with QuadDobl_Complex_Laurentials_io;    use QuadDobl_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_SysFun;       use DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Laur_Poly_Convertors;      use DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_SysFun;       use QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Laur_Poly_Convertors;      use QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Black_Box_Solvers;                  use Black_Box_Solvers;
with Standard_Complex_Series;
with DoblDobl_Complex_Series;
with QuadDobl_Complex_Series;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Standard_Newton_Matrix_Series;
with DoblDobl_Newton_Matrix_Series;
with QuadDobl_Newton_Matrix_Series;

package body Regular_Solution_Curves_Series is

  procedure Mixed_Cell_Tropisms
              ( report : in boolean;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 ) is

    dim : constant integer32 := sup'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);

  begin
    if report
     then put_line("The supports : "); put(sup);
    end if;
    mix := (mix'range => 1);
    if report then
      Integer_Create_Mixed_Cells(standard_output,dim,mix,false,sup,mcc);
      Integer_Volume_Computation(standard_output,dim,mix,true,sup,mcc,mv);
    else
      Integer_Create_Mixed_Cells(dim,mix,sup,mcc);
      Integer_Volume_Computation(dim,mix,true,sup,mcc,mv);
    end if;
  end Mixed_Cell_Tropisms;

  procedure Mixed_Cell_Tropisms
              ( file : in file_type;
                sup : in out Array_of_Lists;
                mcc : out Mixed_Subdivision;
                mv : out natural32 ) is

    dim : constant integer32 := sup'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);

  begin
    new_line(file);
    put_line(file,"THE SUPPORTS : ");
    put(file,sup);
    mix := (mix'range => 1);
    Integer_Create_Mixed_Cells(file,dim,mix,false,sup,mcc);
    Integer_Volume_Computation(file,dim,mix,true,sup,mcc,mv);
  end Mixed_Cell_Tropisms;

  function Tropisms ( mcc : Mixed_Subdivision; mixvol : natural32 )
                    return Standard_Integer_VecVecs.VecVec is

    use Mixed_Volume_Computation;

    len : constant integer32 := integer32(mixvol);
    res : Standard_Integer_VecVecs.VecVec(1..len);
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;
    dim : constant integer32 := Head_Of(mcc).nor'last-1;
    mix : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 1);
    idx : integer32 := 0;
    mv : natural32;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      mv := Mixed_Volume(dim,mix,mic);
      for i in 1..mv loop
        idx := idx + 1;
        res(idx) := new Standard_Integer_Vectors.Vector'(mic.nor.all);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Tropisms;

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    q : constant Standard_Complex_Laur_Systems.Laur_Sys(p'range)
      := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(1.0);
    rc : natural32;

  begin
    put_line(file,"The initial form system with t :"); put(file,q);
    psub := Eval(q,one,idx);
    put_line(file,"The initial form system with t = 1 :"); put(file,psub);
    Solve(psub,false,rc,sols);
    put(file,"Computed ");
    put(file,Length_Of(sols),1); put_line(file," solutions.");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Initial_Coefficients;

  procedure Initial_Coefficients
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                report : in boolean ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    q : constant Standard_Complex_Laur_Systems.Laur_Sys(p'range)
      := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(1.0);
    rc : natural32;

  begin
    if report then
      put_line("The initial form system with t :"); put(q);
    end if;
    psub := Eval(q,one,idx);
    if report then
      put_line("The initial form system with t = 1 :"); put(psub);
    end if;
    Solve(psub,false,rc,sols);
    if report then
      put("Computed "); put(Length_Of(sols),1); put_line(" solutions.");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Initial_Coefficients;

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    q : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range)
      := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(integer32(1));
    rc : natural32;

  begin
    put_line(file,"The initial form system with t :"); put(file,q);
    psub := Eval(q,one,idx);
    put_line(file,"The initial form system with t = 1 :"); put(file,psub);
    Solve(psub,false,rc,sols);
    put(file,"Computed ");
    put(file,Length_Of(sols),1); put_line(file," solutions.");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Initial_Coefficients;

  procedure Initial_Coefficients
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    q : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range)
      := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(integer32(1));
    rc : natural32;

  begin
    if report then
      put_line("The initial form system with t :"); put(q);
    end if;
    psub := Eval(q,one,idx);
    if report then
      put_line("The initial form system with t = 1 :"); put(psub);
    end if;
    Solve(psub,false,rc,sols);
    if report then
      put("Computed "); put(Length_Of(sols),1); put_line(" solutions.");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Initial_Coefficients;

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    q : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range)
      := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(integer32(1));
    rc : natural32;

  begin
    put_line(file,"The initial form system with t :"); put(file,q);
    psub := Eval(q,one,idx);
    put_line(file,"The initial form system with t = 1 :"); put(file,psub);
    Solve(psub,false,rc,sols);
    put(file,"Computed ");
    put(file,Length_Of(sols),1); put_line(file," solutions.");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Initial_Coefficients;

  procedure Initial_Coefficients
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                psub : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    q : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range)
      := Select_Terms(p,mic.pts.all);
    idx : constant integer32 := p'last+1;
    one : constant Complex_Number := Create(integer32(1));
    rc : natural32;

  begin
    if report then
      put_line("The initial form system with t :"); put(q);
    end if;
    psub := Eval(q,one,idx);
    if report then
      put_line("The initial form system with t = 1 :"); put(psub);
    end if;
    Solve(psub,false,rc,sols);
    if report then
      put("Computed "); put(Length_Of(sols),1); put_line(" solutions.");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Initial_Coefficients;

  procedure Shift ( file : in file_type;
                    p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    put(file,"The minimal degrees : ");
    put(file,mindeg.all); new_line(file);
    put_line(file,"The polynomial before the shift :");
    put(file,p); new_line(file);
    t.cf := Standard_Complex_Numbers.Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    put_line(file,"The polynomial after the shift :");
    put(file,p); new_line(file);
  end Shift;

  procedure Shift ( file : in file_type;
                    p : in out DoblDobl_Complex_Laurentials.Poly ) is

    use DoblDobl_Complex_Laurentials;

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    put(file,"The minimal degrees : ");
    put(file,mindeg.all); new_line(file);
    put_line(file,"The polynomial before the shift :");
    put(file,p); new_line(file);
    t.cf := DoblDobl_Complex_Numbers.Create(integer32(1));
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    put_line(file,"The polynomial after the shift :");
    put(file,p); new_line(file);
  end Shift;

  procedure Shift ( file : in file_type;
                    p : in out QuadDobl_Complex_Laurentials.Poly ) is

    use QuadDobl_Complex_Laurentials;

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    put(file,"The minimal degrees : ");
    put(file,mindeg.all); new_line(file);
    put_line(file,"The polynomial before the shift :");
    put(file,p); new_line(file);
    t.cf := QuadDobl_Complex_Numbers.Create(integer32(1));
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    put_line(file,"The polynomial after the shift :");
    put(file,p); new_line(file);
  end Shift;

  procedure Shift ( p : in out Standard_Complex_Laurentials.Poly;
                    verbose : in boolean ) is

    use Standard_Complex_Laurentials;

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    if verbose then
      put("The minimal degrees : "); put(mindeg.all); new_line;
      put_line("The polynomial before the shift :");
      put(p); new_line;
    end if;
    t.cf := Standard_Complex_Numbers.Create(1.0);
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    if verbose then
      put_line("The polynomial after the shift :");
      put(p); new_line;
    end if;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Complex_Laurentials.Poly;
                    verbose : in boolean ) is

    use DoblDobl_Complex_Laurentials;

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    if verbose then
      put("The minimal degrees : "); put(mindeg.all); new_line;
      put_line("The polynomial before the shift :");
      put(p); new_line;
    end if;
    t.cf := DoblDobl_Complex_Numbers.Create(integer32(1));
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    if verbose then
      put_line("The polynomial after the shift :");
      put(p); new_line;
    end if;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Complex_Laurentials.Poly;
                    verbose : in boolean ) is

    use QuadDobl_Complex_Laurentials;

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    if verbose then
      put("The minimal degrees : "); put(mindeg.all); new_line;
      put_line("The polynomial before the shift :");
      put(p); new_line;
    end if;
    t.cf := QuadDobl_Complex_Numbers.Create(integer32(1));
    t.dg := new Standard_Integer_Vectors.Vector(mindeg'range);
    for i in mindeg'range loop
      t.dg(i) := -mindeg(i);
    end loop;
    Mul(p,t);
    if verbose then
      put_line("The polynomial after the shift :");
      put(p); new_line;
    end if;
  end Shift;

  procedure Shift ( file : in file_type;
                    p : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in p'range loop
      Shift(file,p(i));
    end loop;
  end Shift;

  procedure Shift ( file : in file_type;
                    p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in p'range loop
      Shift(file,p(i));
    end loop;
  end Shift;

  procedure Shift ( file : in file_type;
                    p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in p'range loop
      Shift(file,p(i));
    end loop;
  end Shift;

  procedure Shift ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                    verbose : in boolean ) is
  begin
    for i in p'range loop
      Shift(p(i),verbose);
    end loop;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    verbose : in boolean ) is
  begin
    for i in p'range loop
      Shift(p(i),verbose);
    end loop;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    verbose : in boolean ) is
  begin
    for i in p'range loop
      Shift(p(i),verbose);
    end loop;
  end Shift;

  function Transform ( d,v : Standard_Integer_Vectors.Vector )
                     return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(d'range) := d;

  begin
    res(res'last) := 0;
    for k in res'range loop
      res(res'last) := res(res'last) + v(k)*d(k);
    end loop;
    return res;
  end Transform;

  function Transform ( t : Standard_Complex_Laurentials.Term;
                       v : Standard_Integer_Vectors.Vector )
                     return Standard_Complex_Laurentials.Term is

    res : Standard_Complex_Laurentials.Term;
    rdg : constant Standard_Integer_Vectors.Vector(t.dg'range)
        := Transform(t.dg.all,v);

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(rdg);
    return res;
  end Transform;

  function Transform ( t : DoblDobl_Complex_Laurentials.Term;
                       v : Standard_Integer_Vectors.Vector )
                     return DoblDobl_Complex_Laurentials.Term is

    res : DoblDobl_Complex_Laurentials.Term;
    rdg : constant Standard_Integer_Vectors.Vector(t.dg'range)
        := Transform(t.dg.all,v);

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(rdg);
    return res;
  end Transform;

  function Transform ( t : QuadDobl_Complex_Laurentials.Term;
                       v : Standard_Integer_Vectors.Vector )
                     return QuadDobl_Complex_Laurentials.Term is

    res : QuadDobl_Complex_Laurentials.Term;
    rdg : constant Standard_Integer_Vectors.Vector(t.dg'range)
        := Transform(t.dg.all,v);

  begin
    res.cf := t.cf;
    res.dg := new Standard_Integer_Vectors.Vector'(rdg);
    return res;
  end Transform;

  function Transform ( p : Standard_Complex_Laurentials.Poly;
                       v : Standard_Integer_Vectors.Vector )
                     return Standard_Complex_Laurentials.Poly is

    res : Standard_Complex_Laurentials.Poly
        := Standard_Complex_Laurentials.Null_Poly;

    procedure Monomial ( t : in Standard_Complex_Laurentials.Term;
                         c : out boolean ) is
   
      rt : constant Standard_Complex_Laurentials.Term := Transform(t,v);

    begin
      Standard_Complex_Laurentials.Add(res,rt);
      c := true;
    end Monomial;
    procedure Monomials is
      new Standard_Complex_Laurentials.Visiting_Iterator(Monomial);

  begin
    Monomials(p);
    return res;
  end Transform;

  function Transform ( p : DoblDobl_Complex_Laurentials.Poly;
                       v : Standard_Integer_Vectors.Vector )
                     return DoblDobl_Complex_Laurentials.Poly is

    res : DoblDobl_Complex_Laurentials.Poly
        := DoblDobl_Complex_Laurentials.Null_Poly;

    procedure Monomial ( t : in DoblDobl_Complex_Laurentials.Term;
                         c : out boolean ) is
   
      rt : constant DoblDobl_Complex_Laurentials.Term := Transform(t,v);

    begin
      DoblDobl_Complex_Laurentials.Add(res,rt);
      c := true;
    end Monomial;
    procedure Monomials is
      new DoblDobl_Complex_Laurentials.Visiting_Iterator(Monomial);

  begin
    Monomials(p);
    return res;
  end Transform;

  function Transform ( p : QuadDobl_Complex_Laurentials.Poly;
                       v : Standard_Integer_Vectors.Vector )
                     return QuadDobl_Complex_Laurentials.Poly is

    res : QuadDobl_Complex_Laurentials.Poly
        := QuadDobl_Complex_Laurentials.Null_Poly;

    procedure Monomial ( t : in QuadDobl_Complex_Laurentials.Term;
                         c : out boolean ) is
   
      rt : constant QuadDobl_Complex_Laurentials.Term := Transform(t,v);

    begin
      QuadDobl_Complex_Laurentials.Add(res,rt);
      c := true;
    end Monomial;
    procedure Monomials is
      new QuadDobl_Complex_Laurentials.Visiting_Iterator(Monomial);

  begin
    Monomials(p);
    return res;
  end Transform;

  function Transform ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                       v : Standard_Integer_Vectors.Vector )
                     return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Transform(p(k),v);
    end loop;
    return res;
  end Transform;

  function Transform ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                       v : Standard_Integer_Vectors.Vector )
                     return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Transform(p(k),v);
    end loop;
    return res;
  end Transform;

  function Transform ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                       v : Standard_Integer_Vectors.Vector )
                     return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Transform(p(k),v);
    end loop;
    return res;
  end Transform;

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    q := Transform(p,v);
    Shift(q,false);
    put_line(file,"The transformed system : "); put_line(file,q);
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    q := Transform(p,v);
    Shift(q,false);
    put_line(file,"The transformed system : "); put_line(file,q);
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    q := Transform(p,v);
    Shift(q,false);
    put_line(file,"The transformed system : "); put_line(file,q);
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is
  begin
    q := Transform(p,v);
    Shift(q,report);
    if report then
      put_line("The transformed system : "); put_line(q);
    end if;
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is
  begin
    q := Transform(p,v);
    Shift(q,report);
    if report then
      put_line("The transformed system : "); put_line(q);
    end if;
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                report : in boolean ) is
  begin
    q := Transform(p,v);
    Shift(q,report);
    if report then
      put_line("The transformed system : "); put_line(q);
    end if;
  end Transform_Coordinates;

  function Initial_Residual
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sol : in Standard_Complex_Vectors.Vector )
              return double_float is

    res : double_float := 0.0;
    ext : Standard_Complex_Vectors.Vector(sol'first..sol'last+1);
    eva : Standard_Complex_Vectors.Vector(p'range);

  begin
    ext(sol'range) := sol;
    ext(sol'last+1) := Standard_Complex_Numbers.Create(0.0);
    eva := Eval(p,ext);
    res := Norm2(eva);
    return res;
  end Initial_Residual;

  function Initial_Residual
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sol : in DoblDobl_Complex_Vectors.Vector )
              return double_float is

    res : double_double;
    ext : DoblDobl_Complex_Vectors.Vector(sol'first..sol'last+1);
    eva : DoblDobl_Complex_Vectors.Vector(p'range);
    zero : constant double_double := create(0.0);

  begin
    ext(sol'range) := sol;
    ext(sol'last+1) := DoblDobl_Complex_Numbers.Create(zero);
    eva := Eval(p,ext);
    res := Norm2(eva);
    return hi_part(res);
  end Initial_Residual;

  function Initial_Residual
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sol : in QuadDobl_Complex_Vectors.Vector )
              return double_float is

    res : quad_double;
    ext : QuadDobl_Complex_Vectors.Vector(sol'first..sol'last+1);
    eva : QuadDobl_Complex_Vectors.Vector(p'range);
    zero : constant quad_double := create(0.0);

  begin
    ext(sol'range) := sol;
    ext(sol'last+1) := QuadDobl_Complex_Numbers.Create(zero);
    eva := Eval(p,ext);
    res := Norm2(eva);
    return hihi_part(res);
  end Initial_Residual;

  function Initial_Residuals
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List )
              return double_float is

    use Standard_Complex_Solutions;

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      put(file,"At solution "); put(file,k,1);
      put(file," : "); put(file,eva); new_line(file);
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  function Initial_Residuals
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                report : in boolean ) return double_float is

    use Standard_Complex_Solutions;

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      if report then
        put("At solution "); put(k,1);
        put(" : "); put(eva); new_line;
      end if;
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  function Initial_Residuals
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List )
              return double_float is

    use DoblDobl_Complex_Solutions;

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      put(file,"At solution "); put(file,k,1);
      put(file," : "); put(file,eva); new_line(file);
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  function Initial_Residuals
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean ) return double_float is

    use DoblDobl_Complex_Solutions;

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      if report then
        put("At solution "); put(k,1);
        put(" : "); put(eva); new_line;
      end if;
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  function Initial_Residuals
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List )
              return double_float is

    use QuadDobl_Complex_Solutions;

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      put(file,"At solution "); put(file,k,1);
      put(file," : "); put(file,eva); new_line(file);
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  function Initial_Residuals
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean ) return double_float is

    use QuadDobl_Complex_Solutions;

    res : double_float := 0.0;
    eva : double_float;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
 
  begin
    for k in 1..Length_Of(tmp) loop
      ls := Head_Of(tmp);
      eva := Initial_Residual(p,ls.v);
      if report then
        put("At solution "); put(k,1);
        put(" : "); put(eva); new_line;
      end if;
      res := res + eva;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Initial_Residuals;

  procedure Initial
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    psub,tvp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    res : double_float;

  begin
    put(file,"-> pretropism ");
    put(file,mic.nor); new_line(file);
    Initial_Coefficients(file,p,mic,psub,sols);
    Transform_Coordinates(file,p,mic.nor.all,tvp);
    res := Initial_Residuals(file,tvp,sols);
    put(file,"The residual :"); put(file,res,3); new_line(file);
    tsq := Positive_Laurent_Polynomial_System(tvp);
  end Initial;

  procedure Initial
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                report : in boolean ) is

    psub,tvp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    res : double_float;

  begin
    if report then
      put("-> pretropism "); put(mic.nor); new_line;
    end if;
    Initial_Coefficients(p,mic,psub,sols,report);
    Transform_Coordinates(p,mic.nor.all,tvp,report);
    res := Initial_Residuals(tvp,sols,report);
    if report then
      put("The residual :"); put(res,3); new_line;
    end if;
    tsq := Positive_Laurent_Polynomial_System(tvp);
  end Initial;

  procedure Initial
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    psub,tvp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    res : double_float;

  begin
    put(file,"-> pretropism ");
    put(file,mic.nor); new_line(file);
    Initial_Coefficients(file,p,mic,psub,sols);
    Transform_Coordinates(file,p,mic.nor.all,tvp);
    res := Initial_Residuals(file,tvp,sols);
    put(file,"The residual :"); put(file,res,3); new_line(file);
    tsq := Positive_Laurent_Polynomial_System(tvp);
  end Initial;

  procedure Initial
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean ) is

    psub,tvp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    res : double_float;

  begin
    if report then
      put("-> pretropism "); put(mic.nor); new_line;
    end if;
    Initial_Coefficients(p,mic,psub,sols,report);
    Transform_Coordinates(p,mic.nor.all,tvp,report);
    res := Initial_Residuals(tvp,sols,report);
    if report then
      put("The residual :"); put(res,3); new_line;
    end if;
    tsq := Positive_Laurent_Polynomial_System(tvp);
  end Initial;

  procedure Initial
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    psub,tvp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    res : double_float;

  begin
    put(file,"-> pretropism ");
    put(file,mic.nor); new_line(file);
    Initial_Coefficients(file,p,mic,psub,sols);
    Transform_Coordinates(file,p,mic.nor.all,tvp);
    res := Initial_Residuals(file,tvp,sols);
    put(file,"The residual :"); put(file,res,3); new_line(file);
    tsq := Positive_Laurent_Polynomial_System(tvp);
  end Initial;

  procedure Initial
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mic : in Mixed_Cell;
                tsq : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean ) is

    psub,tvp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    res : double_float;

  begin
    if report then
      put("-> pretropism "); put(mic.nor); new_line;
    end if;
    Initial_Coefficients(p,mic,psub,sols,report);
    Transform_Coordinates(p,mic.nor.all,tvp,report);
    res := Initial_Residuals(tvp,sols,report);
    if report then
      put("The residual :"); put(res,3); new_line;
    end if;
    tsq := Positive_Laurent_Polynomial_System(tvp);
  end Initial;

  function Series ( file : file_type;
                    p : Standard_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : Standard_Complex_Vectors.Vector;
                    maxdeg,nit : integer32 )
                  return Standard_Complex_Series_Vectors.Vector is

    use Standard_Newton_Matrix_Series;

    res : Standard_Complex_Series_Vectors.Vector(p'range);
    info : integer32;
    deg : integer32 := 1; -- doubles in every step

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Series.Create(xt0(i));
    end loop;
    LU_Newton_Steps(file,p,deg,maxdeg,nit,res,info);
    put_line(file,"The solution series : ");
    Complex_Series_and_Polynomials_io.put(file,res);
    return res;
  end Series;

  function Series ( p : Standard_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : Standard_Complex_Vectors.Vector;
                    maxdeg,nit : integer32; report : boolean )
                  return Standard_Complex_Series_Vectors.Vector is

    use Standard_Newton_Matrix_Series;

    res : Standard_Complex_Series_Vectors.Vector(p'range);
    info : integer32;
    deg : integer32 := 1; -- doubles in every step

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Series.Create(xt0(i));
    end loop;
    if report then
      LU_Newton_Steps(standard_output,p,deg,maxdeg,nit,res,info);
      put_line("The solution series : ");
      Complex_Series_and_Polynomials_io.put(res);
    else
      LU_Newton_Steps(p,deg,maxdeg,nit,res,info);
    end if;
    return res;
  end Series;

  function Series ( file : file_type;
                    p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : DoblDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32 )
                  return DoblDobl_Complex_Series_Vectors.Vector is

    use DoblDobl_Newton_Matrix_Series;

    res : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    info : integer32;
    deg : integer32 := 1; -- doubles in every step

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Series.Create(xt0(i));
    end loop;
    LU_Newton_Steps(file,p,deg,maxdeg,nit,res,info);
    put_line(file,"The solution series : ");
    Complex_Series_and_Polynomials_io.put(file,res);
    return res;
  end Series;

  function Series ( p : DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : DoblDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32; report : boolean )
                  return DoblDobl_Complex_Series_Vectors.Vector is

    use DoblDobl_Newton_Matrix_Series;

    res : DoblDobl_Complex_Series_Vectors.Vector(p'range);
    info : integer32;
    deg : integer32 := 1; -- doubles in every step

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Series.Create(xt0(i));
    end loop;
    if report then
      LU_Newton_Steps(standard_output,p,deg,maxdeg,nit,res,info);
      put_line("The solution series : ");
      Complex_Series_and_Polynomials_io.put(res);
    else
      LU_Newton_Steps(p,deg,maxdeg,nit,res,info);
    end if;
    return res;
  end Series;

  function Series ( file : file_type;
                    p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : QuadDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32 )
                  return QuadDobl_Complex_Series_Vectors.Vector is

    use QuadDobl_Newton_Matrix_Series;

    res : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    info : integer32;
    deg : integer32 := 1; -- doubles in every step

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Series.Create(xt0(i));
    end loop;
    LU_Newton_Steps(file,p,deg,maxdeg,nit,res,info);
    put_line(file,"The solution series : ");
    Complex_Series_and_Polynomials_io.put(file,res);
    return res;
  end Series;

  function Series ( p : QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                    xt0 : QuadDobl_Complex_Vectors.Vector;
                    maxdeg,nit : integer32; report : boolean )
                  return QuadDobl_Complex_Series_Vectors.Vector is

    use QuadDobl_Newton_Matrix_Series;

    res : QuadDobl_Complex_Series_Vectors.Vector(p'range);
    info : integer32;
    deg : integer32 := 1; -- doubles in every step

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Series.Create(xt0(i));
    end loop;
    if report then
      LU_Newton_Steps(standard_output,p,deg,maxdeg,nit,res,info);
      put_line("The solution series : ");
      Complex_Series_and_Polynomials_io.put(res);
    else
      LU_Newton_Steps(p,deg,maxdeg,nit,res,info);
    end if;
    return res;
  end Series;

  function Series ( file : file_type;
                    p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32 )
                  return Standard_Complex_Series_VecVecs.VecVec is

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    res : Standard_Complex_Series_VecVecs.VecVec(1..len);
    s : constant Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
      := Complex_Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for k in res'range loop
      ls := Head_Of(tmp);
      declare
        sersol : Standard_Complex_Series_Vectors.Vector(p'range);
      begin
        sersol := Series(file,s,ls.v,maxdeg,nit);
        res(k) := new Standard_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32; report : boolean )
                  return Standard_Complex_Series_VecVecs.VecVec is

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    res : Standard_Complex_Series_VecVecs.VecVec(1..len);
    s : constant Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
      := Complex_Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for k in res'range loop
      ls := Head_Of(tmp);
      declare
        sersol : Standard_Complex_Series_Vectors.Vector(p'range);
      begin
        sersol := Series(s,ls.v,maxdeg,nit,report);
        res(k) := new Standard_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( file : file_type;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32 )
                  return DoblDobl_Complex_Series_VecVecs.VecVec is

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    res : DoblDobl_Complex_Series_VecVecs.VecVec(1..len);
    s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
      := Complex_Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for k in res'range loop
      ls := Head_Of(tmp);
      declare
        sersol : DoblDobl_Complex_Series_Vectors.Vector(p'range);
      begin
        sersol := Series(file,s,ls.v,maxdeg,nit);
        res(k) := new DoblDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32; report : boolean )
                  return DoblDobl_Complex_Series_VecVecs.VecVec is

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    res : DoblDobl_Complex_Series_VecVecs.VecVec(1..len);
    s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
      := Complex_Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for k in res'range loop
      ls := Head_Of(tmp);
      declare
        sersol : DoblDobl_Complex_Series_Vectors.Vector(p'range);
      begin
        sersol := Series(s,ls.v,maxdeg,nit,report);
        res(k) := new DoblDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( file : file_type;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32 )
                  return QuadDobl_Complex_Series_VecVecs.VecVec is

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    res : QuadDobl_Complex_Series_VecVecs.VecVec(1..len);
    s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
      := Complex_Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for k in res'range loop
      ls := Head_Of(tmp);
      declare
        sersol : QuadDobl_Complex_Series_Vectors.Vector(p'range);
      begin
        sersol := Series(file,s,ls.v,maxdeg,nit);
        res(k) := new QuadDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    maxdeg,nit : integer32; report : boolean )
                  return QuadDobl_Complex_Series_VecVecs.VecVec is

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(sols));
    res : QuadDobl_Complex_Series_VecVecs.VecVec(1..len);
    s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
      := Complex_Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for k in res'range loop
      ls := Head_Of(tmp);
      declare
        sersol : QuadDobl_Complex_Series_Vectors.Vector(p'range);
      begin
        sersol := Series(s,ls.v,maxdeg,nit,report);
        res(k) := new QuadDobl_Complex_Series_Vectors.Vector'(sersol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  procedure Concat ( sto : in out Standard_Complex_Series_VecVecs.VecVec;
                     idx : in out integer32;
                     sfrom : in Standard_Complex_Series_VecVecs.VecVec ) is
  begin
    for k in sfrom'range loop
      idx := idx + 1;
      sto(idx) := sfrom(k);
    end loop;
  end Concat;

  procedure Concat ( sto : in out DoblDobl_Complex_Series_VecVecs.VecVec;
                     idx : in out integer32;
                     sfrom : in DoblDobl_Complex_Series_VecVecs.VecVec ) is
  begin
    for k in sfrom'range loop
      idx := idx + 1;
      sto(idx) := sfrom(k);
    end loop;
  end Concat;

  procedure Concat ( sto : in out QuadDobl_Complex_Series_VecVecs.VecVec;
                     idx : in out integer32;
                     sfrom : in QuadDobl_Complex_Series_VecVecs.VecVec ) is
  begin
    for k in sfrom'range loop
      idx := idx + 1;
      sto(idx) := sfrom(k);
    end loop;
  end Concat;

  function Series ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    mcc : in Mixed_Subdivision; mv : in natural32;
                    maxdeg,nit : in integer32 ) 
                  return Standard_Complex_Series_VecVecs.VecVec is

    res : Standard_Complex_Series_VecVecs.VecVec(1..integer32(mv));
    idx : integer32 := 0;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      put(file,"Mixed cell "); put(file,k); put_line(file," :");
      declare
        q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
        qsols : Standard_Complex_Solutions.Solution_List;
      begin
        Initial(file,p,mic,q,qsols);
        declare
          len : constant integer32
              := integer32(Standard_Complex_Solutions.Length_Of(qsols));
          srs : Standard_Complex_Series_VecVecs.VecVec(1..len);
        begin
          srs := Series(file,q,qsols,maxdeg,nit);
          Concat(res,idx,srs);
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32; report : boolean )
                  return Standard_Complex_Series_VecVecs.VecVec is

    res : Standard_Complex_Series_VecVecs.VecVec(1..integer32(mv));
    idx : integer32 := 0;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      if report then
        put("Mixed cell "); put(k); put_line(" :");
      end if;
      declare
        q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
        qsols : Standard_Complex_Solutions.Solution_List;
      begin
        Initial(p,mic,q,qsols,report);
        declare
          len : constant integer32
              := integer32(Standard_Complex_Solutions.Length_Of(qsols));
          srs : Standard_Complex_Series_VecVecs.VecVec(1..len);
        begin
          srs := Series(q,qsols,maxdeg,nit,report);
          Concat(res,idx,srs);
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : in Mixed_Subdivision; mv : in natural32;
                    maxdeg,nit : in integer32 )
                  return DoblDobl_Complex_Series_VecVecs.VecVec is

    res : DoblDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
    idx : integer32 := 0;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      put(file,"Mixed cell "); put(file,k); put_line(file," :");
      declare
        q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
        qsols : DoblDobl_Complex_Solutions.Solution_List;
      begin
        Initial(file,p,mic,q,qsols);
        declare
          len : constant integer32
              := integer32(DoblDobl_Complex_Solutions.Length_Of(qsols));
          srs : DoblDobl_Complex_Series_VecVecs.VecVec(1..len);
        begin
          srs := Series(file,q,qsols,maxdeg,nit);
          Concat(res,idx,srs);
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32; report : boolean )
                  return DoblDobl_Complex_Series_VecVecs.VecVec is

    res : DoblDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
    idx : integer32 := 0;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      if report then
        put("Mixed cell "); put(k); put_line(" :");
      end if;
      declare
        q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
        qsols : DoblDobl_Complex_Solutions.Solution_List;
      begin
        Initial(p,mic,q,qsols,report);
        declare
          len : constant integer32
              := integer32(DoblDobl_Complex_Solutions.Length_Of(qsols));
          srs : DoblDobl_Complex_Series_VecVecs.VecVec(1..len);
        begin
          srs := Series(q,qsols,maxdeg,nit,report);
          Concat(res,idx,srs);
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( file : file_type;
                    p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32 )
                  return QuadDobl_Complex_Series_VecVecs.VecVec is

    res : QuadDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
    idx : integer32 := 0;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      put(file,"Mixed cell "); put(file,k); put_line(file," :");
      declare
        q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
        qsols : QuadDobl_Complex_Solutions.Solution_List;
      begin
        Initial(file,p,mic,q,qsols);
        declare
          len : constant integer32
              := integer32(QuadDobl_Complex_Solutions.Length_Of(qsols));
          srs : QuadDobl_Complex_Series_VecVecs.VecVec(1..len);
        begin
          srs := Series(file,q,qsols,maxdeg,nit);
          Concat(res,idx,srs);
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

  function Series ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    mcc : Mixed_Subdivision; mv : natural32;
                    maxdeg,nit : integer32; report : boolean ) 
                  return QuadDobl_Complex_Series_VecVecs.VecVec is

    res : QuadDobl_Complex_Series_VecVecs.VecVec(1..integer32(mv));
    idx : integer32 := 0;
    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      if report then
        put("Mixed cell "); put(k); put_line(" :");
      end if;
      declare
        q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
        qsols : QuadDobl_Complex_Solutions.Solution_List;
      begin
        Initial(p,mic,q,qsols,report);
        declare
          len : constant integer32
              := integer32(QuadDobl_Complex_Solutions.Length_Of(qsols));
          srs : QuadDobl_Complex_Series_VecVecs.VecVec(1..len);
        begin
          srs := Series(q,qsols,maxdeg,nit,report);
          Concat(res,idx,srs);
        end;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Series;

end Regular_Solution_Curves_Series;
