with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vector_Norms;      use Standard_Complex_Vector_Norms;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Integer32_Transformations; use Standard_Integer32_Transformations;
with Standard_Integer_Transformations_io;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Drivers_for_Static_Lifting;         use Drivers_for_Static_Lifting;
with Black_Box_Solvers;                  use Black_Box_Solvers;
with Series_and_Polynomials;

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

  procedure Initial_Coefficients
              ( file : in file_type;
                p : in Laur_Sys; mic : in Mixed_Cell;
                psub : out Laur_Sys; sols : out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Terms(p,mic.pts.all);
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
              ( p : in Laur_Sys; mic : in Mixed_Cell;
                psub : out Laur_Sys; sols : out Solution_List;
                report : in boolean ) is

    q : Laur_Sys(p'range) := Select_Terms(p,mic.pts.all);
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

  procedure Shift ( p : in out Poly; verbose : in boolean ) is

    mindeg : constant Degrees := Minimal_Degrees(p);
    t : Term;

  begin
    if verbose then
      put("The minimal degrees : ");
      put(Standard_Integer_Vectors.Vector(mindeg.all)); new_line;
      put_line("The polynomial before the shift :");
      put(p); new_line;
    end if;
    t.cf := Create(1.0);
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

  procedure Shift ( p : in out Laur_Sys; verbose : in boolean ) is
  begin
    for i in p'range loop
      Shift(p(i),verbose);
    end loop;
  end Shift;

  procedure Transform_Coordinates
              ( file : in file_type;
                p : in Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Laur_Sys ) is

    i : constant integer32 := v'last;
    t : Transfo := Build_Transfo(v,i);

  begin
    q := Transform(t,p);
    Shift(q,true);
    put(file,"The transformation defined by ");
    put(file,v); put_line(file," :");
    Standard_Integer_Transformations_io.put(file,t);
    put_line(file,"The transformed system : "); put_line(file,q);
  end Transform_Coordinates;

  procedure Transform_Coordinates
              ( p : in Laur_Sys;
                v : in Standard_Integer_Vectors.Vector;
                q : out Laur_Sys; report : in boolean ) is

    i : constant integer32 := v'last;
    t : Transfo := Build_Transfo(v,i);

  begin
    q := Transform(t,p);
    Shift(q,false);
    if report then
      put("The transformation defined by "); put(v); put_line(" :");
      Standard_Integer_Transformations_io.put(t);
      put_line("The transformed system : "); put_line(q);
    end if;
  end Transform_Coordinates;

  function Initial_Residual
              ( p : in Laur_Sys;
                sol : in Standard_Complex_Vectors.Vector )
              return double_float is

    res : double_float := 0.0;
    ext : Standard_Complex_Vectors.Vector(sol'first..sol'last+1);
    eva : Standard_Complex_Vectors.Vector(p'range);

  begin
    ext(sol'range) := sol;
    ext(sol'last+1) := Create(0.0);
    eva := Eval(p,ext);
    -- res := Norm(eva);
    return res;
  end Initial_Residual;

  function Initial_Residuals
              ( file : in file_type; p : in Laur_Sys;
                sols : in Solution_List ) return double_float is

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
              ( p : in Laur_Sys;
                sols : in Solution_List;
                report : in boolean ) return double_float is

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
                p : in Laur_Sys; mic : in Mixed_Cell;
                tsq : out Poly_Sys; sols : out Solution_List ) is

    psub,tvp : Laur_Sys(p'range);
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
              ( p : in Laur_Sys; mic : in Mixed_Cell;
                tsq : out Poly_Sys; sols : out Solution_List;
                report : in boolean ) is

    psub,tvp : Laur_Sys(p'range);
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

  procedure Series
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                xt0 : in Standard_Complex_Vectors.Vector ) is
  begin
    null;
  end Series;

  procedure Series
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                xt0 : in Standard_Complex_Vectors.Vector;
                report : in boolean ) is
  begin
    null;
  end Series;

  procedure Series
              ( file : in file_type;
                p : in Poly_Sys; sols : in Solution_List;
                report : in boolean ) is

    s : Standard_Series_Poly_Systems.Poly_Sys(p'range)
      := Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Series(file,s,ls.v);
      tmp := Tail_Of(tmp);
    end loop;
  end Series;

  procedure Series
              ( p : in Poly_Sys; sols : in Solution_List;
                report : in boolean ) is

    s : Standard_Series_Poly_Systems.Poly_Sys(p'range)
      := Series_and_Polynomials.System_to_Series_System(p,p'last+1);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Series(s,ls.v,report);
      tmp := Tail_Of(tmp);
    end loop;
  end Series;

  procedure Initials
              ( file : in file_type;
                p : in Laur_Sys; mcc : in Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      put(file,"Mixed cell "); put(file,k); put_line(file," :");
      declare
        q : Poly_Sys(p'range);
        qsols : Solution_List;
      begin
        Initial(file,p,mic,q,qsols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Initials;

  procedure Initials
              ( p : in Laur_Sys; mcc : in Mixed_Subdivision;
                report : in boolean ) is

    tmp : Mixed_Subdivision := mcc;
    mic : Mixed_Cell;

  begin
    for k in 1..Length_Of(mcc) loop
      mic := Head_Of(tmp);
      if report then
        put("Mixed cell "); put(k); put_line(" :");
      end if;
      declare
        q : Poly_Sys(p'range);
        qsols : Solution_List;
      begin
        Initial(p,mic,q,qsols,report);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Initials;

end Regular_Solution_Curves_Series;
