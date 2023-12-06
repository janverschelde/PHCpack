with Standard_Integer_Numbers_io;          use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Polynomials_io;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Polynomials_io;
with PentDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_Systems_io;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Polynomials_io;
with OctoDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_Systems_io;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Polynomials_io;
with DecaDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Poly_Systems_io;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Polynomials_io;
with HexaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Poly_Systems_io;
with Complex_Series_and_Polynomials;       use Complex_Series_and_Polynomials;

package body Complex_Series_and_Polynomials_io is

  procedure get ( s : out Standard_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out DoblDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out TripDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out QuadDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out PentDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out OctoDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out DecaDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out HexaDobl_Complex_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( file : in file_type;
                  s : out Standard_Complex_Series.Series ) is

    p : Standard_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    Standard_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    Standard_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out DoblDobl_Complex_Series.Series ) is

    p : DoblDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    DoblDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    DoblDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out TripDobl_Complex_Series.Series ) is

    p : TripDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    TripDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    TripDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out QuadDobl_Complex_Series.Series ) is

    p : QuadDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    QuadDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    QuadDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out PentDobl_Complex_Series.Series ) is

    p : PentDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    PentDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    PentDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out OctoDobl_Complex_Series.Series ) is

    p : OctoDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    OctoDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    OctoDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out DecaDobl_Complex_Series.Series ) is

    p : DecaDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    DecaDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    DecaDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure get ( file : in file_type;
                  s : out HexaDobl_Complex_Series.Series ) is

    p : HexaDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    HexaDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    HexaDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure put ( s : in Standard_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in DoblDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in TripDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in QuadDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in PentDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in OctoDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in DecaDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in HexaDobl_Complex_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type;
                  s : in Standard_Complex_Series.Series ) is

    p : Standard_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    Standard_Complex_Polynomials_io.put(file,p);
    Standard_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in DoblDobl_Complex_Series.Series ) is

    p : DoblDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    DoblDobl_Complex_Polynomials_io.put(file,p);
    DoblDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in TripDobl_Complex_Series.Series ) is

    p : TripDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    TripDobl_Complex_Polynomials_io.put(file,p);
    TripDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in QuadDobl_Complex_Series.Series ) is

    p : QuadDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    QuadDobl_Complex_Polynomials_io.put(file,p);
    QuadDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in PentDobl_Complex_Series.Series ) is

    p : PentDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    PentDobl_Complex_Polynomials_io.put(file,p);
    PentDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in OctoDobl_Complex_Series.Series ) is

    p : OctoDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    OctoDobl_Complex_Polynomials_io.put(file,p);
    OctoDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in DecaDobl_Complex_Series.Series ) is

    p : DecaDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    DecaDobl_Complex_Polynomials_io.put(file,p);
    DecaDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in HexaDobl_Complex_Series.Series ) is

    p : HexaDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    HexaDobl_Complex_Polynomials_io.put(file,p);
    HexaDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure get ( lv : out Standard_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Standard_Complex_Polynomials;
    use Standard_Complex_Series;

  begin
    Standard_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      Standard_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new Standard_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out DoblDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Series;

  begin
    DoblDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      DoblDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new DoblDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out TripDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use TripDobl_Complex_Polynomials;
    use TripDobl_Complex_Series;

  begin
    TripDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      TripDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new TripDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out QuadDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Series;

  begin
    QuadDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      QuadDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new QuadDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out PentDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use PentDobl_Complex_Polynomials;
    use PentDobl_Complex_Series;

  begin
    PentDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      PentDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new PentDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out OctoDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use OctoDobl_Complex_Polynomials;
    use OctoDobl_Complex_Series;

  begin
    OctoDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      OctoDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new OctoDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out DecaDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use DecaDobl_Complex_Polynomials;
    use DecaDobl_Complex_Series;

  begin
    DecaDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      DecaDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new DecaDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure get ( lv : out HexaDobl_Complex_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use HexaDobl_Complex_Polynomials;
    use HexaDobl_Complex_Series;

  begin
    HexaDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      HexaDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new HexaDobl_Complex_Series_Vectors.Vector(lp'range);
    for i in lp'range loop
      if verbose 
       then put("considering series "); put(i,1); put_line(" ...");
      end if;
      declare
        s : constant Series := Polynomial_to_Series(lp(i),idx);
      begin
        if verbose
         then put_line("The series :"); put(s); new_line;
        end if;
        lv(i) := new Series'(s);
      end;
    end loop;
  end get;

  procedure put ( v : in Standard_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in DoblDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in TripDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in QuadDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in PentDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in OctoDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in DecaDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in HexaDobl_Complex_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type;
                  v : in Standard_Complex_Series_Vectors.Vector ) is

    use Standard_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in DoblDobl_Complex_Series_Vectors.Vector ) is

    use DoblDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in TripDobl_Complex_Series_Vectors.Vector ) is

    use TripDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in QuadDobl_Complex_Series_Vectors.Vector ) is

    use QuadDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in PentDobl_Complex_Series_Vectors.Vector ) is

    use PentDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in OctoDobl_Complex_Series_Vectors.Vector ) is

    use OctoDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in DecaDobl_Complex_Series_Vectors.Vector ) is

    use DecaDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in HexaDobl_Complex_Series_Vectors.Vector ) is

    use HexaDobl_Complex_Series;

  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      if v(i) /= null then
        put(file,v(i).all);
        new_line(file);
      end if;
    end loop;
  end put;

  procedure get ( p : out Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( file : in file_type;
                  p : out Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : Standard_Complex_Polynomials.Poly;

  begin
    Standard_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    Standard_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : DoblDobl_Complex_Polynomials.Poly;

  begin
    DoblDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    DoblDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : TripDobl_Complex_Polynomials.Poly;

  begin
    TripDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    TripDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : QuadDobl_Complex_Polynomials.Poly;

  begin
    QuadDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    QuadDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : PentDobl_Complex_Polynomials.Poly;

  begin
    PentDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    PentDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : OctoDobl_Complex_Polynomials.Poly;

  begin
    OctoDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    OctoDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : DecaDobl_Complex_Polynomials.Poly;

  begin
    DecaDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    DecaDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : HexaDobl_Complex_Polynomials.Poly;

  begin
    HexaDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    HexaDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure put ( p : in Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( file : in file_type;
                  p : in Standard_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : Standard_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    Standard_Complex_Polynomials_io.put(file,q);
    Standard_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in DoblDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : DoblDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    DoblDobl_Complex_Polynomials_io.put(file,q);
    DoblDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in TripDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : TripDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    TripDobl_Complex_Polynomials_io.put(file,q);
    TripDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in QuadDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : QuadDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    QuadDobl_Complex_Polynomials_io.put(file,q);
    QuadDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in PentDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : PentDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    PentDobl_Complex_Polynomials_io.put(file,q);
    PentDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in OctoDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : OctoDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    OctoDobl_Complex_Polynomials_io.put(file,q);
    OctoDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in DecaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : DecaDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    DecaDobl_Complex_Polynomials_io.put(file,q);
    DecaDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in HexaDobl_CSeries_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : HexaDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    HexaDobl_Complex_Polynomials_io.put(file,q);
    HexaDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure get ( ls : out Standard_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Standard_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant Standard_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new Standard_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out DoblDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    DoblDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new DoblDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out TripDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    TripDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant TripDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new TripDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out QuadDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    QuadDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new QuadDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out PentDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PentDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant PentDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new PentDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out OctoDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    OctoDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant OctoDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new OctoDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out DecaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    DecaDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant DecaDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new DecaDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out HexaDobl_CSeries_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    HexaDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant HexaDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
        := System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new HexaDobl_CSeries_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure put ( s : in Standard_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( file : in file_type;
                  s : in Standard_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : Standard_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    Standard_Complex_Poly_Systems_io.put(file,p);
    Standard_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    DoblDobl_Complex_Poly_Systems_io.put(file,p);
    DoblDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : TripDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    TripDobl_Complex_Poly_Systems_io.put(file,p);
    TripDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    QuadDobl_Complex_Poly_Systems_io.put(file,p);
    QuadDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in PentDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : PentDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    PentDobl_Complex_Poly_Systems_io.put(file,p);
    PentDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in OctoDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : OctoDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    OctoDobl_Complex_Poly_Systems_io.put(file,p);
    OctoDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in DecaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : DecaDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    DecaDobl_Complex_Poly_Systems_io.put(file,p);
    DecaDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in HexaDobl_CSeries_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : HexaDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_System_to_System(s,idx,verbose);

  begin
    HexaDobl_Complex_Poly_Systems_io.put(file,p);
    HexaDobl_Complex_Poly_Systems.Clear(p);
  end put;

end Complex_Series_and_Polynomials_io;
