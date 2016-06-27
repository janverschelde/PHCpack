with Standard_Integer_Numbers_io;           use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;
with Series_and_Polynomials;                use Series_and_Polynomials;

package body Series_and_Polynomials_io is

  procedure get ( s : out Standard_Dense_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out DoblDobl_Dense_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( s : out QuadDobl_Dense_Series.Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( file : in file_type;
                  s : out Standard_Dense_Series.Series ) is

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
                  s : out DoblDobl_Dense_Series.Series ) is

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
                  s : out QuadDobl_Dense_Series.Series ) is

    p : QuadDobl_Complex_Polynomials.Poly;

  begin
    if Symbol_Table.Empty
     then Symbol_Table.Init(1);
    end if;
    QuadDobl_Complex_Polynomials_io.get(file,p);
    s := Polynomial_to_Series(p);
    QuadDobl_Complex_Polynomials.Clear(p);
  end get;

  procedure put ( s : in Standard_Dense_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in DoblDobl_Dense_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( s : in QuadDobl_Dense_Series.Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type;
                  s : in Standard_Dense_Series.Series ) is

    p : Standard_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    Standard_Complex_Polynomials_io.put(file,p);
    Standard_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in DoblDobl_Dense_Series.Series ) is

    p : DoblDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    DoblDobl_Complex_Polynomials_io.put(file,p);
    DoblDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in QuadDobl_Dense_Series.Series ) is

    p : QuadDobl_Complex_Polynomials.Poly := Series_to_Polynomial(s);

  begin
    QuadDobl_Complex_Polynomials_io.put(file,p);
    QuadDobl_Complex_Polynomials.Clear(p);
  end put;

  procedure get ( lv : out Standard_Dense_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

    use Standard_Complex_Polynomials;
    use Standard_Dense_Series;

  begin
    Standard_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      Standard_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new Standard_Dense_Series_Vectors.Vector(lp'range);
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
        lv(i) := s;
      end;
    end loop;
  end get;

  procedure get ( lv : out DoblDobl_Dense_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Dense_Series;

  begin
    DoblDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      DoblDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new DoblDobl_Dense_Series_Vectors.Vector(lp'range);
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
        lv(i) := s;
      end;
    end loop;
  end get;

  procedure get ( lv : out QuadDobl_Dense_Series_Vectors.Link_to_Vector;
                  idx : in integer32 := 1; verbose : in boolean := false ) is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Dense_Series;

  begin
    QuadDobl_Complex_Poly_Systems_io.get(lp);
    if verbose then
      put_line("The series as polynomial system :");
      QuadDobl_Complex_Poly_Systems_io.put(lp.all);
      put("Number of variables : ");
      put(integer32(Number_of_Unknowns(lp(lp'first))),1); new_line;
    end if;
    lv := new QuadDobl_Dense_Series_Vectors.Vector(lp'range);
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
        lv(i) := s;
      end;
    end loop;
  end get;

  procedure put ( v : in Standard_Dense_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in DoblDobl_Dense_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( v : in QuadDobl_Dense_Series_Vectors.Vector ) is
  begin
    put(standard_output,v);
  end put;

  procedure put ( file : in file_type;
                  v : in Standard_Dense_Series_Vectors.Vector ) is
  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      put(file,v(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in DoblDobl_Dense_Series_Vectors.Vector ) is
  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      put(file,v(i));
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  v : in QuadDobl_Dense_Series_Vectors.Vector ) is
  begin
    put(file,v'last,2); put(file,"  ");
    put_line(file,"1");
    for i in v'range loop
      put(file,v(i));
      new_line(file);
    end loop;
  end put;

  procedure get ( p : out Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( p : out QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    get(standard_input,p,idx,verbose);
  end get;

  procedure get ( file : in file_type;
                  p : out Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : Standard_Complex_Polynomials.Poly;

  begin
    Standard_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    Standard_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : DoblDobl_Complex_Polynomials.Poly;

  begin
    DoblDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    DoblDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure get ( file : in file_type;
                  p : out QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : QuadDobl_Complex_Polynomials.Poly;

  begin
    QuadDobl_Complex_Polynomials_io.get(file,q);
    p := Polynomial_to_Series_Polynomial(q,idx,verbose);
    QuadDobl_Complex_Polynomials.Clear(q);
  end get;

  procedure put ( p : in Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( p : in QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,p,idx,verbose);
  end put;

  procedure put ( file : in file_type;
                  p : in Standard_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : Standard_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    Standard_Complex_Polynomials_io.put(file,q);
    Standard_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in DoblDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : DoblDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    DoblDobl_Complex_Polynomials_io.put(file,q);
    DoblDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure put ( file : in file_type;
                  p : in QuadDobl_Series_Polynomials.Poly;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    q : QuadDobl_Complex_Polynomials.Poly
      := Series_Polynomial_to_Polynomial(p,idx,verbose);

  begin
    QuadDobl_Complex_Polynomials_io.put(file,q);
    QuadDobl_Complex_Polynomials.Clear(q);
  end put;

  procedure get ( ls : out Standard_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Standard_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant Standard_Series_Poly_Systems.Poly_Sys(lp'range)
        := Series_and_Polynomials.System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new Standard_Series_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out DoblDobl_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    DoblDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant DoblDobl_Series_Poly_Systems.Poly_Sys(lp'range)
        := Series_and_Polynomials.System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new DoblDobl_Series_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure get ( ls : out QuadDobl_Series_Poly_Systems.Link_to_Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    QuadDobl_Complex_Poly_Systems_io.get(lp);
    declare
      s : constant QuadDobl_Series_Poly_Systems.Poly_Sys(lp'range)
        := Series_and_Polynomials.System_to_Series_System(lp.all,idx,verbose);
    begin
      ls := new QuadDobl_Series_Poly_Systems.Poly_Sys'(s);
    end;
  end get;

  procedure put ( s : in Standard_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( s : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is
  begin
    put(standard_output,s,idx,verbose);
  end put;

  procedure put ( file : in file_type;
                  s : in Standard_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : Standard_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_and_Polynomials.Series_System_to_System(s,idx,verbose);

  begin
    Standard_Complex_Poly_Systems_io.put(file,p);
    Standard_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in DoblDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_and_Polynomials.Series_System_to_System(s,idx,verbose);

  begin
    DoblDobl_Complex_Poly_Systems_io.put(file,p);
    DoblDobl_Complex_Poly_Systems.Clear(p);
  end put;

  procedure put ( file : in file_type;
                  s : in QuadDobl_Series_Poly_Systems.Poly_Sys;
                  idx : in integer32 := 0; verbose : in boolean := false ) is

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(s'range)
      := Series_and_Polynomials.Series_System_to_System(s,idx,verbose);

  begin
    QuadDobl_Complex_Poly_Systems_io.put(file,p);
    QuadDobl_Complex_Poly_Systems.Clear(p);
  end put;

end Series_and_Polynomials_io;
