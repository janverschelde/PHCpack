with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with File_Scanning;                     use File_Scanning;
with Greeting_Banners;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_System_Readers;

procedure compsolve
            ( nt : in natural32; infilename,outfilename : in string ) is

  procedure Main is

    infile : file_type;
    append_sols : boolean := false;
    p : Link_to_Poly_Sys;
    q : Link_to_Laur_Sys;

  begin
    Standard_System_Readers.Read_System(infile,infilename,q);
    if q = null then
      put_line(Greeting_Banners.welcome & ".");
      put("Numerical irreducible decomposition solver");
      if nt = 0
       then put(", no tasking");
       else put(", with "); put(nt,1); put(" tasks");
      end if;
      put_line(", in double precision.");
      new_line; get(q);
    else
      Scan_and_Skip(infile,"SOLUTIONS",append_sols);
      append_sols := not append_sols;
      close(infile);
    end if;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(q.all) then
      put_line("compute the numerical irreducible dcomposition ...");
      -- Solve(q,append_sols);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        t : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q.all);
      begin
        p := new Poly_Sys'(t);
        put_line("Computing the numerical irreducible decomposition ...");
        -- Solve(p,append_sols);
      end;
    end if;
  end Main;

begin
  Main;
end compsolve;
