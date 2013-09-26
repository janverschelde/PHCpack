with text_io;                           use text_io;
with File_Scanning;                     use File_Scanning;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Witness_Sets;                      use Witness_Sets;
with Standard_Deflation_Methods;        use Standard_Deflation_Methods;

procedure ts_locdim is

-- DESCRIPTION :
--   Interactive development of a local dimension test.

  function Local_Embedding
             ( p : Poly_Sys; v : Vector; d : integer32 ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns a locally embedded system, adding d randomly chosen
  --   hyperplanes through v.  The system on return shares the first
  --   p'length equations with p, the last d equations are random planes.

    res : Poly_Sys(p'first..p'last+d);
    h : Vector(0..v'last);

  begin
    res(p'range) := p;
    for i in 1..d loop
      h(v'range) := Random_Vector(v'first,v'last);
      h(0) := -h(v'range)*v;
      res(p'last+i) := Hyperplane(h);
    end loop;
    return Add_Embedding(res,natural32(d));
  end Local_Embedding;

  procedure Interactive_Local_Dimension_Test
              ( p : in Poly_Sys; sols : in Solution_List;
                tol : in double_float ) is

  -- DESCRIPTION :
  --   Runs a local dimension test for all of the solutions of p.

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    d : integer32 := 0;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put("testing dimension of solution "); put(i,1); put_line(" ...");
      put("  give the dimension to test : "); get(d);
      declare
        ep : constant Poly_Sys(p'first..p'last+d) := Local_Embedding(p,ls.v,d);
        v : Vector(ls.v'first..ls.v'last+d);
      begin
        v(ls.v'range) := ls.v;
        for i in 1..d loop
          v(ls.v'last+i) := Create(0.0);
        end loop;
        Interactive_Algorithmic_Deflation(standard_output,ep,v,tol);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Interactive_Local_Dimension_Test;

  procedure Main is

    infile : file_type;
    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    found : boolean;
    tol : constant double_float := 1.0E-5;

  begin
    new_line;
    put_line("Local dimension test of a solution of a polynomial system.");
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(infile);
    get(infile,lp);
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found then
      get(infile,sols);
     -- if not Is_Null(sols) then
     --   put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
     -- end if;
    end if;
    if Is_Null(sols) then
      new_line;
      Read(sols);
     -- if not Is_Null(sols) then
     --   put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
     -- end if;
    end if;
    Interactive_Local_Dimension_Test(lp.all,sols,tol);
  end Main;

begin
  Main;
end ts_locdim;
