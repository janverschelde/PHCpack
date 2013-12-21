with text_io,integer_io;                use text_io,integer_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Solution_Filters;         use Standard_Solution_Filters;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Witness_Points;                    use Witness_Points;

procedure ts_sqem is

-- DESCRIPTION :
--   Test on finding the index of a Differential Algebraic Equation.

  function Embed ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                   nv,dimsys : natural )
                 return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..nv);
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..nv => 0);
    for i in p'range loop
      res(i) := Add_Variables(p(i),dimsys);
      for j in 1..dimsys loop
        t.cf := Random1;
        t.dg(nv-j+1) := 1;
        Add(res(i),t);
        t.dg(nv-j+1) := 0;
      end loop;
    end loop;
    t.cf := Create(1.0);
    for i in 1..dimsys loop
      t.dg(nv-dimsys+i) := 1;
      res(p'last+i) := Create(t);
      t.dg(nv-dimsys+i) := 0;
    end loop;
    for i in dimj1+1..dimj1+dimsys loop
      declare
        slcff : Standard_Complex_Vectors.Vector(0..nv)
              := Random_Vector(0,nv);
        slice : Poly := Hyperplane(slcff);
      begin
        res(i) := slice;
      end;
    end loop;
    return res;
  end Embed;

  procedure Embed ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    nbequ : constant natural := p'length;
    nbunk : constant natural := Number_of_Unknowns(p(p'first));
    diff : constant integer := nbunk - nbequ;
    dim : natural;
    ans : character;
    file : file_type;

  begin
    put("Number of equations : "); put(nbequ,1); new_line;
    put("Number of unknowns  : "); put(nbunk,1); new_line;
    if diff < 0
     then put_line("The system is overdetermined");
     else put("I am assuming the top dimension is at least ");
          put(diff,1); put_line("...");
          put("Type "); put(dimsys,1); put(" or higher number : ");
          get(dim);
          Add_Embed_Symbols(dim);
          nv := nbunk + dim;
          declare
            ep : Standard_Complex_Poly_Systems.Poly_Sys(1..nv);
          begin
            ep := Embed(p,nv,dim);
            put_line("The squared and embedded polynomial system : ");
            put_line(ep);
            put("Do you want to save the embedded system on file ? (y/n) ");
            Ask_Yes_or_No(ans);
            if ans = 'y'
             then new_line;
                  put_line("Reading the name of the output file.");
                  Read_Name_and_Create_File(file);
                  put_line(file,ep);
            end if;
          end;
    end if;
  end Embed;

  procedure Main is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Slicing and Embedding a General Polynomial System.");
    new_line;
    get(lp);
    new_line;
    Embed(lp.all);
  end Main;

begin
  Main;
end ts_sqem;
