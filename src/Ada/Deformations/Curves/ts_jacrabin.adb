with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Embed_Polynomials;
with Standard_Deflate_Singularities;
with Standard_Deflation_Trees_io;

procedure ts_jacrabin is

-- DESSCRIPTION :
--   Interactive test on the development of the formulation of the
--   polynomial system along the lines of the trick of Rabinowitsch
--   to move singular solutions to infinity.

  function Identity_Matrix
              ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for i in 1..n loop
      for j in 1..n loop
        if i = j
         then res(i,j) := Standard_Complex_Numbers.Create(1.0);
         else res(i,j) := Standard_Complex_Numbers.Create(0.0);
        end if;
      end loop;
    end loop;
    return res;
  end Identity_Matrix;

  procedure Add_Last_Multiplier
              ( p : in out Standard_Complex_Polynomials.Poly;
                d : in integer32 ) is

  -- DESCRIPTION :
  --   Updates the last equation as p*y - 1,
  --   where y is the variable at index d.

    t : Standard_Complex_Polynomials.Term;

  begin
    t.cf := Standard_Complex_Numbers.Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(1..d => 0);
    t.dg(d) := 1;
    Standard_Complex_Polynomials.Mul(p,t);
    t.dg(d) := 0;
    Standard_Complex_Polynomials.Sub(p,t);
    Standard_Complex_Polynomials.Clear(t);
  end Add_Last_Multiplier;

  function Jacobian_Rabinowitsch
              ( p : Standard_Complex_Poly_Systems.Poly_Sys )
              return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system augmented with the Jacobian matrix which
  --   puts singular solutions at infinity.

    nvar : constant natural32
         := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
    dim : constant integer32 := integer32(nvar);
    idm : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Identity_Matrix(dim);
    cff : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    defres : Standard_Complex_Poly_Systems.Poly_Sys
           := Standard_Deflate_Singularities.Deflate(p,nvar,idm,cff);
    res : Standard_Complex_Poly_Systems.Poly_Sys(defres'range);
    resdim : constant integer32 := 2*dim+1;

  begin
    res := Standard_Embed_Polynomials.Add_Variables(defres,1);
    Add_Last_Multiplier(res(res'last),resdim);
    return res;
  end Jacobian_Rabinowitsch;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and then makes
  --   the system augmented with the Jacobian matrix.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    nvar : natural32;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put_line("Your system :"); put(lp.all);
    new_line;
    put_line("The Rabinowitsch trick applied to the Jacobian :");
    nvar := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Standard_Deflation_Trees_io.Add_Multiplier_Symbols(1,nvar);
    Symbol_Table.Enlarge(1);
    declare
      sb : Symbol_Table.Symbol;
      jacrbp : constant Standard_Complex_Poly_Systems.Poly_Sys
             := Jacobian_Rabinowitsch(lp.all);
    begin
      sb := (sb'range => ' ');
      sb(1) := 'y';
      sb(2) := 'r';
      sb(3) := 'b';
      Symbol_Table.add(sb);
      put(jacrbp);
    end;
  end Main;

begin
  Main;
end ts_jacrabin;
