with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;          use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;         use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;     use Standard_Laur_Poly_Convertors;
with Standard_Power_Transformations;
with Transformation_of_Supports;

procedure ts_tropelim is

-- DESCRIPTION :
--   Elimination of a variable via a unimodular coordinate transformation,
--   defined by a tropism.

  function Map ( v : Vector ) return Matrix is

  -- DESCRIPTION :
  --   Returns the monomial map defined by the vector v.

    res : Matrix(v'range,v'range);
    piv : constant integer32 := Standard_Power_Transformations.Pivot(v); 

  begin
    if piv = v'first then
      if v(piv) = 1 then
        for j in v'range loop
          res(v'first,j) := v(j);
        end loop;
        for i in res'first(1)+1..res'last(1) loop
          for j in res'range(2) loop
            if i = j
             then res(i,j) := 1;
             else res(i,j) := 0;
            end if;
          end loop;
        end loop;
      end if;
    else
      res := Standard_Power_Transformations.Eliminate(v,piv);
    end if;
    return res;
  end Map;

  procedure Transform ( p : in Laur_Sys; v : in Vector ) is

    elm : constant Matrix(v'range,v'range) := Map(v);
    trp : Laur_Sys(p'range);
    q : Poly_Sys(p'range);

  begin
    new_line;
    put_line("An eliminating power transformation :"); put(elm);
    trp := Transformation_of_Supports.Transform(p,elm);
    new_line;
    put_line("The transformed system :"); put(trp);
    q := Laurent_to_Polynomial_System(trp);
    put_line("The polynomial system :"); put(q);
  end Transform;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a Laurent polynomial system.

    lp : Link_to_Laur_Sys;
    nv : integer32;
    tp : Link_to_Vector;

  begin
    new_line;
    put_line("Transforming a system, given a tropism.");
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of variables : "); put(nv,1); new_line;
    tp := new Vector'(1..nv => 0);
    new_line;
    put_line("Reading the tropism ...");
    put("Give "); put(nv,1); put(" integers : ");
    for i in tp'range loop
      get(tp(i));
    end loop;
    put("The tropism : "); put(tp); new_line;
    Transform(lp.all,tp.all);
  end Main;

begin
  Main;
end ts_tropelim;
