with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Supports_of_Polynomial_Systems;
with Standard_Lattice_Polygons;
with Standard_Binomial_Factors;          use Standard_Binomial_Factors;
with Standard_Binomial_Factors_io;       use Standard_Binomial_Factors_io;
with Standard_Puiseux_Certificates;      use Standard_Puiseux_Certificates;
with Standard_Puiseux_Certificates_io;   use Standard_Puiseux_Certificates_io;

procedure Driver_for_Common_Factor is

-- DESCRIPTION :
--   Development of finding common factor of two multivariate polynomials
--   using polyhedral methods.

  procedure Residuals ( f,g : in Poly; r : in List_of_Terms ) is

  -- DESCRITPION :
  --   Evaluates the list of terms at the polynomials f and g
  --   to decide whether there is a binomial factor.

    tmp : List_of_Terms := r;
    vf,vg : Poly;
    tol : constant double_float := 1.0E-8;
    fres,gres : double_float;

  begin
    Initialize_Symbol_Table("t");
    while not Is_Null(tmp) loop
      vf := Evaluate(f,Head_Of(tmp)); fres := Residual(vf);
      vg := Evaluate(g,Head_Of(tmp)); gres := Residual(vg);
      put("value of "); 
      Standard_Binomial_Factors_io.put(Head_Of(tmp));
     -- put_line(" at f :"); put_line(vf);
     -- put_line(" at g :"); put_line(vg);
      put("Residuals at (f,g) : (");
      put(fres,3); put(","); put(gres,3); put(" )");
      if fres + gres < tol
       then put_line("  Found a binomial factor.");
       else put_line("  No binomial factor found.");
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Residuals;

  function Orders_of_Evaluation 
             ( p : Poly; g : Germ ) return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a pair of lower powers of t, obtained after evaluation
  --   of g.t and g in p.

    res : Standard_Integer_Vectors.Vector(1..2);
    e0 : Poly := Evaluate(p,g.t);
    e1 : Poly := Evaluate(p,g);
    tol : constant double_float := 1.0E-8;

  begin
   -- new_line; put("e0 : "); put_line(e0);
   -- new_line; put("e1 : "); put_line(e1);
    res(1) := Order(e0,tol);
    res(2) := Order(e1,tol);
    Clear(e0); Clear(e1);
    return res;
  end Orders_of_Evaluation;

  procedure Residuals ( f,g : in Poly; r : in List_of_Germs ) is

  -- DESCRIPTION :
  --   To verify whether the list of germs (initials and second terms)
  --   properly indicates a regular common factor, we compare the orders
  --   after evaluation (with and without the second term).

    tmp : List_of_Germs := r;
    h : Germ;

  begin
    Initialize_Symbol_Table("t");
    while not Is_Null(tmp) loop
      h := Head_Of(tmp);
      put("Orders of evaluation for tropism"); put(h.t.dg.all);
      put(" at (f,g) : (");
      put(Orders_of_Evaluation(f,h)); put(",");
      put(Orders_of_Evaluation(g,h)); put_line(")");
      tmp := Tail_Of(tmp);
    end loop; 
  end Residuals;

  procedure Common_Factor ( f,g : in Poly; r : in Poly ) is

  -- DESCRIPTION :
  --   Computes initial terms of Puiseux series expansion of common
  --   factors of the polynomials f and g.
  --   In case of randomly generated input, the common factor r should
  --   be provided on input for testing purposes.

    A,B,V,W,N,M,T,rS,rT : List;
    init,binf,rest : List_of_Terms;
    fail : boolean;
    tol : constant double_float := 1.0E-8;
    result : List_of_Germs;

  begin
    Common_Inner_Normals(f,g,false,A,B,V,W,N,M,T,fail);
    if r /= Null_Poly then
      rS := Supports_of_Polynomial_Systems.Create(r);
      rT := Standard_Lattice_Polygons.Inner_Normals(rS);
      put_line("The inner normals to the common factor :"); put(rT);
    end if;
    Initial_Terms(f,g,t,false,20,1.0E-12,1.0E-8,init,fail);
    if Is_Null(init) then
      put_line("No initial terms => no common factor.");
    else
      put_line("The initial terms :"); put(init);
      Residuals(f,g,init);
      Split(f,g,init,1.0E-8,binf,rest);
      if not Is_Null(binf)
       then put_line("Binomial factors : "); put(binf);
      end if;
      if not Is_Null(rest) then
        put("Number of remainder terms is "); put(Length_Of(rest),1);
        put_line(" : "); put(rest);
        Second_Terms(f,g,rest,tol,false,result,fail);
        if not fail then
          put("Number of second terms is "); put(Length_Of(result),1);
          put_line(" : "); put(result);
          Residuals(f,g,result);
        else
          put_line("No second terms found.");
        end if;
      end if;
    end if;
  end Common_Factor;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for two bivariate polynomials
  --   and calls the common factor procedure.

    f,g : Poly;

  begin
    Initialize_Symbol_Table("x","y");
    new_line;
    put_line("Reading two polynomials, f and g, in x and y, end with ;  ...");
    put("Give f : "); get(f);
    put("-> f = "); put(f); new_line;
    put("Give g : "); get(g);
    put("-> g = "); put(g); new_line;
    Common_Factor(f,g,Null_Poly);
  end Main;

begin
  Main;
end Driver_for_Common_Factor;
