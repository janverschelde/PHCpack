with unchecked_deallocation;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package body Generic_Jacobian_Matrices is

-- CREATORS :

  function Create ( p : Poly_Sys ) return Jaco_Mat is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    res : Jaco_Mat(p'range,1..integer32(n));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Diff(p(i),j);
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( j : Jaco_Mat ) return Eval_Jaco_Mat is

    res : Eval_Jaco_Mat(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        res(k,l) := Create(j(k,l));
      end loop;
    end loop;
    return res;
  end Create;

  procedure Create ( p : Poly_Sys;
                     j : out Eval_Coeff_Jaco_Mat; m : out Mult_Factors ) is
  
    nb : constant natural32 := Number_of_Unknowns(p(p'first));
    nbk : natural32;

  begin
    for k in p'range loop
      nbk := Number_of_Terms(p(k));
      for l in 1..integer32(nb) loop
        declare
          mkl : Vector(1..integer32(nbk));
        begin
          Diff(p(k),l,j(k,l),mkl);
          m(k,l) := new Vectors.Vector'(mkl);
        end;
      end loop;
    end loop;
  end Create;

-- EVALUATORS :

  function Eval ( j : Jaco_Mat; x : Vector ) return Matrix is

    m : Matrix(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        m(k,l) := Eval(j(k,l),x);
      end loop;
    end loop;
    return m;
  end Eval;

  function Eval ( j : Eval_Jaco_Mat; x : Vector ) return Matrix is

    m : Matrix(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        m(k,l) := Eval(j(k,l),x);
      end loop;
    end loop;
    return m;
  end Eval;

  function Eval ( j : Eval_Coeff_Jaco_Mat; m : Mult_Factors;
                  c : VecVec; x : Vector ) return Matrix is
 
    res : Matrix(j'range(1),j'range(2));
    ck,mkl : Link_to_Vector;

  begin
    for k in j'range(1) loop
      ck := c(k);
      declare
        cm : Vector(ck'range);
      begin
        for ell in j'range(2) loop
          mkl := m(k,ell);
          for i in cm'range loop
            cm(i) := mkl(i)*ck(i);
          end loop;
          res(k,ell) := Eval(j(k,ell),cm,x);
          for i in cm'range loop  -- clean up auxiliary coefficients!
            Clear(cm(i));
          end loop;
        end loop;
      end;
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation(Jaco_Mat,Link_to_Jaco_Mat);
  procedure free is 
    new unchecked_deallocation(Eval_Jaco_Mat,Link_to_Eval_Jaco_Mat);
  procedure free is 
    new unchecked_deallocation(Eval_Coeff_Jaco_Mat,
                               Link_to_Eval_Coeff_Jaco_Mat);
  procedure free is
    new unchecked_deallocation(Mult_Factors,Link_to_Mult_Factors);

  procedure Clear ( j : in out Jaco_Mat ) is
  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        Clear(j(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( j : in out Link_to_Jaco_Mat ) is
  begin
    if j /= null
     then Clear(j.all);
          free(j);
    end if;
  end Clear;

  procedure Shallow_Clear ( j : in out Link_to_Jaco_Mat ) is
  begin
    if j /= null
     then free(j);
    end if;
  end Shallow_Clear;

  procedure Clear ( j : in out Eval_Jaco_Mat ) is
  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        Clear(j(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( j : in out Link_to_Eval_Jaco_Mat ) is
  begin
    if j /= null
     then Clear(j.all); free(j);
    end if;
  end Clear;

  procedure Shallow_Clear ( j : in out Link_to_Eval_Jaco_Mat ) is
  begin
    if j /= null
     then free(j);
    end if;
  end Shallow_Clear;

  procedure Clear ( j : in out Eval_Coeff_Jaco_Mat ) is
  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        Clear(j(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( j : in out Link_to_Eval_Coeff_Jaco_Mat ) is
  begin
    if j /= null
     then Clear(j.all);
          free(j);
    end if;
  end Clear;

  procedure Shallow_Clear ( j : in out Link_to_Eval_Coeff_Jaco_Mat ) is
  begin
    if j /= null
     then free(j);
    end if;
  end Shallow_Clear;

  procedure Clear ( m : in out Mult_Factors ) is
  begin
    for k in m'range(1) loop
      for l in m'range(2) loop
        Clear(m(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( m : in out Link_to_Mult_Factors ) is
  begin
    if m /= null
     then Clear(m.all); free(m);
    end if;
  end Clear;
  
  procedure Shallow_Clear ( m : in out Link_to_Mult_Factors ) is
  begin
    if m /= null
     then free(m);
    end if;
  end Shallow_Clear;

end Generic_Jacobian_Matrices;
