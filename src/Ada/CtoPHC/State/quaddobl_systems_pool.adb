with unchecked_deallocation;
with QuadDobl_Complex_Polynomials;      use QuadDobl_Complex_Polynomials;

package body QuadDobl_Systems_Pool is

-- DATA STRUCTURES :

  type Link_to_Array_of_Eval_Poly_Sys is access Array_of_Eval_Poly_Sys;
  type Link_to_Array_of_Jaco_Mat is access Array_of_Jaco_Mat;
  type Link_to_Array_of_Eval_Jaco_Mat is access Array_of_Eval_Jaco_Mat;

-- INTERNAL DATA :

  size_pool : integer32;
  sp : Link_to_Array_of_Poly_Sys;
  ep : Link_to_Array_of_Eval_Poly_Sys;
  jm : Link_to_Array_of_Jaco_Mat;
  jf : Link_to_Array_of_Eval_Jaco_Mat;

-- CREATORS :

  procedure Initialize ( n : in integer32 ) is
  begin
    size_pool := n;
    sp := new Array_of_Poly_Sys(1..n);
    ep := new Array_of_Eval_Poly_Sys(1..n);
    jm := new Array_of_Jaco_Mat(1..n);
    jf := new Array_of_Eval_Jaco_Mat(1..n);
  end Initialize;

  procedure Initialize ( k : in integer32; p : in Poly_Sys ) is
  begin
    if k > 0 and k <= size_pool then 
      declare
        q : Poly_Sys(p'range);
      begin
        Copy(p,q);
        sp(k) := new Poly_Sys'(q);
      end;
    end if;
  end Initialize;

-- CONSTRUCTORS :

  procedure Create ( k : in integer32; p : in Poly_Sys ) is
  begin
    Initialize(k,p);
    Create_Evaluator(k);
    Create_Jacobian_Matrix(k);
    Create_Jacobian_Evaluator(k);
  end Create;

  procedure Create_Evaluator ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool then
      declare
        p : constant Link_to_Poly_Sys := sp(k);
        f : constant Eval_Poly_Sys(p'range) := Create(p.all);
      begin
        ep(k) := new Eval_Poly_Sys'(f);
      end;
    end if;
  end Create_Evaluator;

  procedure Create_Jacobian_Matrix ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool then
      declare
        p : constant Link_to_Poly_Sys := sp(k);
        nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
        m : constant Jaco_Mat(p'range,1..nv) := Create(p.all);
      begin
        jm(k) := new Jaco_Mat'(m);
      end;
    end if;
  end Create_Jacobian_Matrix;

  procedure Create_Jacobian_Evaluator ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool then
      declare
        m : constant Link_to_Jaco_Mat := jm(k);
        f : constant Eval_Jaco_Mat(m'range(1),m'range(2)) := Create(m.all);
      begin
        jf(k) := new Eval_Jaco_Mat'(f);
      end;
    end if;
  end Create_Jacobian_Evaluator;

-- SELECTORS :

  function Size return natural32 is
  begin
    return natural32(size_pool);
  end Size;

  function Retrieve ( k : integer32 ) return Link_to_Poly_Sys is

    res : Link_to_Poly_Sys := null;

  begin
    if k > 0 and k <= size_pool
     then res := sp(k);
    end if;
    return res;
  end Retrieve;

  function Evaluator ( k : integer32 ) return Link_to_Eval_Poly_Sys is

    res : Link_to_Eval_Poly_Sys := null;

  begin
    if k > 0 and k <= size_pool
     then res := ep(k);
    end if;
    return res;
  end Evaluator;

  function Jacobian_Matrix ( k : integer32 ) return Link_to_Jaco_Mat is

    res : Link_to_Jaco_Mat := null;

  begin
    if k > 0 and k <= size_pool
     then res := jm(k);
    end if;
    return res;
  end Jacobian_Matrix;

  function Jacobian_Evaluator ( k : integer32 ) return Link_to_Eval_Jaco_Mat is

    res : Link_to_Eval_Jaco_Mat := null;

  begin
    if k > 0 and k <= size_pool
     then res := jf(k);
    end if;
    return res;
  end Jacobian_Evaluator;

-- DESTRUCTORS :

  procedure Clear_System ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool
     then QuadDobl_Complex_Poly_Systems.Clear(sp(k));
    end if;
  end Clear_System;

  procedure Clear_Evaluator ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool
     then QuadDobl_Complex_Poly_SysFun.Clear(ep(k));
    end if;
  end Clear_Evaluator;

  procedure Clear_Jacobian_Matrix ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool
     then QuadDobl_Complex_Jaco_Matrices.Clear(jm(k));
    end if;
  end Clear_Jacobian_Matrix;

  procedure Clear_Jacobian_Evaluator ( k : in integer32 ) is
  begin
    if k > 0 and k <= size_pool
     then QuadDobl_Complex_Jaco_Matrices.Clear(jf(k));
    end if;
  end Clear_Jacobian_Evaluator;

  procedure Clear ( k : in integer32 ) is
  begin
    Clear_System(k);
    Clear_Evaluator(k);
    Clear_Jacobian_Matrix(k);
    Clear_Jacobian_Evaluator(k);
  end Clear;

  procedure Clear is

    procedure free is 
      new unchecked_deallocation(Array_of_Poly_Sys,Link_to_Array_of_Poly_Sys);
    procedure free is 
      new unchecked_deallocation(Array_of_Eval_Poly_Sys,
                                 Link_to_Array_of_Eval_Poly_Sys);
    procedure free is 
      new unchecked_deallocation(Array_of_Jaco_Mat,Link_to_Array_of_Jaco_Mat);
    procedure free is 
      new unchecked_deallocation(Array_of_Eval_Jaco_Mat,
                                 Link_to_Array_of_Eval_Jaco_Mat);
 
  begin
    if size_pool > 0 then
      for k in 1..size_pool loop
        Clear(k);
      end loop;
    end if;
    free(sp); free(ep); free(jm); free(jf);
    size_pool := 0;
  end Clear;

begin
  size_pool := 0;
end QuadDobl_Systems_Pool;
