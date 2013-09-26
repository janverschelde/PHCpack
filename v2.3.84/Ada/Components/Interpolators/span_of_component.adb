with unchecked_deallocation;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Standard_Complex_Poly_SysFun;
with Multprec_Complex_Poly_SysFun;
with Standard_Linear_Spaces;             use Standard_Linear_Spaces;
with Multprec_Linear_Spaces;             use Multprec_Linear_Spaces;
with Planes_and_Polynomials;             use Planes_and_Polynomials;

package body Span_of_Component is

-- DATA STRUCTURES :

  type Standard_Span_Rep ( nd,d : integer32 ) is record    -- nd = n-d
    frv : Standard_Integer_Vectors.Vector(1..d);
    equ : Standard_Complex_Poly_Systems.Poly_Sys(1..nd);
    lin : Standard_Complex_VecVecs.VecVec(1..nd);
  end record;
  type Multprec_Span_Rep ( nd,d : integer32 ) is record    -- nd = n-d
    frv : Standard_Integer_Vectors.Vector(1..d);
    equ : Multprec_Complex_Poly_Systems.Poly_Sys(1..nd);
    lin : Multprec_Complex_VecVecs.VecVec(1..nd);
  end record;

-- AUXILIARY FOR UNIFORM OUTPUT FORMAT :

  procedure Write_Conclusion
                ( file : in file_type;
                  res : in boolean; eva,tol : in double_float ) is

  -- DESCRIPTION :
  --   This procedure ensures a uniform output format in writing the
  --   conclusion of a stop or membership test.

  begin
    put(file,"Residual test at span of component : "); put(file,eva,3);
    if res
     then put(file,"  <"); put(file,tol,3); put_line(file,"  success.");
     else put(file,"  >"); put(file,tol,3); put_line(file,"  failure.");
    end if;
  end Write_Conclusion;

-- AUXILIARIES TO CREATOR :

  function Create ( L,n : integer32; s : Standard_Sample_List;
                    tol : double_float ) return Standard_Span is

  -- DESCRIPTION :
  --   Returns the span of a nonempty list s with l points of length n.

    res : Standard_Span;
    mat : Standard_Complex_Matrices.Matrix(1..L,1..n);
    vec : Standard_Complex_Matrices.Matrix(1..L-1,1..n);
    rnk : natural32;
    tmp : Standard_Sample_List := s;

  begin
    for i in 1..L loop
      declare
        v : constant Standard_Complex_Vectors.Vector(1..n)
          := Sample_Point(Head_Of(tmp)).v;
      begin
        for j in v'range loop
          mat(i,j) := v(j);
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Rank(L-1,n,mat,tol,vec,rnk);
    if integer32(rnk) < L-1 then
      declare
        res_rep : Standard_Span_Rep(n-integer32(rnk),integer32(rnk));
      begin
        res_rep.frv := Pivots(l-1,n,vec,rnk,tol);
        res_rep.lin := Kernel(mat,vec,rnk,res_rep.frv,tol);
        for i in res_rep.equ'range loop
          res_rep.equ(i) := Hyperplane(res_rep.lin(i).all,tol);
        end loop;
        res := new Standard_Span_Rep'(res_rep);
       end;
    end if;
    return res;
  end Create;

  function Create ( L,n : integer32; size : natural32;
                    s : Multprec_Sample_List;
                    tol : double_float ) return Multprec_Span is

  -- DESCRIPTION :
  --   Returns the span of a nonempty list s with l points of length n.

    res : Multprec_Span;
    mat : Multprec_Complex_Matrices.Matrix(1..L,1..n);
    vec : Multprec_Complex_Matrices.Matrix(1..L-1,1..n);
    rnk : natural32;
    tmp : Multprec_Sample_List := s;

  begin
    for i in 1..L loop
      declare
        v : constant Multprec_Complex_Vectors.Vector(1..n)
          := Sample_Point(Head_Of(tmp)).v;
      begin
        for j in v'range loop
          mat(i,j) := v(j);
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Rank(L-1,n,size,mat,tol,vec,rnk);
    if integer32(rnk) < L-1 then
      declare
        res_rep : Multprec_Span_Rep(n-integer32(rnk),integer32(rnk));
      begin
        res_rep.frv := Pivots(L-1,n,vec,rnk,tol);
        res_rep.lin := Kernel(mat,vec,rnk,res_rep.frv,tol);
        for i in res_rep.equ'range loop
          res_rep.equ(i) := Hyperplane(res_rep.lin(i).all,tol);
        end loop;
        res := new Multprec_Span_Rep'(res_rep);
      end;
    end if;
    return res;
  end Create;

-- CREATORS :

  function Create ( s : Standard_Sample_List; tol : double_float )
                  return Standard_Span is

    res : Standard_Span;
    L,n : integer32;
    first : Standard_Sample;

  begin
    if not Is_Null(s) then
      L := integer32(Length_Of(s));
      first := Head_Of(s);
      n := Number_of_Variables(first) - Number_of_Slices(first);
      res := Create(L,n,s,tol);
    end if;
    return res;
  end Create;

  function Create ( s : Multprec_Sample_List; size : natural32;
                    tol : double_float ) return Multprec_Span is

    res : Multprec_Span;
    L,n : integer32;
    first : Multprec_Sample;

  begin
    if not Is_Null(s) then
      L := integer32(Length_Of(s));
      first := Head_Of(s);
      n := Number_of_Variables(first) - Number_of_Slices(first);
      res := Create(L,n,size,s,tol);
    end if;
    return res;
  end Create;

  function Create ( n,d : natural32; frv : Standard_Integer_Vectors.Vector;
                    equ : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Standard_Span is

    res : Standard_Span;
    res_rep : Standard_Span_Rep(integer32(n-d),integer32(d));

  begin
    res_rep.equ := equ;
    res_rep.frv := frv;
    res := new Standard_Span_Rep'(res_rep);
    return res;
  end Create;

  function Create ( n,d : natural32; frv : Standard_Integer_Vectors.Vector;
                    equ : Multprec_Complex_Poly_Systems.Poly_Sys )
                  return Multprec_Span is

    res : Multprec_Span;
    res_rep : Multprec_Span_Rep(integer32(n-d),integer32(d));

  begin
    res_rep.equ := equ;
    res_rep.frv := frv;
    res := new Multprec_Span_Rep'(res_rep);
    return res;
  end Create;

-- SELECTORS :

  function Empty ( sp : Standard_Span ) return boolean is
  begin
    return (sp = null);
  end Empty;

  function Empty ( sp : Multprec_Span ) return boolean is
  begin
    return (sp = null);
  end Empty;

  function Ambient_Dimension ( sp : Standard_Span ) return natural32 is
  begin
    if sp = null
     then return 0;
     else return natural32(sp.nd + sp.d);
    end if;
  end Ambient_Dimension;

  function Ambient_Dimension ( sp : Multprec_Span ) return natural32 is
  begin
    if sp = null
     then return 0;
     else return natural32(sp.nd + sp.d);
    end if;
  end Ambient_Dimension;

  function Dimension ( sp : Standard_Span ) return natural32 is
  begin
    if sp = null
     then return 0;
     else return natural32(sp.d);
    end if;
  end Dimension;

  function Dimension ( sp : Multprec_Span ) return natural32 is
  begin
    if sp = null
     then return 0;
     else return natural32(sp.d);
    end if;
  end Dimension;

  function Free_Variables ( sp : Standard_Span )
                          return Standard_Integer_Vectors.Vector is
  begin
    return sp.frv;
  end Free_Variables;

  function Free_Variables ( sp : Multprec_Span )
                          return Standard_Integer_Vectors.Vector is
  begin
    return sp.frv;
  end Free_Variables;

  function Equations ( sp : Standard_Span )
                     return Standard_Complex_Poly_Systems.Poly_Sys is
  begin
    return sp.equ;
  end Equations;

  function Equations ( sp : Multprec_Span )
                     return Multprec_Complex_Poly_Systems.Poly_Sys is
  begin
    return sp.equ;
  end Equations;

  function In_Span ( sp : Standard_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean is
  begin
    if sp = null then
      return true;
    else
      declare
        eva : constant Standard_Complex_Vectors.Vector(sp.equ'range)
            := Standard_Complex_Poly_SysFun.Eval(sp.equ,x);
        val : constant double_float := Max_Norm(eva);
      begin
        return (val <= tol);
      end;
    end if;
  end In_Span;

  function In_Span ( file : file_type;
                     sp : Standard_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean is

    res : boolean;

  begin
    if sp = null then
      put_line(file,"Residual test at empty span : success.");
      res := true;
    else
      declare
        eva : constant Standard_Complex_Vectors.Vector(sp.equ'range)
            := Standard_Complex_Poly_SysFun.Eval(sp.equ,x);
        val : constant double_float := Max_Norm(eva);
      begin
        res := (val <= tol);
        Write_Conclusion(file,res,val,tol);
      end;
    end if;
    return res;
  end In_Span;

  function In_Span ( sp : Multprec_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean is
  begin
    if sp = null then
      return true;
    else
      declare
        mpx : Multprec_Complex_Vectors.Vector := Create(x);
        res : constant boolean := In_Span(sp,tol,mpx);
      begin
        Multprec_Complex_Vectors.Clear(mpx);
        return res;
      end;
    end if;
  end In_Span;

  function In_Span ( file : file_type;
                     sp : Multprec_Span; tol : double_float;
                     x : Standard_Complex_Vectors.Vector ) return boolean is

    res : boolean;

  begin
    if sp = null then
      put_line(file,"Residual test at empty span : success.");
      res := true;
    else
      declare
        mpx : Multprec_Complex_Vectors.Vector := Create(x);
      begin
        res := In_Span(file,sp,tol,mpx);
        Multprec_Complex_Vectors.Clear(mpx);
      end;
    end if;
    return res;
  end In_Span;

  function In_Span ( sp : Multprec_Span; tol : double_float;
                     x : Multprec_Complex_Vectors.Vector ) return boolean is
  begin
    if sp = null then
      return true;
    else
      declare
        eva : Multprec_Complex_Vectors.Vector(sp.equ'range)
            := Multprec_Complex_Poly_SysFun.Eval(sp.equ,x);
        val : Floating_Number := Max_Norm(eva);
        fltval : constant double_float := Trunc(val);
        res : constant boolean := (fltval <= tol);
      begin
        Multprec_Complex_Vectors.Clear(eva);
        Clear(val);
        return res;
      end;
    end if;
  end In_Span;

  function In_Span ( file : file_type;
                     sp : Multprec_Span; tol : double_float;
                     x : Multprec_Complex_Vectors.Vector ) return boolean is

    res : boolean;

  begin
    if sp = null then
      put_line(file,"Residual test at empty span : success.");
      res := true;
    else
      declare
        eva : Multprec_Complex_Vectors.Vector(sp.equ'range)
            := Multprec_Complex_Poly_SysFun.Eval(sp.equ,x);
        val : Floating_Number := Max_Norm(eva);
        fltval : constant double_float := Trunc(val);
      begin
        res := (fltval <= tol);
        Write_Conclusion(file,res,fltval,tol);
        Multprec_Complex_Vectors.Clear(eva);
        Clear(val);
      end;
    end if;
    return res;
  end In_Span;

-- SUBSPACE RESTRICTIONS :

  function Restrict_Hyperplane
              ( sp : Standard_Span; L : natural32; tol : double_float;
                hyp : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(0..sp.d);
    p,rp : Standard_Complex_Poly_Systems.Poly_Sys(1..1);

  begin
    p(1) := Hyperplane(hyp);
    rp := Restrict_to_Linear_Space(p,integer32(L),sp.frv,sp.lin,tol);
    res := Polynomial(rp(1));
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_Systems.Clear(rp);
    return res;
  end Restrict_Hyperplane;

  function Restrict_Hyperplane
              ( sp : Multprec_Span; L : natural32; tol : double_float;
                hyp : Multprec_Complex_Vectors.Vector )
              return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(0..sp.d);
    p,rp : Multprec_Complex_Poly_Systems.Poly_Sys(1..1);

  begin
    p(1) := Hyperplane(hyp);
    rp := Restrict_to_Linear_Space(p,integer32(L),sp.frv,sp.lin,tol);
    res := Polynomial(rp(1));
    Multprec_Complex_Poly_Systems.Clear(p);
    Multprec_Complex_Poly_Systems.Clear(rp);
    return res;
  end Restrict_Hyperplane;

  function Restrict_Hyperplane
              ( sp : Standard_Span; L : natural32; tol : double_float;
                hyp : Standard_Complex_VecVecs.VecVec )
              return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(hyp'range);
    p,rp : Standard_Complex_Poly_Systems.Poly_Sys(hyp'range);

  begin
    for i in hyp'range loop
      p(i) := Hyperplane(hyp(i).all);
    end loop;
    rp := Restrict_to_Linear_Space(p,integer32(L),sp.frv,sp.lin,tol);
    for i in hyp'range loop
      res(i) := new Standard_Complex_Vectors.Vector'(Polynomial(rp(i)));
    end loop;
    Standard_Complex_Poly_Systems.Clear(p);
    Standard_Complex_Poly_Systems.Clear(rp);
    return res;
  end Restrict_Hyperplane;

  function Restrict_Hyperplane
              ( sp : Multprec_Span; L : natural32; tol : double_float;
                hyp : Multprec_Complex_VecVecs.VecVec )
              return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(hyp'range);
    p,rp : Multprec_Complex_Poly_Systems.Poly_Sys(hyp'range);

  begin
    for i in hyp'range loop
      p(i) := Hyperplane(hyp(i).all);
    end loop;
    rp := Restrict_to_Linear_Space(p,integer32(L),sp.frv,sp.lin,tol);
    for i in hyp'range loop
      res(i) := new Multprec_Complex_Vectors.Vector'(Polynomial(rp(i)));
    end loop;
    Multprec_Complex_Poly_Systems.Clear(p);
    Multprec_Complex_Poly_Systems.Clear(rp);
    return res;
  end Restrict_Hyperplane;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      tol : double_float;
                      p : Standard_Complex_Poly_Systems.Poly_Sys )
                    return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    res := Restrict_to_Linear_Space(p,integer32(L),sp.frv,sp.lin,tol);
    return res;
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      tol : double_float;
                      p : Multprec_Complex_Poly_Systems.Poly_Sys )
                    return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    res := Restrict_to_Linear_Space(p,integer32(L),sp.frv,sp.lin,tol);
    return res;
  end Restrict;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is
  begin
    return Remove_Variables(x,integer32(L),sp.frv'length+integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is
  begin
    return Remove_Variables(x,integer32(L),sp.frv'length+integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      x : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Vectors.Vector is
  begin
    return Remove_Variables(x,integer32(L),sp.frv'length+integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution )
                    return Standard_Complex_Solutions.Solution is
  begin
    return Remove_Variables(sol,integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution )
                    return Standard_Complex_Solutions.Solution is
  begin
    return Remove_Variables(sol,integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Multprec_Complex_Solutions.Solution )
                    return Multprec_Complex_Solutions.Solution is
  begin
    return Remove_Variables(sol,integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      spt : Standard_Sample ) return Standard_Sample is

    res : Standard_Sample;
    tol : constant double_float := 1.0E-8;
    dim : constant natural32 := Dimension(sp)+L;
    sol : constant Standard_Complex_Solutions.Solution(integer32(dim))
        := Restrict(sp,L,Sample_Point(spt));
    hyp : constant Standard_Complex_VecVecs.VecVec
        := Restrict_Hyperplane(sp,L,tol,Hyperplane_Sections(spt));

  begin
    res := Create(sol,hyp);
    return res;
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      spt : Standard_Sample ) return Standard_Sample is

    res : Standard_Sample;
    tol : constant double_float := 1.0E-8;
    dim : constant natural32 := Dimension(sp)+L;
    sol : constant Standard_Complex_Solutions.Solution(integer32(dim))
        := Restrict(sp,L,Sample_Point(spt));
    hyp : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(spt);
    mphyp,rshyp : Multprec_Complex_VecVecs.VecVec(hyp'range);
    sthyp : Standard_Complex_VecVecs.VecVec(hyp'range);

  begin
    mphyp := Create(hyp);
    rshyp := Restrict_Hyperplane(sp,L,tol,mphyp);
    for i in rshyp'range loop
      sthyp(i) := new Standard_Complex_Vectors.Vector(rshyp(i)'range);
      for j in rshyp(i)'range loop
        sthyp(i)(j) := Round(rshyp(i)(j));
      end loop;
    end loop;
    Multprec_Complex_VecVecs.Clear(mphyp);
    Multprec_Complex_VecVecs.Clear(rshyp);
    res := Create(sol,sthyp);
    return res;
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      spt : Multprec_Sample ) return Multprec_Sample is

    res : Multprec_Sample;
    tol : constant double_float := 1.0E-8;
    dim : constant natural32 := Dimension(sp)+L;
    sol : constant Multprec_Complex_Solutions.Solution(integer32(dim))
        := Restrict(sp,L,Sample_Point(spt));
    hyp : constant Multprec_Complex_VecVecs.VecVec
        := Restrict_Hyperplane(sp,l,tol,Hyperplane_Sections(spt));

  begin
    res := Create(sol,hyp);
    return res;
  end Restrict;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      sps : Standard_Sample_List )
                    return Standard_Sample_List is

    res,res_last : Standard_Sample_List;
    tmp : Standard_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Restrict(sp,L,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sps : Multprec_Sample_List )
                    return Multprec_Sample_List is

    res,res_last : Multprec_Sample_List;
    tmp : Multprec_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Restrict(sp,L,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Restrict;

  function Restrict ( sp : Standard_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution_List )
                    return Standard_Complex_Solutions.Solution_List is
  begin
    return Remove_Variables(sol,integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Standard_Complex_Solutions.Solution_List )
                    return Standard_Complex_Solutions.Solution_List is
  begin
    return Remove_Variables(sol,integer32(L),sp.frv);
  end Restrict;

  function Restrict ( sp : Multprec_Span; L : in natural32;
                      sol : Multprec_Complex_Solutions.Solution_List )
                    return Multprec_Complex_Solutions.Solution_List is
  begin
    return Remove_Variables(sol,integer32(L),sp.frv);
  end Restrict;

-- DESTRUCTORS :

  procedure Clear ( sp : in out Standard_Span ) is

    procedure free is
      new unchecked_deallocation(Standard_Span_Rep,Standard_Span);

  begin
    if sp /= null then
      Standard_Complex_Poly_Systems.Clear(sp.equ);
      free(sp);
    end if;
  end Clear;

  procedure Clear ( sp : in out Multprec_Span ) is

    procedure free is
      new unchecked_deallocation(Multprec_Span_Rep,Multprec_Span);

  begin
    if sp /= null then
      Multprec_Complex_Poly_Systems.Clear(sp.equ);
      free(sp);
    end if;
  end Clear;

end Span_of_Component;
