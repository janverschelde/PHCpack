with unchecked_deallocation;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Multprec_Complex_Poly_Functions;    use Multprec_Complex_Poly_Functions;
with Standard_Polynomial_Interpolators;  use Standard_Polynomial_Interpolators;
with Multprec_Polynomial_interpolators;  use Multprec_Polynomial_Interpolators;
with Interpolation_Points;               use Interpolation_Points;

package body Interpolation_Filters is

-- DATA STRUCTURES :

  type Standard_Filter_Rep is record
    s,s_last : Standard_Sample_Node_List;            -- interpolation points
    pr : Standard_Projector;                         -- projection operator
    fp : Standard_Complex_Polynomials.Poly;          -- filtering polynomial
  end record;

  type Multprec_Filter_Rep is record
    s,s_last : Multprec_Sample_Node_List;            -- interpolation points
    pr : Multprec_Projector;                         -- projection operator
    fp : Multprec_Complex_Polynomials.Poly;          -- filtering polynomial
  end record;

-- AUXILIARY FOR UNIFORM OUTPUT FORMAT :

  procedure Write_Conclusion
                ( file : in file_type; d : in natural32;
                  res : in boolean; eva,tol : in double_float ) is

  -- DESCRIPTION :
  --   This procedure ensures a uniform output format in writing the
  --   conclusion of a stop or membership test.

  begin
    put(file,"Residual test at filter of degree ");
    put(file,d,1); put(file," : ");
    put(file,eva,3);
    if res
     then put(file,"  <"); put(file,tol,3); put_line(file,"  success.");
     else put(file,"  >"); put(file,tol,3); put_line(file,"  failure.");
    end if;
  end Write_Conclusion;

  procedure Write_Conclusion
                ( file : in file_type; d : in natural32;
                  res : in boolean; eva : Floating_Number;
                  tol : in double_float ) is

  -- DESCRIPTION :
  --   This procedure ensures a uniform output format in writing the
  --   conclusion of a stop or membership test.

  begin
    put(file,"Residual test at filter of degree ");
    put(file,d,1); put(file," : ");
    put(file,eva,3);
    if res
     then put(file,"  <"); put(file,tol,3); put_line(file,"  success.");
     else put(file,"  >"); put(file,tol,3); put_line(file,"  failure.");
    end if;
  end Write_Conclusion;

-- AUXILIARIES TO CREATORS :

  function Interpolate ( file : file_type;
                         s : Standard_Sample_Node_List; d,n,nb : natural32 )
                       return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the interpolating polynomial of degree d in n variables
  --   with nb terms through the interpolation points in the list s.

    res : Standard_Complex_Polynomials.Poly;
    ip : Standard_Complex_Polynomials.Poly := Create(d,n,1);
    vv : Standard_Complex_VecVecs.VecVec(1..integer32(nb));
    tmp : Standard_Sample_Node_List := s;
    snd : Standard_Sample_Node;
    invcond : double_float;

  begin
    for i in 1..integer32(nb) loop
      snd := Head_Of(tmp);
      vv(i) := new Standard_Complex_Vectors.Vector'(Sample_Node(snd));
      tmp := Tail_Of(tmp);
    end loop;
   -- res := Interpolate(ip,vv);
    Interpolate(ip,vv,res,invcond);
    put(file,"Estimate for inverse condition number : ");
    put(file,invcond,3); new_line(file);
   -- Standard_Complex_VecVecs.Clear(vv);
    Standard_Complex_Polynomials.Clear(ip);
    return res;
  end Interpolate;

  function Interpolate ( file : file_type;
                         s : Multprec_Sample_Node_List; d,n,nb : natural32 )
                       return Multprec_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the interpolating polynomial of degree d in n variables
  --   with nb terms through the interpolation points in the list s.

    res : Multprec_Complex_Polynomials.Poly;
    ip : Multprec_Complex_Polynomials.Poly := Create(d,n,1);
    vv : Multprec_Complex_VecVecs.VecVec(1..integer32(nb));
    tmp : Multprec_Sample_Node_List := s;
    snd : Multprec_Sample_Node;
    invcond : Floating_Number;

  begin
    for i in 1..integer32(nb) loop
      snd := Head_Of(tmp);
      vv(i) := new Multprec_Complex_Vectors.Vector'(Sample_Node(snd));
      tmp := Tail_Of(tmp);
    end loop;
   -- res := Interpolate(ip,vv);
    Interpolate(ip,vv,res,invcond);
    put(file,"Estimate for inverse condition number : ");
    put(file,invcond,3); new_line(file);
    Clear(invcond);
   -- Multprec_Complex_VecVecs.Clear(vv);
    Multprec_Complex_Polynomials.Clear(ip);
    return res;
  end Interpolate;

-- CREATORS :

  function Create ( p : Standard_Projector ) return Standard_Filter is

    res : Standard_Filter;
    res_rep : Standard_Filter_Rep;

  begin
    res_rep.pr := p;
    res := new Standard_Filter_Rep'(res_rep);
    return res;
  end Create;

  function Create ( p : Multprec_Projector ) return Multprec_Filter is

    res : Multprec_Filter;
    res_rep : Multprec_Filter_Rep;

  begin
    res_rep.pr := p;
    res := new Multprec_Filter_Rep'(res_rep);
    return res;
  end Create;

  procedure Sample_Update
              ( file : in file_type; f : in out Standard_Filter;
                s : in Standard_Sample_List; d : in natural32 ) is

    n,nb : natural32;

  begin
    if f /= null then
      n := Target_Dimension(f.pr);
      Project(f.pr,s,f.s,f.s_last);
      nb := Standard_Polynomial_Interpolators.Number_of_Terms(d,n)-1;
      if nb <= Length_Of(f.s) then
        Standard_Complex_Polynomials.Clear(f.fp);
        f.fp := Interpolate(file,f.s,d,n,nb);
      end if;
    end if;
  end Sample_Update;

  procedure Sample_Update
              ( file : in file_type; f : in out Multprec_Filter;
                s : in Multprec_Sample_List; d : in natural32 ) is

    n,nb : natural32;

  begin
    if f /= null then
      n := Target_Dimension(f.pr);
      Project(f.pr,s,f.s,f.s_last);
      nb := Multprec_Polynomial_Interpolators.Number_of_Terms(d,n)-1;
      if nb <= Length_Of(f.s) then
        Multprec_Complex_Polynomials.Clear(f.fp);
        f.fp := Interpolate(file,f.s,d,n,nb);
      end if;
    end if;
  end Sample_Update;

  procedure Central_Update                                                                     ( f : in out Standard_Filter; basept : in Standard_Sample;
                 basehyp : in Standard_Complex_Vectors.Vector ) is

    tmp : Standard_Sample_Node_List;

  begin
    if f /= null then
      Update(f.pr,basept,basehyp);
      tmp := f.s;
      while not Is_Null(tmp) loop
        declare
          snd : Standard_Sample_Node := Head_Of(tmp);
        begin
          Update(snd,Project(f.pr,Sample_Vector(snd)));
          Set_Head(tmp,snd);
        end;
        tmp := Tail_Of(tmp);
       end loop;
    end if;
  end Central_Update;

  procedure Central_Update
               ( f : in out Multprec_Filter; basept : in Multprec_Sample;
                 basehyp : in Multprec_Complex_Vectors.Vector ) is

    tmp : Multprec_Sample_Node_List;

  begin
    if f /= null then
      Update(f.pr,basept,basehyp);
      tmp := f.s;
      while not Is_Null(tmp) loop
        declare
          snd : Multprec_Sample_Node := Head_Of(tmp);
        begin
          Update(snd,Project(f.pr,Sample_Vector(snd)));
          Set_Head(tmp,snd);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Central_Update;
   
-- SELECTORS :

  function Dimension ( f : Standard_Filter ) return natural32 is
  begin
    if f = null
     then return 0;
     else return Target_Dimension(f.pr)-1;
    end if;
  end Dimension;

  function Dimension ( f : Multprec_Filter ) return natural32 is
  begin
    if f = null
     then return 0;
     else return Target_Dimension(f.pr)-1;
    end if;
  end Dimension;

  function Degree ( f : Standard_Filter ) return natural32 is
  begin
    if f = null
     then return 0;
     else return natural32(Standard_Complex_Polynomials.Degree(f.fp));
    end if;
  end Degree;

  function Degree ( f : Multprec_Filter ) return natural32 is
  begin
    if f = null
     then return 0;
     else return natural32(Multprec_Complex_Polynomials.Degree(f.fp));
    end if;
  end Degree;

  function Sample_Nodes ( f : Standard_Filter )
                        return Standard_Sample_Node_List is
  begin
    return f.s;
  end Sample_Nodes;

  function Sample_Nodes ( f : Multprec_Filter )
                        return Multprec_Sample_Node_List is
  begin
    return f.s;
  end Sample_Nodes;

  function Projector ( f : Standard_Filter ) return Standard_Projector is
  begin
    return f.pr;
  end Projector;

  function Projector ( f : Multprec_Filter ) return Multprec_Projector is
  begin
    return f.pr;
  end Projector;

  function Centrality ( f : Standard_Filter ) return natural32 is
  begin
    if f = null
     then return 0;
     else return Centrality(f.pr);
    end if;
  end Centrality;

  function Centrality ( f : Multprec_Filter ) return natural32 is
  begin
    if f = null
     then return 0;
     else return Centrality(f.pr);
    end if;
  end Centrality;

  function Interpolator ( f : Standard_Filter )
                        return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly;

  begin
    if f /= null
     then res := f.fp;
    end if;
    return res;
  end Interpolator;

  function Interpolator ( f : Multprec_Filter )
                        return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly;

  begin
    if f /= null
     then res := f.fp;
    end if;
    return res;
  end Interpolator;

  function Evaluate ( f : Standard_Filter;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;

  begin
    if f = null then
      res := Standard_Complex_Numbers.Create(0.0);
    else
      declare
        px : constant Standard_Complex_Vectors.Vector := Project(f.pr,x);
      begin
        res := Eval(f.fp,px);
      end;
    end if;
    return res;
  end Evaluate;

  function Evaluate ( f : Multprec_Filter;
                      x : Standard_Complex_Vectors.Vector )
                    return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;
    mpx : Multprec_Complex_Vectors.Vector(x'range) := Create(x);

  begin
    res := Evaluate(f,mpx);
    Multprec_Complex_Vectors.Clear(mpx);   
    return res;
  end Evaluate;

  function Evaluate ( f : Multprec_Filter;
                      x : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Numbers.Complex_Number is

    res : Multprec_Complex_Numbers.Complex_Number;

  begin
    if f = null then
      res := Multprec_Complex_Numbers.Create(integer(0));
    else
      declare
        td : constant integer32 := integer32(Target_Dimension(f.pr));
        px : Multprec_Complex_Vectors.Vector(1..td) := Project(f.pr,x);
      begin
        res := Eval(f.fp,px);
        Multprec_Complex_Vectors.Clear(px);
      end;
    end if;
    return res;
  end Evaluate;

  function On_Component ( f : Standard_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean is

    eva : double_float;

  begin
    if f = null then
      return true;
    else
      eva := Standard_Complex_Numbers.AbsVal(Evaluate(f,x));
      return (eva <= tol);
    end if;
  end On_Component;

  function On_Component ( file : file_type;
                          f : Standard_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean is

    res : boolean;
    eva : double_float;

  begin
    if f = null then
      put_line(file,"Residual test at empty filter : success.");
      res := true;
    else
      eva := Standard_Complex_Numbers.AbsVal(Evaluate(f,x));
      res := (eva <= tol);
      Write_Conclusion(file,Degree(f),res,eva,tol);
    end if;
    return res;
  end On_Component;

  function On_Component ( f : Multprec_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean is
  begin
    if f = null then
      return true;
    else
      declare
        res : boolean;        
        mpx : Multprec_Complex_Vectors.Vector(x'range) := Create(x);
      begin
        res := On_Component(f,tol,mpx);
        Multprec_Complex_Vectors.Clear(mpx);
        return res;
      end;
    end if;
  end On_Component;

  function On_Component ( file : file_type;
                          f : Multprec_Filter; tol : double_float;
                          x : Standard_Complex_Vectors.Vector ) 
                        return boolean is

    res : boolean;

  begin
    if f = null then
      put_line(file,"Residual test at empty filter : success.");
      res := true;
    else
      declare
        mpx : Multprec_Complex_Vectors.Vector(x'range) := Create(x);
      begin
        res := On_Component(file,f,tol,mpx);
        Multprec_Complex_Vectors.Clear(mpx);
      end;
    end if;
    return res;
  end On_Component;

  function On_Component ( f : Multprec_Filter; tol : double_float;
                          x : Multprec_Complex_Vectors.Vector ) 
                        return boolean is

    res : boolean;

  begin
    if f = null then
      res := true;
    else
      declare
        eva : Multprec_Complex_Numbers.Complex_Number := Evaluate(f,x);
        abseva : Floating_Number := Multprec_Complex_Numbers.AbsVal(eva);
        -- flteva : double_float := Trunc(abseva);
      begin
        res := (abseva < tol);
        Multprec_Complex_Numbers.Clear(eva);
        Multprec_Floating_Numbers.Clear(abseva);
        -- return (flteva <= tol);
      end;
    end if;
    return res;
  end On_Component;

  function On_Component ( file : file_type;
                          f : Multprec_Filter; tol : double_float;
                          x : Multprec_Complex_Vectors.Vector ) 
                        return boolean is

    res : boolean;

  begin
    if f = null then
      put_line("Residual test at empty filter : success.");
      res := true;
    else
      declare
        eva : Multprec_Complex_Numbers.Complex_Number := Evaluate(f,x);
        abseva : Floating_Number := Multprec_Complex_Numbers.AbsVal(eva);
        -- flteva : double_float := Trunc(abseva);
      begin
        res := (abseva < tol);
        -- res := (flteva <= tol);
        -- Write_Conclusion(file,Degree(f),res,flteva,tol);
        Write_Conclusion(file,Degree(f),res,abseva,tol);
        Multprec_Complex_Numbers.Clear(eva);
        Multprec_Floating_Numbers.Clear(abseva);
      end;
    end if;
    return res;
  end On_Component;

-- DESTRUCTORS :

  procedure Shallow_Clear ( f : in out Standard_Filter ) is

    procedure free is 
      new unchecked_deallocation(Standard_Filter_Rep,Standard_Filter);

  begin
    if f /= null
     then Standard_Complex_Polynomials.Clear(f.fp); free(f);
    end if;
  end Shallow_Clear;

  procedure Shallow_Clear ( f : in out Multprec_Filter ) is

    procedure free is
      new unchecked_deallocation(Multprec_Filter_Rep,Multprec_Filter);

  begin
    if f /= null
     then Multprec_Complex_Polynomials.Clear(f.fp); free(f);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( f : in out Standard_Filter ) is
  begin
    if f /= null then
      Deep_Clear(f.s);
      Deep_Clear(f.pr);
      Shallow_Clear(f);
    end if;
  end Deep_Clear;

  procedure Deep_Clear ( f : in out Multprec_Filter ) is
  begin
    if f /= null then
      Deep_Clear(f.s);
      Deep_Clear(f.pr);
      Shallow_Clear(f);
    end if;
  end Deep_Clear;

end Interpolation_Filters;
