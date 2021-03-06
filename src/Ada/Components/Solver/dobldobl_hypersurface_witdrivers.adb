with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_Polar;    use DoblDobl_Complex_Numbers_Polar;
with DoblDobl_Random_Numbers;           use DoblDobl_Random_Numbers;
with DoblDobl_Random_Vectors;           use DoblDobl_Random_Vectors;
with Standard_Natural_Vectors;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;   use DoblDobl_Complex_Poly_Functions;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_Poly_Laur_Convertors;
with DoblDobl_Embed_Polynomials;        use DoblDobl_Embed_Polynomials;
with DoblDobl_Plane_Representations;    use DoblDobl_Plane_Representations;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Witness_Sets;                      use Witness_Sets;
with Witness_Sets_Formats;              use Witness_Sets_Formats;
with DoblDobl_Hypersurface_Witsets;     use DoblDobl_Hypersurface_Witsets;
with DoblDobl_Hypersurface_Witsets_io;  use DoblDobl_Hypersurface_Witsets_io;

package body DoblDobl_Hypersurface_Witdrivers is

  function Embedded_System 
               ( n : in integer32;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,t : in Vector )
               return Poly_Sys is

    use DoblDobl_Complex_Polynomials;

    dim : constant integer32 := 2*n-1;
    ep : constant Poly := Add_Embedding(p,natural32(n-1));
    eqs : constant DoblDobl_Complex_VecVecs.VecVec(1..n-1) := Equations1(b,v);
    hyp : Poly;
    res : Poly_Sys(1..dim);
    et : Term;

  begin
    res(1) := ep;
    et.cf := Create(integer(1));
    et.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
    for i in eqs'range loop
      hyp := Hyperplane(eqs(i).all);
      res(n+i) := Add_Variables(hyp,natural32(n-1));
      et.dg(n+i) := 1;
      Add(res(n+i),et);
      res(i+1) := Create(et);
      et.dg(n+i) := 0;
    end loop;
    Clear(et);
    return res;
  end Embedded_System;

  procedure Call_Root_Finder ( p : in DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    ep : constant Eval_Poly := Create(p);
    b : constant Vector(1..n) := Random_Vector(1,n);
    v : constant Vector(1..n) := Random_Vector(1,n);
    d : constant integer32 := Degree(p);
    max : constant natural32 := 10*natural32(d);
    t,err,res : Vector(1..d);
    eps : constant double_double := create(1.0E-16);
    fail : boolean;
    ans : character;
    nrm,f : double_double;

    function Eval ( x : Vector ) return Complex_Number is
    begin
      return Eval(ep,x);
    end Eval;

    procedure S_Root_Finder is new Silent_Root_Finder1(Eval);
    procedure R_Root_Finder is new Reporting_Root_Finder0(Eval);

  begin
    new_line;
    put("The polynomial : "); put(p); new_line;
    new_line;
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      R_Root_Finder(Standard_Output,natural32(d),eps,max,fail,b,v,t,nrm);
    else
      S_Root_Finder(natural32(d),eps,max,fail,b,v,t,err,res,nrm);
      put_line("A witness set for the polynomial : ");
      for i in t'range loop
        put(t(i)); put(" : "); f := AbsVal(err(i)); put(f,3);
                   put(" : "); f := AbsVal(res(i)); put(f,3);
        new_line;
      end loop;
      put("The residual norm : "); put(nrm,3); new_line;
    end if;
    if fail
     then put_line("The required accuracy is not reached.");
     else put_line("Witness set computation succeeded.");
          Write_Witness_Set(n,p,b,v,t);
    end if;
  end Call_Root_Finder;

  procedure Call_Root_Finder
               ( file : in file_type;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 output : in boolean;
                 eps : in double_double; fail : out boolean ) is

    use DoblDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    ep : constant Eval_Poly := Create(p);
    b : constant Vector(1..n) := Random_Vector(1,n);
    v : constant Vector(1..n) := Random_Vector(1,n);
    d : constant integer32 := Degree(p);
    max : constant natural32 := 10*natural32(d);
    t,err,res : Vector(1..d);
    nrm : double_double;

    function Eval ( x : Vector ) return Complex_Number is
    begin
      return Eval(ep,x);
    end Eval;

    procedure S_Root_Finder is new Silent_Root_Finder1(Eval);
    procedure R_Root_Finder is new Reporting_Root_Finder0(Eval);

  begin
    if output then
      R_Root_Finder(file,natural32(d),eps,max,fail,b,v,t,nrm);
    else
      S_Root_Finder(natural32(d),eps,max,fail,b,v,t,err,res,nrm);
    end if;
    Write_Witness_Set(file,n,p,b,v,t);
    new_line(file);
    put(file,"The input polynomial : "); put(file,p); new_line(file);
    put(file,"The residual norm : "); put(file,nrm,3); new_line(file);
    if fail then
      put_line(file,"The required accuracy is not reached.");
    else
      put_line(file,"Witness set computation succeeded.");
    end if;
  end Call_Root_Finder;

  procedure Silent_Root_Finder
               ( p : in DoblDobl_Complex_Polynomials.Poly;
                 eps : in double_double; fail : out boolean;
                 e : out Link_to_Poly_Sys; esols : out Solution_List ) is

    use DoblDobl_Complex_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p));
    ep : constant Eval_Poly := Create(p);
    b : constant Vector(1..n) := Random_Vector(1,n);
    v : constant Vector(1..n) := Random_Vector(1,n);
    d : constant integer32 := Degree(p);
    t,err,res : Vector(1..d);
    nrm : double_double;

    function Eval ( x : Vector ) return Complex_Number is
    begin
      return Eval(ep,x);
    end Eval;

    procedure S_Root_Finder is new Silent_Root_Finder1(Eval);

  begin
    S_Root_Finder(natural32(d),eps,natural32(d)*10,fail,b,v,t,err,res,nrm);
    e := new Poly_Sys'(Embedded_System(n,p,b,v,t));
    esols := Embedded_Extrinsic_Solutions(n,b,v,t);
  end Silent_Root_Finder;

  procedure Silent_Root_Finder
               ( p : in DoblDobl_Complex_Laurentials.Poly;
                 eps : in double_double; fail : out boolean;
                 e : out Link_to_Laur_Sys; esols : out Solution_List ) is

    n : constant natural32
      := DoblDobl_Complex_Laurentials.Number_of_Unknowns(p);
    q : DoblDobl_Complex_Polynomials.Poly;
    eq : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    laurent : constant boolean
            := DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(p);

  begin
    if laurent then
      q := DoblDobl_Laur_Poly_Convertors.Laurent_Polynomial_to_Polynomial(p);
    else
      q := DoblDobl_Laur_Poly_Convertors.Positive_Laurent_Polynomial(p);
    end if;
    Silent_Root_Finder(q,eps,fail,eq,esols);
    declare
      use DoblDobl_Poly_Laur_Convertors;
      s : Laur_Sys(eq'range) := Polynomial_to_Laurent_System(eq.all);
    begin
      e := new Laur_Sys'(s);
      if laurent then
        DoblDobl_Complex_Laurentials.Clear(e(1));
        e(1) := Witness_Sets.Add_Embedding(p,natural32(e'last)-n);
      end if;
    end;
    DoblDobl_Complex_Polynomials.Clear(q);
    Clear(eq);
  end Silent_Root_Finder;

end DoblDobl_Hypersurface_Witdrivers;
