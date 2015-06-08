with Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Gradient_Evaluations;
with DoblDobl_Complex_Numbers;
with Double_Double_Vectors;
with DoblDobl_Gradient_Evaluations;
with QuadDobl_Complex_Numbers;
with Quad_Double_Vectors;
with QuadDobl_Gradient_Evaluations;
with Multprec_Complex_Numbers;
with Multprec_Floating_Vectors;
with Multprec_Gradient_Evaluations;

package body VarbPrec_Gradient_Evaluations is

-- PART I : ordinary polynomials

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_float;
               maxng,mindg,rcogd : out double_float ) is

    use Standard_Complex_Numbers;
    use Standard_Gradient_Evaluations;

    numcnd : Standard_Floating_Vectors.Vector(0..x'last);
    val : double_float;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    fxdrc := numcnd(0);
    fxrco := fxnrc/fxdrc;
    maxng := AbsVal(ydx(1));
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      if val > maxng
       then maxng := val;
      end if;
    end loop;
    mindg := numcnd(1);
    for k in 2..numcnd'last loop
      if numcnd(k) < mindg
       then mindg := numcnd(k);
      end if;
    end loop;
    rcogd := maxng/mindg;
  end Gradient_with_Inverse_Condition;

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.VecVec;
               ydx : out DoblDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_double;
               maxng,mindg,rcogd : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Gradient_Evaluations;

    numcnd : Double_Double_Vectors.Vector(0..x'last);
    val : double_double;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    fxdrc := numcnd(0);
    fxrco := fxnrc/fxdrc;
    maxng := AbsVal(ydx(1));
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      if val > maxng
       then maxng := val;
      end if;
    end loop;
    mindg := numcnd(1);
    for k in 2..numcnd'last loop
      if numcnd(k) < mindg
       then mindg := numcnd(k);
      end if;
    end loop;
    rcogd := maxng/mindg;
  end Gradient_with_Inverse_Condition;

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.VecVec;
               ydx : out QuadDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out quad_double;
               maxng,mindg,rcogd : out quad_double ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Gradient_Evaluations;

    numcnd : Quad_Double_Vectors.Vector(0..x'last);
    val : quad_double;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    fxdrc := numcnd(0);
    fxrco := fxnrc/fxdrc;
    maxng := AbsVal(ydx(1));
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      if val > maxng
       then maxng := val;
      end if;
    end loop;
    mindg := numcnd(1);
    for k in 2..numcnd'last loop
      if numcnd(k) < mindg
       then mindg := numcnd(k);
      end if;
    end loop;
    rcogd := maxng/mindg;
  end Gradient_with_Inverse_Condition;

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               maxng,mindg,rcogd : out Floating_Number ) is

    use Multprec_Complex_Numbers;
    use Multprec_Gradient_Evaluations;

    numcnd : Multprec_Floating_Vectors.Vector(0..x'last);
    val : Floating_Number;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    Copy(numcnd(0),fxdrc);
    fxrco := fxnrc/fxdrc;
    maxng := AbsVal(ydx(1));
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      if val > maxng
       then Copy(val,maxng);
      end if;
      Clear(val);
    end loop;
    mindg := numcnd(1);
    for k in 2..numcnd'last loop
      if numcnd(k) < mindg
       then Copy(numcnd(k),mindg);
      end if;
    end loop;
    rcogd := maxng/mindg;
    Multprec_Floating_Vectors.Clear(numcnd);
  end Gradient_with_Inverse_Condition;

-- PART II : ordinary polynomial systems

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in Standard_Complex_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.Array_of_VecVecs;
               ydx : in Standard_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out double_float;
               maxng,mindg,rcogd : out double_float ) is

    kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd : double_float;

  begin
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,maxng,mindg,rcogd);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd);
      if kxnrc > fxnrc
       then fxnrc := kxnrc;
      end if;
      if kxdrc < fxdrc
       then fxdrc := kxdrc;
      end if;
      if kmxng > maxng
       then maxng := kmxng;
      end if;
      if kmndg < mindg
       then mindg := kmndg;
      end if;
    end loop;
    fxrco := fxnrc/fxdrc;
    rcogd := maxng/mindg;
  end Jacobian_with_Inverse_Condition;

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in DoblDobl_Complex_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.Array_of_VecVecs;
               ydx : in DoblDobl_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out double_double;
               maxng,mindg,rcogd : out double_double ) is

    kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd : double_double;

  begin
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,maxng,mindg,rcogd);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd);
      if kxnrc > fxnrc
       then fxnrc := kxnrc;
      end if;
      if kxdrc < fxdrc
       then fxdrc := kxdrc;
      end if;
      if kmxng > maxng
       then maxng := kmxng;
      end if;
      if kmndg < mindg
       then mindg := kmndg;
      end if;
    end loop;
    fxrco := fxnrc/fxdrc;
    rcogd := maxng/mindg;
  end Jacobian_with_Inverse_Condition;

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in QuadDobl_Complex_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.Array_of_VecVecs;
               ydx : in QuadDobl_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out quad_double;
               maxng,mindg,rcogd : out quad_double ) is

    kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd : quad_double;

  begin
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,maxng,mindg,rcogd);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd);
      if kxnrc > fxnrc
       then fxnrc := kxnrc;
      end if;
      if kxdrc < fxdrc
       then fxdrc := kxdrc;
      end if;
      if kmxng > maxng
       then maxng := kmxng;
      end if;
      if kmndg < mindg
       then mindg := kmndg;
      end if;
    end loop;
    fxrco := fxnrc/fxdrc;
    rcogd := maxng/mindg;
  end Jacobian_with_Inverse_Condition;

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in Multprec_Complex_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
               ydx : in Multprec_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               maxng,mindg,rcogd : out Floating_Number ) is

    kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd : Floating_Number;

  begin
    Clear(fxnrc); Clear(fxdrc); Clear(fxrco);
    Clear(maxng); Clear(mindg); Clear(rcogd);
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,maxng,mindg,rcogd);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kxnrc,kxdrc,kxrco,kmxng,kmndg,krcgd);
      if kxnrc > fxnrc
       then Copy(kxnrc,fxnrc);
      end if;
      if kxdrc < fxdrc
       then Copy(kxdrc,fxdrc);
      end if;
      if kmxng > maxng
       then Copy(kmxng,maxng);
      end if;
      if kmndg < mindg
       then Copy(kmndg,mindg);
      end if;
      Clear(kxnrc); Clear(kxdrc); Clear(kxrco);
      Clear(maxng); Clear(mindg); Clear(rcogd);
    end loop;
    fxrco := fxnrc/fxdrc;
    rcogd := maxng/mindg;
  end Jacobian_with_Inverse_Condition;

end VarbPrec_Gradient_Evaluations;
