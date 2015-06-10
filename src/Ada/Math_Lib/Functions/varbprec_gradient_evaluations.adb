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
               gxnrc,gxdrc,gxrco : out double_float ) is

    use Standard_Complex_Numbers;
    use Standard_Gradient_Evaluations;

    numcnd : Standard_Floating_Vectors.Vector(0..x'last);
    val,rco : double_float;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    fxdrc := numcnd(0);
    fxrco := fxnrc/fxdrc;
    gxnrc := AbsVal(ydx(1));
    gxdrc := numcnd(1);
    gxrco := gxnrc/gxdrc;
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      rco := val/numcnd(k);
      if rco < gxrco then
        gxnrc := val;
        gxdrc := numcnd(k);
        gxrco := rco;
      end if;
    end loop;
  end Gradient_with_Inverse_Condition;

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.VecVec;
               ydx : out DoblDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_double;
               gxnrc,gxdrc,gxrco : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Gradient_Evaluations;

    numcnd : Double_Double_Vectors.Vector(0..x'last);
    val,rco : double_double;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    fxdrc := numcnd(0);
    fxrco := fxnrc/fxdrc;
    gxnrc := AbsVal(ydx(1));
    gxdrc := numcnd(1);
    gxrco := gxnrc/gxdrc;
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      rco := val/numcnd(k);
      if rco < gxrco then
        gxnrc := val;
        gxdrc := numcnd(k);
        gxrco := rco;
      end if;
    end loop;
  end Gradient_with_Inverse_Condition;

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.VecVec;
               ydx : out QuadDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out quad_double;
               gxnrc,gxdrc,gxrco : out quad_double ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Gradient_Evaluations;

    numcnd : Quad_Double_Vectors.Vector(0..x'last);
    val,rco : quad_double;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    fxdrc := numcnd(0);
    fxrco := fxnrc/fxdrc;
    gxnrc := AbsVal(ydx(1));
    gxdrc := numcnd(1);
    gxrco := gxnrc/gxdrc;
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      rco := val/numcnd(k);
      if rco < gxrco then
        gxnrc := val;
        gxdrc := numcnd(k);
        gxrco := rco;
      end if;
    end loop;
  end Gradient_with_Inverse_Condition;

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               gxnrc,gxdrc,gxrco : out Floating_Number ) is

    use Multprec_Complex_Numbers;
    use Multprec_Gradient_Evaluations;

    numcnd : Multprec_Floating_Vectors.Vector(0..x'last);
    val,rco : Floating_Number;

  begin
    Conditioned_Gradient_of_Polynomial(f,b,c,x,wrk,ydx,numcnd);
    fxnrc := AbsVal(ydx(0));
    Copy(numcnd(0),fxdrc);
    fxrco := fxnrc/fxdrc;
    gxnrc := AbsVal(ydx(1));
    Copy(numcnd(1),gxdrc);
    gxrco := gxnrc/gxdrc;
    for k in 2..ydx'last loop
      val := AbsVal(ydx(k));
      rco := val/numcnd(k);
      if rco < gxrco then
        Copy(val,gxnrc);
        Copy(numcnd(k),gxdrc);
        Copy(rco,gxrco);
      end if;
      Clear(val); Clear(rco);
    end loop;
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
               gxnrc,gxdrc,gxrco : out double_float ) is

    kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco : double_float;

  begin
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco);
      if kfxrco < fxrco then
        fxrco := kfxrco;
        fxnrc := kfxnrc;
        fxdrc := kfxdrc;
      end if;
      if kgxrco < gxrco then
        gxrco := kgxrco;
        gxnrc := kgxnrc;
        gxdrc := kgxdrc;
      end if;
    end loop;
  end Jacobian_with_Inverse_Condition;

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in DoblDobl_Complex_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.Array_of_VecVecs;
               ydx : in DoblDobl_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out double_double;
               gxnrc,gxdrc,gxrco : out double_double ) is

    kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco : double_double;

  begin
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco);
      if kfxrco < fxrco then
        fxrco := kfxrco;
        fxnrc := kfxnrc;
        fxdrc := kfxdrc;
      end if;
      if kgxrco < gxrco then
        gxrco := kgxrco;
        gxnrc := kgxnrc;
        gxdrc := kgxdrc;
      end if;
    end loop;
  end Jacobian_with_Inverse_Condition;

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in QuadDobl_Complex_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.Array_of_VecVecs;
               ydx : in QuadDobl_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out quad_double;
               gxnrc,gxdrc,gxrco : out quad_double ) is

    kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco : quad_double;

  begin
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco);
      if kfxrco < fxrco then
        fxrco := kfxrco;
        fxnrc := kfxnrc;
        fxdrc := kfxdrc;
      end if;
      if kgxrco < gxrco then
        gxrco := kgxrco;
        gxnrc := kgxnrc;
        gxdrc := kgxdrc;
      end if;
    end loop;
  end Jacobian_with_Inverse_Condition;

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in Multprec_Complex_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
               ydx : in Multprec_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               gxnrc,gxdrc,gxrco : out Floating_Number ) is

    kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco : Floating_Number;

  begin
    Clear(fxnrc); Clear(fxdrc); Clear(fxrco);
    Clear(gxnrc); Clear(gxdrc); Clear(gxrco);
    Gradient_with_Inverse_Condition
      (f(1).all,b(1).all,c(1).all,x,wrk(1).all,ydx(1).all,
       fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco);
    for k in 2..b'last loop
      Gradient_with_Inverse_Condition
        (f(k).all,b(k).all,c(k).all,x,wrk(k).all,ydx(k).all,
         kfxnrc,kfxdrc,kfxrco,kgxnrc,kgxdrc,kgxrco);
      if kfxrco < fxrco then
        Copy(kfxrco,fxrco);
        Copy(kfxnrc,fxnrc);
        Copy(kfxdrc,fxdrc);
      end if;
      if kgxrco < gxrco then
        Copy(kgxrco,gxrco);
        Copy(kgxnrc,gxnrc);
        Copy(kgxdrc,gxdrc);
      end if;
      Clear(kfxnrc); Clear(kfxdrc); Clear(kfxrco);
      Clear(kgxnrc); Clear(kgxdrc); Clear(kgxrco);
    end loop;
  end Jacobian_with_Inverse_Condition;

end VarbPrec_Gradient_Evaluations;
