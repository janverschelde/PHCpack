with unchecked_deallocation;
with Standard_Mathematical_Functions;
with Standard_Complex_Singular_Values;
with Exponent_Indices;
with Standard_Hessian_Updaters;

package body Standard_Coefficient_Circuits is

-- CONSTRUCTORS :

  function Allocate ( nbr,dim : integer32 ) return Circuit is

    res : Circuit(nbr);
    rforward : constant Standard_Floating_Vectors.Vector(1..dim-1)
             := (1..dim-1 => 0.0);
    iforward : constant Standard_Floating_Vectors.Vector(1..dim-1)
             := (1..dim-1 => 0.0);
    rbackward : constant Standard_Floating_Vectors.Vector(1..dim-2)
              := (1..dim-2 => 0.0);
    ibackward : constant Standard_Floating_Vectors.Vector(1..dim-2)
              := (1..dim-2 => 0.0);
    rcross : constant Standard_Floating_Vectors.Vector(1..dim-2)
           := (1..dim-2 => 0.0);
    icross : constant Standard_Floating_Vectors.Vector(1..dim-2)
           := (1..dim-2 => 0.0);

  begin
    res.dim := dim;
    res.pdg := -1;
    res.rfwd := new Standard_Floating_Vectors.Vector'(rforward);
    res.ifwd := new Standard_Floating_Vectors.Vector'(iforward);
    res.rbck := new Standard_Floating_Vectors.Vector'(rbackward);
    res.ibck := new Standard_Floating_Vectors.Vector'(ibackward);
    res.rcrs := new Standard_Floating_Vectors.Vector'(rcross);
    res.icrs := new Standard_Floating_Vectors.Vector'(icross);
    return res;
  end Allocate;

  function Exponent_Maxima
             ( c : Circuits; dim : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(c(c'first).xps);

  begin
    for k in c'first+1..c'last loop
      declare
        mxe : constant Standard_Integer_Vectors.Vector(1..dim)
            := Exponent_Indices.Maxima(c(k).xps);
      begin
        for i in mxe'range loop
          if mxe(i) > res(i)
           then res(i) := mxe(i);
          end if;
        end loop;
      end;
    end loop;
    return res;
  end Exponent_Maxima;

  function Create ( c : Circuits; dim : integer32 ) return System is

    res : System(c'last,dim);

  begin
    res.crc := c;
    res.mxe := Exponent_Maxima(c,dim);
    res.rpwt := Allocate(res.mxe);
    res.ipwt := Allocate(res.mxe);
    res.ryd := new Standard_Floating_Vectors.Vector'(0..dim => 0.0);
    res.iyd := new Standard_Floating_Vectors.Vector'(0..dim => 0.0);
    Allocate_Jacobian_Space(res.neq,res.dim,res.jrc,res.jic);
    return res;
  end Create;

  procedure Allocate_Jacobian_Space
              ( neq,dim : in integer32;
                jrc,jic : out Standard_Floating_VecVecs.VecVec ) is
  begin
    for k in 1..dim loop
      declare
        row : constant Standard_Floating_Vectors.Vector(1..dim)
            := (1..neq => 0.0);
      begin
        jrc(k) := new Standard_Floating_Vectors.Vector'(row);
        jic(k) := new Standard_Floating_Vectors.Vector'(row);
      end;
    end loop;
  end Allocate_Jacobian_Space;

  procedure Allocate_Jacobian_Space
              ( neq,dim : in integer32;
                jrc,jic : out Standard_Floating_VecVecs.Link_to_VecVec ) is

    jrcols : Standard_Floating_VecVecs.VecVec(1..dim);
    jicols : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    for k in 1..dim loop
      declare
        row : constant Standard_Floating_Vectors.Vector(1..dim)
            := (1..neq => 0.0);
      begin
        jrcols(k) := new Standard_Floating_Vectors.Vector'(row);
        jicols(k) := new Standard_Floating_Vectors.Vector'(row);
      end;
    end loop;
    jrc := new Standard_Floating_VecVecs.VecVec'(jrcols);
    jic := new Standard_Floating_VecVecs.VecVec'(jicols);
  end Allocate_Jacobian_Space;

  procedure Allocate_Hessian_Space
              ( dim : in integer32;
                hrp,hip : out Standard_Floating_VecVecs.VecVec ) is
  begin
    for k in 1..dim loop
      declare
        row : constant Standard_Floating_Vectors.Vector(1..dim)
            := (1..dim => 0.0);
      begin
        hrp(k) := new Standard_Floating_Vectors.Vector'(row);
        hip(k) := new Standard_Floating_Vectors.Vector'(row);
      end;
    end loop;
  end Allocate_Hessian_Space;

  procedure Allocate_Hessian_Space
              ( dim : in integer32;
                hrp,hip : out Standard_Floating_VecVecs.Link_to_VecVec ) is

    hrprows : Standard_Floating_VecVecs.VecVec(1..dim);
    hiprows : Standard_Floating_VecVecs.VecVec(1..dim);

  begin
    for k in 1..dim loop
      declare
        row : constant Standard_Floating_Vectors.Vector(1..dim)
            := (1..dim => 0.0);
      begin
        hrprows(k) := new Standard_Floating_Vectors.Vector'(row);
        hiprows(k) := new Standard_Floating_Vectors.Vector'(row);
      end;
    end loop;
    hrp := new Standard_Floating_VecVecs.VecVec'(hrprows);
    hip := new Standard_Floating_VecVecs.VecVec'(hiprows);
  end Allocate_Hessian_Space;
 
  procedure Allocate_Jacobian_Space ( s : in out System ) is
  begin
    Allocate_Jacobian_Space(s.neq,s.dim,s.jrc,s.jic);
  end Allocate_Jacobian_Space;

  procedure Allocate_Jacobian_Space ( s : in Link_to_System ) is
  begin
    if s /= null
     then Allocate_Jacobian_Space(s.all);
    end if;
  end Allocate_Jacobian_Space;
 
  procedure Allocate_Hessian_Space ( s : in out System ) is
  begin
    Allocate_Hessian_Space(s.dim,s.hrp,s.hip);
  end Allocate_Hessian_Space;

  procedure Allocate_Hessian_Space ( s : in Link_to_System ) is
  begin
    if s /= null
     then Allocate_Hessian_Space(s.all);
    end if;
  end Allocate_Hessian_Space;

  function Merge ( hrp,hip : Standard_Floating_VecVecs.VecVec )
                 return Standard_Complex_Matrices.Matrix is

    dim : constant integer32 := hrp'last;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in 1..dim loop
      hrprow := hrp(i); hiprow := hip(i);
      for j in 1..dim loop
        res(i,j) := Standard_Complex_Numbers.Create(hrprow(j),hiprow(j));
      end loop;
    end loop;
    return res;
  end Merge;

  function Copy ( c : Circuit ) return Circuit is

    res : Circuit(c.nbr) := Allocate(c.nbr,c.dim);

  begin
    res.pdg := c.pdg;
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    Standard_Floating_Vectors.Copy(c.rcf,res.rcf);
    Standard_Floating_Vectors.Copy(c.icf,res.icf);
    res.rcst := c.rcst;
    res.icst := c.icst;
    return res;
  end Copy;

  function Copy ( c : Link_to_Circuit ) return Link_to_Circuit is

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        cpc : constant Circuit(c.nbr) := Copy(c.all);
      begin
        res := new Circuit'(cpc);
      end;
    end if;
    return res;
  end Copy;

  function Copy ( c : Circuits ) return Circuits is

    res : Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := Copy(c(k));
    end loop;
    return res;
  end Copy;

  function Copy ( s : System ) return System is

    res : System(s.neq,s.dim);
    crc : constant Circuits(s.crc'range) := Copy(s.crc);

  begin
    res := Create(crc,s.dim);
    Allocate_Jacobian_Space(res);
    Allocate_Hessian_Space(res);
    return res;
  end Copy;

  function Copy ( s : Link_to_System ) return Link_to_System is

    res : Link_to_System;

  begin
    if s /= null then
      declare
        cps : constant System(s.neq,s.dim) := Copy(s.all);
      begin
        res := new System'(cps);
      end;
    end if;
    return res;
  end Copy;

-- RADIUS COEFFICIENTS :

  procedure AbsVal ( c : in out Circuit ) is

    use Standard_Mathematical_Functions;

  begin
    c.rcst := SQRT(c.rcst**2 + c.icst**2);
    c.icst := 0.0;
    for k in 1..c.nbr loop
      c.rcf(k) := SQRT(c.rcf(k)**2 + c.icf(k)**2);
      c.icf(k) := 0.0;
    end loop;
  end AbsVal;

  procedure AbsVal ( c : in Link_to_Circuit ) is
  begin
    if c /= null
     then AbsVal(c.all);
    end if;
  end AbsVal;

  procedure AbsVal ( c : in Circuits ) is
  begin
    for k in c'range loop
      AbsVal(c(k));
    end loop;
  end AbsVal;

  procedure AbsVal ( s : in System ) is
  begin
    AbsVal(s.crc);
  end AbsVal;

  procedure AbsVal ( s : in Link_to_System ) is
  begin
    if s /= null
     then AbsVal(s.all);
    end if;
  end AbsVal;

-- EVALUATION OF CIRCUITS :

  function Eval ( c : in Circuit;
                  xr : in Standard_Floating_Vectors.Link_to_Vector;
                  rpwt : in Standard_Floating_VecVecs.VecVec )
                return double_float is

    res : double_float := c.rcst;
    rkcff,prd : double_float;
    idx : integer32;

    use Standard_Integer_Vectors;
 
  begin
    for k in c.xps'range loop
      declare
        xpk : constant Standard_Integer_Vectors.Link_to_Vector := c.xps(k);
        idk : constant Standard_Integer_Vectors.Vector := c.idx(k).all;
        fck : constant Standard_Integer_Vectors.Link_to_Vector := c.fac(k);
      begin
        rkcff := c.rcf(k);
        if fck = null then        -- no common factor
          if idk'last = 1 then    -- y := y + kcff*x(idx1);
            idx := idk(1);
            res := res + rkcff*xr(idx);
          else -- multiply coefficient with forward product
            prd := rkcff;
            for i in idk'range loop
              idx := idk(i);
              prd := prd*xr(idx);
            end loop;
            res := res + prd;
          end if;
        else -- compute the common factor first
          Multiply_Factor(xpk,fck,xr,rkcff,rpwt,prd);
          for i in idk'range loop -- multiply with forward product
            idx := idk(i);
            prd := prd*xr(idx);
          end loop;
          res := res + prd;
        end if;
      end;
    end loop;
    return res;
  end Eval;

  function Eval ( c : in Circuit;
                  xr : in Standard_Floating_Vectors.Link_to_Vector;
                  xi : in Standard_Floating_Vectors.Link_to_Vector;
                  rpwt : in Standard_Floating_VecVecs.VecVec;
                  ipwt : in Standard_Floating_VecVecs.VecVec ) 
                return Standard_Complex_Numbers.Complex_Number is

    res : Standard_Complex_Numbers.Complex_Number;
    rpres : double_float := c.rcst; -- real part of result
    ipres : double_float := c.icst; -- imaginary part of result
    rkcff,ikcff,zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;

    use Standard_Integer_Vectors;
 
  begin
    for k in c.xps'range loop
      declare
        xpk : constant Standard_Integer_Vectors.Link_to_Vector := c.xps(k);
        idk : constant Standard_Integer_Vectors.Vector := c.idx(k).all;
        fck : constant Standard_Integer_Vectors.Link_to_Vector := c.fac(k);
      begin
        rkcff := c.rcf(k); ikcff := c.icf(k);
        if fck = null then        -- no common factor
          if idk'last = 1 then    -- y := y + kcff*x(idx1);
            idx := idk(1);
            pr := xr(idx); pi := xi(idx);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            rpres := rpres + zr; ipres := ipres + zi;
          else -- multiply coefficient with forward product
            zr := rkcff; zi := ikcff;
            for i in idk'range loop
              pr := zr; pi := zi;
              idx := idk(i); qr := xr(idx); qi := xi(idx);
              zr := pr*qr - pi*qi; zi := pr*qi + pi*qr;
            end loop;
            rpres := rpres + zr; ipres := ipres + zi;
          end if;
        else -- compute the common factor first
          Multiply_Factor(xpk,fck,xr,xi,rkcff,ikcff,rpwt,ipwt,zr,zi);
          for i in idk'range loop -- multiply with forward product
            pr := zr; pi := zi;
            idx := idk(i); qr := xr(idx); qi := xi(idx);
            zr := pr*qr - pi*qi; zi := pr*qi + pi*qr;
          end loop;
          rpres := rpres + zr; ipres := ipres + zi;
        end if;
      end;
    end loop;
    res := Standard_Complex_Numbers.Create(rpres,ipres);
    return res;
  end Eval;

  procedure Eval ( s : in out System;
                   xr : in Standard_Floating_Vectors.Link_to_Vector;
                   xi : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    for i in s.crc'range loop
      s.fx(i) := Eval(s.crc(i).all,xr,xi,s.rpwt,s.ipwt);
    end loop;
  end Eval;

  procedure Eval ( s : in Link_to_System;
                   xr : in Standard_Floating_Vectors.Link_to_Vector;
                   xi : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    for i in s.crc'range loop
      s.fx(i) := Eval(s.crc(i).all,xr,xi,s.rpwt,s.ipwt);
    end loop;
  end Eval;

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUITS :

  procedure EvalDiff
              ( s : in out System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    EvalDiff(s.crc,xr,xi,s.ryd,s.iyd,s.rpwt,s.ipwt,s.jrc.all,s.jic.all,
             s.fx,s.jm);
  end EvalDiff;

  procedure EvalDiff
              ( s : in Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    EvalDiff(s.crc,xr,xi,s.ryd,s.iyd,s.rpwt,s.ipwt,s.jrc.all,s.jic.all,
             s.fx,s.jm);
  end EvalDiff;

  procedure EvalDiff2
              ( s : in out System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat ) is
  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    EvalDiff2(s.crc,xr,xi,s.ryd,s.iyd,s.rpwt,s.ipwt,s.jrc.all,s.jic.all,
              s.hrp.all,s.hip.all,s.fx,s.jm,vh);
  end EvalDiff2;

  procedure EvalDiff2
              ( s : in Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat ) is
  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    EvalDiff2(s.crc,xr,xi,s.ryd,s.iyd,s.rpwt,s.ipwt,s.jrc.all,s.jic.all,
              s.hrp.all,s.hip.all,s.fx,s.jm,vh);
  end EvalDiff2;

  procedure EvalDiff
              ( c : in Circuits;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                jrc : in Standard_Floating_VecVecs.VecVec;
                jic : in Standard_Floating_VecVecs.VecVec;
                fx : out Standard_Complex_Vectors.Vector;
                jm : out Standard_Complex_Matrices.Matrix ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in c'range loop
      Speel(c(i).all,xr,xi,ryd,iyd,rpwt,ipwt);
      fx(i) := Standard_Complex_Numbers.Create(ryd(0),iyd(0));
      for j in jm'range(2) loop
        jm(i,j) := Standard_Complex_Numbers.Create(ryd(j),iyd(j));
        rlnk := jrc(j);    ilnk := jic(j);
        rlnk(i) := ryd(j); ilnk(i) := iyd(j);
        ryd(j) := 0.0;     iyd(j) := 0.0;
      end loop;
    end loop;
  end EvalDiff;

  procedure EvalDiff2
              ( c : in Circuits;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                jrc : in Standard_Floating_VecVecs.VecVec;
                jic : in Standard_Floating_VecVecs.VecVec;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec;
                fx : out Standard_Complex_Vectors.Vector;
                jm : out Standard_Complex_Matrices.Matrix;
                vh : in Standard_Complex_VecMats.VecMat ) is

    mat : Standard_Complex_Matrices.Link_to_Matrix;
    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in c'range loop
      mat := vh(i);
      Speel(c(i).all,xr,xi,ryd,iyd,rpwt,ipwt,hrp,hip);
      fx(i) := Standard_Complex_Numbers.Create(ryd(0),iyd(0));
      for j in jm'range(2) loop
        jm(i,j) := Standard_Complex_Numbers.Create(ryd(j),iyd(j));
        rlnk := jrc(j);    ilnk := jic(j);
        rlnk(i) := ryd(j); ilnk(i) := iyd(j);
        ryd(j) := 0.0;     iyd(j) := 0.0;
      end loop;
      for j in hrp'range loop
        hrprow := hrp(j); hiprow := hip(j);
        for k in hrprow'range loop
          mat(j,k) := Standard_Complex_Numbers.Create(hrprow(k),hiprow(k));
        end loop;
      end loop;
    end loop;
  end EvalDiff2;

-- SINGULAR VALUE DECOMPOSITIONS :

  procedure Singular_Values
              ( s : in out System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                svls : in Standard_Complex_VecVecs.VecVec ) is

    info : integer32;

  begin
    Power_Table(s.mxe,xr,xi,s.rpwt,s.ipwt);
    EvalDiff2(s.crc,xr,xi,s.ryd,s.iyd,s.rpwt,s.ipwt,s.jrc.all,s.jic.all,
              s.hrp.all,s.hip.all,s.fx,s.jm,vh);
    Standard_Complex_Singular_Values.SVD
      (s.jm,s.dim,s.dim,svls(0).all,e,U,V,0,info);
    for k in vh'range loop
      Standard_Complex_Singular_Values.SVD
        (vh(k).all,s.dim,s.dim,svls(k).all,e,U,V,0,info);
    end loop;
  end Singular_Values;

  procedure Singular_Values
              ( s : in Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                svls : in Standard_Complex_VecVecs.VecVec ) is
  begin
    if s /= null
     then Singular_Values(s.all,xr,xi,vh,U,V,e,svls);
    end if;
  end Singular_Values;

  procedure Singular_Values
              ( c : in Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Vectors.Vector ) is

    info : integer32;

  begin
    if c.pdg <= 1 then
      s := (s'range => Standard_Complex_Numbers.Create(0.0));
    else
      Speel(c,xr,xi,ryd,iyd,rpwt,ipwt,hrp,hip);
      A := Merge(hrp,hip);
      Standard_Complex_Singular_Values.SVD(A,c.dim,c.dim,s,e,U,V,0,info);
    end if;
  end Singular_Values;

  procedure Singular_Values
              ( c : in Link_to_Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Vectors.Vector ) is
  begin
    if c /= null then
      Singular_Values(c.all,xr,xi,ryd,iyd,rpwt,ipwt,hrp,hip,A,U,V,e,s);
    end if;
  end Singular_Values;

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUITS :

  procedure Indexed_Speel 
              ( c : in Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec ) is
  begin
    Indexed_Speel(c.xps,c.rcf,c.icf,c.rcst,c.icst,xr,xi,ryd,iyd,
                  c.rfwd,c.ifwd,c.rbck,c.ibck,c.rcrs,c.icrs,hrp,hip);
  end Indexed_Speel;

  procedure Indexed_Speel
              ( c : in Circuit;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    Indexed_Speel(c.xps,c.rcf,c.icf,c.rcst,c.icst,xr,xi,ryd,iyd,
                  c.rfwd,c.ifwd,c.rbck,c.ibck,c.rcrs,c.icrs);
  end Indexed_Speel;

  procedure Speel ( c : in Circuit;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec ) is
  begin
    Speel(c.xps,c.idx,c.fac,c.rcf,c.icf,c.rcst,c.icst,xr,xi,ryd,iyd,
          c.rfwd,c.ifwd,c.rbck,c.ibck,c.rcrs,c.icrs,rpwt,ipwt);
  end Speel;

  procedure Speel ( c : in Circuit;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec; 
                    hrp : in Standard_Floating_VecVecs.VecVec;
                    hip : in Standard_Floating_VecVecs.VecVec ) is
  begin
    Speel(c.xps,c.idx,c.fac,c.rcf,c.icf,c.rcst,c.icst,xr,xi,ryd,iyd,
          c.rfwd,c.ifwd,c.rbck,c.ibck,c.rcrs,c.icrs,rpwt,ipwt,hrp,hip);
    null;
  end Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_Vectors.Vector;
                rcff,icff : in double_float;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                rbck : in Standard_Floating_Vectors.Link_to_Vector;
                ibck : in Standard_Floating_Vectors.Link_to_Vector;
                rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                icrs : in Standard_Floating_Vectors.Link_to_Vector;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec ) is

    sz : constant integer32 := idx'last;
    idx1,idx2,idx3,idx4,idx5,idx6 : integer32;
    zr,zi,pr,pi,racc,iacc : double_float;
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Fused_Forward_Backward_Cross
      (idx,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
   -- cff := c.cff(k); yd(0) := yd(0) + cff*c.fwd(sz-1); 
    idx1 := sz-1; pr := rfwd(idx1); pi := ifwd(idx1);
    zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
    ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
    if sz = 2 then
      idx1 := idx(1); idx2 := idx(2);
     -- yd(idx2) := yd(idx2) + cff*x(idx1);
      pr := xr(idx1); pi := xi(idx1);
      zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
      ryd(idx2) := ryd(idx2) + zr; iyd(idx2) := iyd(idx2) + zi;
     -- yd(idx1) := yd(idx1) + cff*x(idx2);
      pr := xr(idx2); pi := xi(idx2);
      zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
      ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
     -- h(idx1,idx2) := h(idx1,idx2) + cff;
      hrprow := hrp(idx1); hiprow := hip(idx1);
      hrprow(idx2) := hrprow(idx2) + rcff;
      hiprow(idx2) := hiprow(idx2) + icff;
    else -- sz'last > 2
      idx1 := idx(1); 
     -- yd(idx1) := yd(idx1) + cff*c.bck(sz-2);
      idx2 := sz-2; pr := rbck(idx2); pi := ibck(idx2);
      zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
      ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
      for j in idx'first+1..sz-1 loop
        idx1 := idx(j);
       -- yd(idx1) := yd(idx1) + cff*c.crs(j-1);
        idx2 := j-1; pr := rcrs(idx2); pi := icrs(idx2);
        zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
        ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
      end loop;
      idx1 := idx(sz);
     -- yd(idx1) := yd(idx1) + cff*c.fwd(sz-2);
      idx2 := sz-2; pr := rfwd(idx2); pi := ifwd(idx2);
      zr := pr*rcff - pi*icff; zi := pr*icff + pi*rcff;
      ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
      if sz = 3 then
        idx1 := idx(1); idx2 := idx(2); idx3 := idx(3);
       -- h(idx1,idx2) := h(idx1,idx2) + cff*x(idx3);
        hrprow := hrp(idx1); hiprow := hip(idx1);
        pr := xr(idx3); pi := xi(idx3);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx2) := hrprow(idx2) + zr;
        hiprow(idx2) := hiprow(idx2) + zi;
       -- h(idx1,idx3) := h(idx1,idx3) + cff*x(idx2);
        pr := xr(idx2); pi := xi(idx2);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx3) := hrprow(idx3) + zr;
        hiprow(idx3) := hiprow(idx3) + zi;
       -- h(idx2,idx3) := h(idx2,idx3) + cff*x(idx1);
        hrprow := hrp(idx2); hiprow := hip(idx2);
        pr := xr(idx1); pi := xi(idx1);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx3) := hrprow(idx3) + zr;
        hiprow(idx3) := hiprow(idx3) + zi;
      else -- sz > 3
        idx5 := idx(sz-1); idx6 := idx(sz); idx3 := sz-3;
       -- last element is copy of fwd(sz-3), multiplied with cff
       -- h(idx5,idx6) := h(idx5,idx6) + cff*fwd(idx3);
        hrprow := hrp(idx5); hiprow := hip(idx5);
        pr := rfwd(idx3); pi := ifwd(idx3);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx6) := hrprow(idx6) + zr;
        hiprow(idx6) := hiprow(idx6) + zi;
       -- first element is copy of bck(sz-3), multiplied with cff
        idx1 := idx(1); idx2 := idx(2);
       -- h(idx1,idx2) := h(idx1,idx2) + cff*bck(idx3);
        hrprow := hrp(idx1); hiprow := hip(idx1);
        pr := rbck(idx3); pi := ibck(idx3);
        zr := rcff*pr - icff*pi; zi := rcff*pi + icff*pr;
        hrprow(idx2) := hrprow(idx2) + zr;
        hiprow(idx2) := hiprow(idx2) + zi;
        if sz = 4 then -- special case for all rows
          idx3 := idx(3); idx4 := idx(4);
         -- acc := cff*x(idx2);
          pr := xr(idx2); pi := xi(idx2);
          racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
         -- h(idx1,idx3) := h(idx1,idx3) + acc*x(idx6);
          hrprow := hrp(idx1); hiprow := hip(idx1);
          pr := xr(idx6); pi := xi(idx6);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx3) := hrprow(idx3) + zr;
          hiprow(idx3) := hiprow(idx3) + zi;
         -- h(idx1,idx4) := h(idx1,idx4) + acc*x(idx5);
          pr := xr(idx5); pi := xi(idx5);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx4) := hrprow(idx4) + zr;
          hiprow(idx4) := hiprow(idx4) + zi;
         -- acc := cff*x(idx1);
          pr := xr(idx1); pi := xi(idx1);
          racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
         -- h(idx2,idx3) := h(idx2,idx3) + acc*x(idx6);
          hrprow := hrp(idx2); hiprow := hip(idx2);
          pr := xr(idx6); pi := xi(idx6);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx3) := hrprow(idx3) + zr;
          hiprow(idx3) := hiprow(idx3) + zi;
         -- h(idx2,idx4) := h(idx2,idx4) + acc*x(idx5);
          pr := xr(idx5); pi := xi(idx5);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx4) := hrprow(idx4) + zr;
          hiprow(idx4) := hiprow(idx4) + zi;
        else -- sz > 4
         -- first row is special, starts with x(idx(2)) after diagonal
          idx3 := idx(3); idx4 := sz-4;
         -- acc := cff*x(idx2);
          pr := xr(idx2); pi := xi(idx2);
          racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
         -- h(idx1,idx3) := h(idx1,idx3) + acc*bck(idx4);
          hrprow := hrp(idx1); hiprow := hip(idx1);
          pr := rbck(idx4); pi := ibck(idx4);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx3) := hrprow(idx3) + zr;
          hiprow(idx3) := hiprow(idx3) + zi;
          for k in 4..sz-2 loop
            idx4 := idx(k-1); idx5 := idx(k); idx6 := sz-k-1;
           -- acc := acc*x(idx4);
            pr := xr(idx4); pi := xi(idx4);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            racc := zr; iacc := zi;
           -- h(idx1,idx5) := h(idx1,idx5) + acc*bck(idx6);
            pr := rbck(idx6); pi := ibck(idx6);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            hrprow(idx5) := hrprow(idx5) + zr;
            hiprow(idx5) := hiprow(idx5) + zi;
          end loop;
          idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
         -- acc := acc*x(idx4);
          pr := xr(idx4); pi := xi(idx4);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          racc := zr; iacc := zi;
         -- h(idx1,idx5) := h(idx1,idx5) + acc*x(idx6);
          pr := xr(idx6); pi := xi(idx6);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx5) := hrprow(idx5) + zr;
          hiprow(idx5) := hiprow(idx5) + zi;
         -- h(idx1,idx6) := h(idx1,idx6) + acc*x(idx5);
          pr := xr(idx5); pi := xi(idx5);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx6) := hrprow(idx6) + zr;
          hiprow(idx6) := hiprow(idx6) + zi;
         -- second row is special, starts with x(idx(1)) after diagonal
         -- acc := cff*x(idx1);
          pr := xr(idx1); pi := xi(idx1);
          racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
          idx4 := sz-4;
         -- h(idx2,idx3) := h(idx2,idx3) + acc*bck(idx4);
          pr := rbck(idx4); pi := ibck(idx4);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow := hrp(idx2); hiprow := hip(idx2);
          hrprow(idx3) := hrprow(idx3) + zr;
          hiprow(idx3) := hiprow(idx3) + zi;
          for k in 4..sz-2 loop
            idx4 := idx(k-1); idx5 := idx(k); idx6 := sz-k-1;
           -- acc := acc*x(idx4);
            pr := xr(idx4); pi := xi(idx4);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            racc := zr; iacc := zi;
           -- h(idx2,idx5) := h(idx2,idx5) + acc*bck(idx6);
            pr := rbck(idx6); pi := ibck(idx6);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            hrprow(idx5) := hrprow(idx5) + zr;
            hiprow(idx5) := hiprow(idx5) + zi;
          end loop;
          idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
         -- acc := acc*x(idx4);
          pr := xr(idx4); pi := xi(idx4);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          racc := zr; iacc := zi;
         -- h(idx2,idx5) := h(idx2,idx5) + acc*x(idx6);
          pr := xr(idx6); pi := xi(idx6);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx5) := hrprow(idx5) + zr;
          hiprow(idx5) := hiprow(idx5) + zi;
         -- h(idx2,idx6) := h(idx2,idx6) + acc*x(idx5);
          pr := xr(idx5); pi := xi(idx5);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx6) := hrprow(idx6) + zr;
          hiprow(idx6) := hiprow(idx6) + zi;
         -- the row with index sz-2 has a general formula
          idx3 := sz-4;
         -- acc := cff*fwd(idx3);
          pr := rfwd(idx3); pi := ifwd(idx3);
          racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
         -- h(idx4,idx5) := h(idx4,idx5) + acc*x(idx6);
          pr := xr(idx6); pi := xi(idx6);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow := hrp(idx4); hiprow := hip(idx4);
          hrprow(idx5) := hrprow(idx5) + zr;
          hiprow(idx5) := hiprow(idx5) + zi;
         -- h(idx4,idx6) := h(idx4,idx6) + acc*x(idx5);
          pr := xr(idx5); pi := xi(idx5);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          hrprow(idx6) := hrprow(idx6) + zr;
          hiprow(idx6) := hiprow(idx6) + zi;
          for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
            idx1 := idx(rw); idx2 := idx(rw+1);
           -- acc := cff*fwd(rw-2);
            pr := rfwd(rw-2); pi := ifwd(rw-2);
            racc := rcff*pr - icff*pi; iacc := rcff*pi + icff*pr;
           -- h(idx1,idx2) := h(idx1,idx2) + acc*bck(sz-rw-2);
            hrprow := hrp(idx1); hiprow := hip(idx1);
            pr := rbck(sz-rw-2); pi := ibck(sz-rw-2);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            hrprow(idx2) := hrprow(idx2) + zr;
            hiprow(idx2) := hiprow(idx2) + zi;
            for k in rw+2..sz-2 loop
              idx2 := idx(k-1); idx3 := idx(k); idx4 := sz-k-1;
             -- acc := acc*x(idx2);
              pr := xr(idx2); pi := xi(idx2);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              racc := zr; iacc := zi;
             -- h(idx1,idx3) := h(idx1,idx3) + acc*bck(idx4);
              pr := rbck(idx4); pi := ibck(idx4);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              hrprow(idx3) := hrprow(idx3) + zr;
              hiprow(idx3) := hiprow(idx3) + zi;
            end loop;
            idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
           -- acc := acc*x(idx4);
            pr := xr(idx4); pi := xi(idx4);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            racc := zr; iacc := zi;
           -- h(idx1,idx5) := h(idx1,idx5) + acc*x(idx6);
            hrprow := hrp(idx1); hiprow := hip(idx1);
            pr := xr(idx6); pi := xi(idx6);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            hrprow(idx5) := hrprow(idx5) + zr;
            hiprow(idx5) := hiprow(idx5) + zi;
           -- h(idx1,idx6) := h(idx1,idx6) + acc*x(idx5);
            pr := xr(idx5); pi := xi(idx5);
            zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
            hrprow(idx6) := hrprow(idx6) + zr;
            hiprow(idx6) := hiprow(idx6) + zi;
          end loop;
        end if;
      end if;
    end if;
  end Indexed_Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                rcf : in Standard_Floating_Vectors.Vector;
                icf : in Standard_Floating_Vectors.Vector;
                rcst,icst : in double_float;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                rbck : in Standard_Floating_Vectors.Link_to_Vector;
                ibck : in Standard_Floating_Vectors.Link_to_Vector;
                rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                icrs : in Standard_Floating_Vectors.Link_to_Vector;
                hrp : in Standard_Floating_VecVecs.VecVec;
                hip : in Standard_Floating_VecVecs.VecVec ) is

    dim : constant integer32 := xr'last;
    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1 : integer32;
    rkcff,ikcff : double_float;
    zr,zi,pr,pi : double_float;
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

    use Standard_Integer_Vectors;

  begin
    for i in 1..dim loop
      hrprow := hrp(i); hiprow := hip(i);
      for j in 1..dim loop
        hrprow(j) := 0.0; hiprow(j) := 0.0;
      end loop;
    end loop;
    ryd(0) := rcst; iyd(0) := icst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        rkcff := rcf(k); ikcff := icf(k);
        if idk'last = 1 then
          idx1 := idk(1);
         -- yd(0) := yd(0) + cff*x(idx1);
          pr := xr(idx1); pi := xi(idx1);
          zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
          ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
         -- yd(idx1) := yd(idx1) + cff;
          ryd(idx1) := ryd(idx1) + rkcff; iyd(idx1) := iyd(idx1) + ikcff;
        else
          Indexed_Speel(idk.all,rkcff,ikcff,xr,xi,ryd,iyd,
                        rfwd,ifwd,rbck,ibck,rcrs,icrs,hrp,hip);
        end if;
      end if;
    end loop;
    for i in 2..dim loop
      hrprow := hrp(i); hiprow := hip(i);
      for j in 1..(i-1) loop
        hrprow(j) := hrp(j)(i); hiprow(j) := hip(j)(i);
      end loop;
    end loop;
  end Indexed_Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                rcf : in Standard_Floating_Vectors.Vector;
                icf : in Standard_Floating_Vectors.Vector;
                rcst,icst : in double_float;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                ryd : in Standard_Floating_Vectors.Link_to_Vector;
                iyd : in Standard_Floating_Vectors.Link_to_Vector;
                rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                rbck : in Standard_Floating_Vectors.Link_to_Vector;
                ibck : in Standard_Floating_Vectors.Link_to_Vector;
                rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                icrs : in Standard_Floating_Vectors.Link_to_Vector ) is

    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    rkcff,ikcff : double_float;
    zr,zi,pr,pi : double_float;

    use Standard_Integer_Vectors;

  begin
    ryd(0) := rcst; iyd(0) := icst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          idx1 := idk(1); rkcff := rcf(k); ikcff := icf(k);
         -- yd(0) := yd(0) + cff*x(idx1);
          pr := xr(idx1); pi := xi(idx1);
          zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
          ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
         -- yd(idx1) := yd(idx1) + cff;
          ryd(idx1) := ryd(idx1) + rkcff; iyd(idx1) := iyd(idx1) + ikcff;
        else
          Fused_Forward_Backward_Cross
            (idk.all,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
         -- cff := c.cff(k); yd(0) := yd(0) + cff*c.fwd(idk'last-1); 
          rkcff := rcf(k); ikcff := icf(k);
          idx1 := idk'last-1; pr := rfwd(idx1); pi := ifwd(idx1);
          zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
          ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
          if idk'last = 2 then
            idx1 := idk(1); idx2 := idk(2);
           -- yd(idx2) := yd(idx2) + cff*x(idx1);
            pr := xr(idx1); pi := xi(idx1);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx2) := ryd(idx2) + zr; iyd(idx2) := iyd(idx2) + zi;
           -- yd(idx1) := yd(idx1) + cff*x(idx2);
            pr := xr(idx2); pi := xi(idx2);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
          else -- idk'last > 2
            idx1 := idk(1); 
           -- yd(idx1) := yd(idx1) + cff*c.bck(idk'last-2);
            idx2 := idk'last-2; pr := rbck(idx2); pi := ibck(idx2);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
            for j in idk'first+1..idk'last-1 loop
              idx1 := idk(j);
             -- yd(idx1) := yd(idx1) + cff*c.crs(j-1);
              idx2 := j-1; pr := rcrs(idx2); pi := icrs(idx2);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
            end loop;
            idx1 := idk(idk'last);
           -- yd(idx1) := yd(idx1) + cff*c.fwd(idk'last-2);
            idx2 := idk'last-2; pr := rfwd(idx2); pi := ifwd(idx2);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
          end if;
        end if;
      end if;
    end loop;
  end Indexed_Speel;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcf : in Standard_Floating_Vectors.Vector;
                    icf : in Standard_Floating_Vectors.Vector;
                    rcst,icst : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec ) is

    idk,fck,xpk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    rkcff,ikcff,racc,iacc,factor : double_float;
    zr,zi,pr,pi : double_float;

    use Standard_Integer_Vectors;

  begin
    ryd(0) := rcst; iyd(0) := icst;
    for k in idx'range loop
      idk := idx(k);
      rkcff := rcf(k); ikcff := icf(k);
      if idk = null then -- we have an extra constant coefficient
        ryd(0) := ryd(0) + rkcff;
        iyd(0) := iyd(0) + ikcff;
      else
        if idk'last < idk'first then -- deal with stray 0 exponent
          ryd(0) := ryd(0) + rkcff;
          iyd(0) := iyd(0) + ikcff;
        else
          idx1 := idk(1); 
          fck := fac(k);
          if fck = null then
            if idk'last = 1 then
             -- yd(0) := yd(0) + cff*x(idx1);
              pr := xr(idx1); pi := xi(idx1);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
             -- yd(idx1) := yd(idx1) + cff;
              ryd(idx1) := ryd(idx1) + rkcff; iyd(idx1) := iyd(idx1) + ikcff;
            else
              Fused_Forward_Backward_Cross
                (idk.all,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
             -- cff := c.cff(k); yd(0) := yd(0) + cff*c.fwd(idk'last-1); 
              rkcff := rcf(k); ikcff := icf(k);
              idx1 := idk'last-1; pr := rfwd(idx1); pi := ifwd(idx1);
              zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
              ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
              if idk'last = 2 then
                idx1 := idk(1); idx2 := idk(2);
               -- yd(idx2) := yd(idx2) + cff*x(idx1);
                pr := xr(idx1); pi := xi(idx1);
                zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
                ryd(idx2) := ryd(idx2) + zr; iyd(idx2) := iyd(idx2) + zi;
               -- yd(idx1) := yd(idx1) + cff*x(idx2);
                pr := xr(idx2); pi := xi(idx2);
                zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
                ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
              else -- idk'last > 2
                idx1 := idk(1); 
               -- yd(idx1) := yd(idx1) + cff*c.bck(idk'last-2);
                idx2 := idk'last-2; pr := rbck(idx2); pi := ibck(idx2);
                zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
                ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
                for j in idk'first+1..idk'last-1 loop
                  idx1 := idk(j);
                 -- yd(idx1) := yd(idx1) + cff*c.crs(j-1);
                  idx2 := j-1; pr := rcrs(idx2); pi := icrs(idx2);
                  zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
                  ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
                end loop;
                idx1 := idk(idk'last);
               -- yd(idx1) := yd(idx1) + cff*c.fwd(idk'last-2);
                idx2 := idk'last-2; pr := rfwd(idx2); pi := ifwd(idx2);
                zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
                ryd(idx1) := ryd(idx1) + zr; iyd(idx1) := iyd(idx1) + zi;
              end if;
            end if;
          else
            xpk := xps(k);
            if idk'last = 1 then
              Multiply_Factor(xpk,fck,xr,xi,rkcff,ikcff,rpwt,ipwt,racc,iacc);
             -- yd(0) := yd(0) + acc*x(idx1);
              pr := xr(idx1); pi := xi(idx1);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
             -- yd(idx1) := yd(idx1) + factor*acc;
              factor := Create(xpk(idx1));
              ryd(idx1) := ryd(idx1) + factor*racc;
              iyd(idx1) := iyd(idx1) + factor*iacc;
            else
              Fused_Forward_Backward_Cross
                (idk.all,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
              Multiply_Factor(xpk,fck,xr,xi,rkcff,ikcff,rpwt,ipwt,racc,iacc);
              pr := rfwd(idk'last-1); pi := ifwd(idk'last-1);
              zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
              ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
              if idk'last = 2 then
                idx2 := idk(2);
               -- wrk := acc*x(idx1); factor := Create(xpk(idx2));
               -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
                pr := xr(idx1); pi := xi(idx1);
                zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
                factor := Create(xpk(idx2));
                ryd(idx2) := ryd(idx2) + factor*zr;
                iyd(idx2) := iyd(idx2) + factor*zi;
               -- wrk := acc*x(idx2); factor := Create(xpk(idx1));
               -- wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
                pr := xr(idx2); pi := xi(idx2);
                zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
                factor := Create(xpk(idx1));
                ryd(idx1) := ryd(idx1) + factor*zr;
                iyd(idx1) := iyd(idx1) + factor*zi;
              else -- idk'last > 2
               -- wrk := acc*bck(idk'last-2);
                pr := rbck(idk'last-2); pi := ibck(idk'last-2);
                zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
                factor := Create(xpk(idx1));
               -- wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
                ryd(idx1) := ryd(idx1) + factor*zr;
                iyd(idx1) := iyd(idx1) + factor*zi;
                for j in idk'first+1..idk'last-1 loop
                 -- wrk := acc*crs(j-1);
                  pr := rcrs(j-1); pi := icrs(j-1);
                  zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
                  idx2 := idk(j); factor := Create(xpk(idx2));
                 -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
                  ryd(idx2) := ryd(idx2) + factor*zr;
                  iyd(idx2) := iyd(idx2) + factor*zi;
                end loop;
               -- wrk := acc*fwd(idk'last-2);
                idx2 := idk(idk'last);
                pr := rfwd(idk'last-2); pi := ifwd(idk'last-2);
                zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
                factor := Create(xpk(idx2));
               -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
                ryd(idx2) := ryd(idx2) + factor*zr;
                iyd(idx2) := iyd(idx2) + factor*zi;
              end if; -- if idx(k)'last is 2
            end if; -- if idx(k)'last is 1
          end if; -- if fac(k) is null
        end if; -- if idx(k)'last < idx(k)'first
      end if; -- if idx(k) is null
    end loop;
  end Speel;

  procedure Speel ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector;
                    rcff,icff : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec ) is

    idx1,idx2 : integer32;
    racc,iacc,factor : double_float;
    zr,zi,pr,pi : double_float;

  begin
    idx1 := idx(1);
    if idx'last = 1 then
      Multiply_Factor(xps,fac,xr,xi,rcff,icff,rpwt,ipwt,racc,iacc);
     -- yd(0) := yd(0) + acc*x(idx1);
      pr := xr(idx1); pi := xi(idx1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
     -- yd(idx1) := yd(idx1) + factor*acc;
      factor := Create(xps(idx1));
      ryd(idx1) := ryd(idx1) + factor*racc;
      iyd(idx1) := iyd(idx1) + factor*iacc;
    else
      Fused_Forward_Backward_Cross
        (idx,xr,xi,rfwd,ifwd,rbck,ibck,rcrs,icrs);
      Multiply_Factor(xps,fac,xr,xi,rcff,icff,rpwt,ipwt,racc,iacc);
      pr := rfwd(idx'last-1); pi := ifwd(idx'last-1);
      zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
      ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
      if idx'last = 2 then
        idx2 := idx(2);
       -- wrk := acc*x(idx1); factor := Create(xps(idx2));
       -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
        pr := xr(idx1); pi := xi(idx1);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        factor := Create(xps(idx2));
        ryd(idx2) := ryd(idx2) + factor*zr;
        iyd(idx2) := iyd(idx2) + factor*zi;
       -- wrk := acc*x(idx2); factor := Create(xps(idx1));
       -- wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
        pr := xr(idx2); pi := xi(idx2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        factor := Create(xps(idx1));
        ryd(idx1) := ryd(idx1) + factor*zr;
        iyd(idx1) := iyd(idx1) + factor*zi;
      else -- idx'last > 2
       -- wrk := acc*bck(idx'last-2);
        pr := rbck(idx'last-2); pi := ibck(idx'last-2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        factor := Create(xps(idx1));
       -- wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
        ryd(idx1) := ryd(idx1) + factor*zr;
        iyd(idx1) := iyd(idx1) + factor*zi;
        for j in idx'first+1..idx'last-1 loop
         -- wrk := acc*crs(j-1);
          pr := rcrs(j-1); pi := icrs(j-1);
          zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
          idx2 := idx(j); factor := Create(xps(idx2));
         -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
          ryd(idx2) := ryd(idx2) + factor*zr;
          iyd(idx2) := iyd(idx2) + factor*zi;
        end loop;
       -- wrk := acc*fwd(idx'last-2);
        idx2 := idx(idx'last);
        pr := rfwd(idx'last-2); pi := ifwd(idx'last-2);
        zr := racc*pr - iacc*pi; zi := racc*pi + iacc*pr;
        factor := Create(xps(idx2));
       -- wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
        ryd(idx2) := ryd(idx2) + factor*zr;
        iyd(idx2) := iyd(idx2) + factor*zi;
      end if; -- if idx(k)'last is 2
    end if; -- if idx(k)'last is 1
  end Speel;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    rcf : in Standard_Floating_Vectors.Vector;
                    icf : in Standard_Floating_Vectors.Vector;
                    rcst,icst : in double_float;
                    xr : in Standard_Floating_Vectors.Link_to_Vector;
                    xi : in Standard_Floating_Vectors.Link_to_Vector;
                    ryd : in Standard_Floating_Vectors.Link_to_Vector;
                    iyd : in Standard_Floating_Vectors.Link_to_Vector;
                    rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                    ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                    rbck : in Standard_Floating_Vectors.Link_to_Vector;
                    ibck : in Standard_Floating_Vectors.Link_to_Vector;
                    rcrs : in Standard_Floating_Vectors.Link_to_Vector;
                    icrs : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_VecVecs.VecVec;
                    ipwt : in Standard_Floating_VecVecs.VecVec;
                    hrp : in Standard_Floating_VecVecs.VecVec;
                    hip : in Standard_Floating_VecVecs.VecVec ) is

    dim : constant integer32 := xr'last;
    idx1,size : integer32;
    rkcff,ikcff : double_float;
    zr,zi,pr,pi : double_float;
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

    use Standard_Integer_Vectors;

  begin
    for i in 1..dim loop
      hrprow := hrp(i); hiprow := hip(i);
      for j in 1..dim loop
        hrprow(j) := 0.0; hiprow(j) := 0.0;
      end loop;
    end loop;
    ryd(0) := rcst; iyd(0) := icst;
    for k in xps'range loop
      declare
        xpk : constant Standard_Integer_Vectors.Vector := xps(k).all;
        idk : constant Standard_Integer_Vectors.Vector := idx(k).all;
        fck : constant Standard_Integer_Vectors.Link_to_Vector := fac(k);
      begin
        rkcff := rcf(k); ikcff := icf(k);
        if fck = null then
          if idk'last = 1 then
            idx1 := idk(1);
           -- yd(0) := yd(0) + kcff*x(idx1);
            pr := xr(idx1); pi := xi(idx1);
            zr := pr*rkcff - pi*ikcff; zi := pr*ikcff + pi*rkcff;
            ryd(0) := ryd(0) + zr; iyd(0) := iyd(0) + zi;
           -- yd(idx1) := yd(idx1) + kcff;
            ryd(idx1) := ryd(idx1) + rkcff; iyd(idx1) := iyd(idx1) + ikcff;
          else
            Indexed_Speel(idk,rkcff,ikcff,xr,xi,ryd,iyd,rfwd,ifwd,
                          rbck,ibck,rcrs,icrs,hrp,hip);
          end if;
        else
          Speel(xps(k),fck,idk,rkcff,ikcff,xr,xi,ryd,iyd,
                rfwd,ifwd,rbck,ibck,rcrs,icrs,rpwt,ipwt);
          size := idk'last;   -- number of participating variables
          if size = 1 then    -- one variable special case
            Standard_Hessian_Updaters.Speel1
              (hrp,hip,rkcff,ikcff,xpk,idk,fck.all,xr.all,xi.all,rpwt,ipwt);
          elsif size = 2 then -- two variable special case
            Standard_Hessian_Updaters.Speel2
              (hrp,hip,rkcff,ikcff,xpk,idk,fck.all,xr.all,xi.all,rpwt,ipwt);
          elsif size = 3 then -- three variable special case
            Standard_Hessian_Updaters.Speel3
              (hrp,hip,rkcff,ikcff,xpk,idk,fck.all,xr.all,xi.all,rpwt,ipwt);
          elsif size = 4 then -- four variable special case
            Standard_Hessian_Updaters.Speel4
              (hrp,hip,rkcff,ikcff,xpk,idk,fck.all,xr.all,xi.all,rpwt,ipwt);
          else
            Standard_Hessian_Updaters.SpeelN
              (hrp,hip,rkcff,ikcff,xpk,idk,fck.all,xr.all,xi.all,
               rfwd,ifwd,rbck,ibck,rpwt,ipwt);
          end if;
        end if;
      end;
    end loop;
    for i in 2..dim loop
      hrprow := hrp(i); hiprow := hip(i);
      for j in 1..(i-1) loop
        hrprow(j) := hrp(j)(i); hiprow(j) := hip(j)(i);
      end loop;
    end loop;
  end Speel;

-- AUXILIARY PROCEDURES :

  procedure Forward ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                      xi : in Standard_Floating_Vectors.Link_to_Vector;
                      fr : in Standard_Floating_Vectors.Link_to_Vector;
                      fi : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
  end Forward;

  procedure Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr := xr(xr'last); pi := xi(xr'last);
    idx := xi'last-1;
    qr := xr(idx); qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    br(1) := zr; bi(1) := zi;
    for k in 2..xr'last-2 loop
     -- b(k) := b(k-1)*x(x'last-k);
      pr := zr; pi := zi;
      idx := xr'last-k;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      br(k) := zr; bi(k) := zi;
    end loop;
  end Forward_Backward;

  procedure Fused_Forward_Backward
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr1,zi1,pr1,pi1,qr1,qi1 : double_float;
    zr2,zi2,pr2,pi2,qr2,qi2 : double_float;
    idx1,idx2 : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr1 := xr(1); pi1 := xi(1);
    idx1 := xr'first+1;
    qr1 := xr(idx1);  qi1 := xi(idx1);
    zr1 := pr1*qr1 - pi1*qi1;
    zi1 := pr1*qi1 + pi1*qr1;
    fr(1) := zr1; fi(1) := zi1;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr2 := xr(xr'last); pi2 := xi(xr'last);
    idx2 := xi'last-1;
    qr2 := xr(idx2); qi2 := xi(idx2);
    zr2 := pr2*qr2 - pi2*qi2;
    zi2 := pr2*qi2 + pi2*qr2;
    br(1) := zr2; bi(1) := zi2;
    for k in 2..dim-1 loop 
     -- f(k) := f(k-1)*x(k+1);
      pr1 := zr1; pi1 := zi1;
      idx1 := k+1;
      qr1 := xr(idx1); qi1 := xi(idx1);
      zr1 := pr1*qr1 - pi1*qi1;
      zi1 := pr1*qi1 + pi1*qr1;
      fr(k) := zr1; fi(k) := zi1;
     -- b(k) := b(k-1)*x(x'last-k);
      pr2 := zr2; pi2 := zi2;
      idx2 := xr'last-k;
      qr2 := xr(idx2); qi2 := xi(idx2);
      zr2 := pr2*qr2 - pi2*qi2;
      zi2 := pr2*qi2 + pi2*qr2;
      br(k) := zr2; bi(k) := zi2;
    end loop;
    if dim > 1 then
      pr1 := zr1; pi1 := zi1;
      idx1 := dim+1;
      qr1 := xr(idx1); qi1 := xi(idx1);
      zr1 := pr1*qr1 - pi1*qi1;
      zi1 := pr1*qi1 + pi1*qr1;
      fr(dim) := zr1; fi(dim) := zi1;
    end if;
  end Fused_Forward_Backward;

  procedure Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr,zi,pr,pi,qr,qi : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;

  begin
   -- f(f'first) := x(x'first)*x(x'first+1);
    pr := xr(1); pi := xi(1);
    idx := xr'first+1;
    qr := xr(idx);  qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    fr(1) := zr; fi(1) := zi;
    for k in 2..dim loop 
     -- f(k) := f(k-1)*x(k+1);
      pr := zr; pi := zi;
      idx := k+1;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      fr(k) := zr; fi(k) := zi;
    end loop;
   -- b(b'first) := x(x'last)*x(x'last-1);
    pr := xr(xr'last); pi := xi(xr'last);
    idx := xi'last-1;
    qr := xr(idx); qi := xi(idx);
    zr := pr*qr - pi*qi;
    zi := pr*qi + pi*qr;
    br(1) := zr; bi(1) := zi;
    for k in 2..xr'last-2 loop
     -- b(k) := b(k-1)*x(x'last-k);
      pr := zr; pi := zi;
      idx := xr'last-k;
      qr := xr(idx); qi := xi(idx);
      zr := pr*qr - pi*qi;
      zi := pr*qi + pi*qr;
      br(k) := zr; bi(k) := zi;
    end loop;
    if xr'last > 2 then
      if xr'last = 3 then
       -- c(1) := x(1)*x(3)
        pr := xr(1); pi := xi(1);
        qr := xr(3); qi := xi(3);
        zr := pr*qr - pi*qi;
        zi := pr*qi + pi*qr;
        cr(1) := zr; ci(1) := zi;
      else
       -- c(1) := x(1)*b(x'last-3);
        pr := xr(1); pi := xi(1);
        idx := xr'last-3;
        qr := br(idx); qi := bi(idx);
        zr := pr*qr - pi*qi;
        zi := pr*qi + pi*qr;
        cr(1) := zr; ci(1) := zi;
        for k in 2..xr'last-3 loop
         -- c(k) := f(k-1)*b(x'last-2-k);
          idx := k-1;
          pr := fr(idx); pi := fi(idx);
          idx := xr'last-2-k;
          qr := br(idx); qi := bi(idx);
          zr := pr*qr - pi*qi;
          zi := pr*qi + pi*qr;
          cr(k) := zr; ci(k) := zi;
        end loop;
        pr := xr(xr'last); pi := xi(xi'last);
        idx := xr'last-3;
        qr := fr(idx); qi := fi(idx);
        zr := pr*qr - pi*qi;
        zi := pr*qi + pi*qr;
        idx := xr'last-2;
        cr(idx) := zr; ci(idx) := zi;
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr1,zi1,pr1,pi1,qr1,qi1 : double_float;
    zr2,zi2,pr2,pi2,qr2,qi2 : double_float;
    zr3,zi3,pr3,pi3,qr3,qi3 : double_float;
    idx : integer32;
    dim : constant integer32 := xr'last-1;
    firstend,lastend,plusidx,minidx : integer32;

  begin
    if xr'last >= 8 then
      if xr'last mod 2 = 0 then
        lastend := xr'last-4;
        firstend := lastend/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        pr1 := xr(1); pi1 := xi(1);
        idx := xr'first+1;
        qr1 := xr(idx);  qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        pr2 := xr(xr'last); pi2 := xi(xr'last);
        idx := xi'last-1;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          idx := plusidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-plusidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
          plusidx := plusidx + 1;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          idx := minidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-minidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
          minidx := minidx - 1;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        idx := dim+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        pr3 := xr(1); pi3 := xi(1);
        idx := xr'last-3;
        qr3 := br(idx); qi3 := bi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        pr3 := xr(xr'last); pi3 := xi(xi'last);
        idx := xr'last-3;
        qr3 := fr(idx); qi3 := fi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        idx := xr'last-2;
        cr(idx) := zr3; ci(idx) := zi3;
      else
        lastend := xr'last-4;
        firstend := (xr'last-3)/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        pr1 := xr(1); pi1 := xi(1);
        idx := xr'first+1;
        qr1 := xr(idx);  qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        pr2 := xr(xr'last); pi2 := xi(xr'last);
        idx := xi'last-1;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        plusidx := firstend+1;
       -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
       -- idx := plusidx-1;
       -- pr3 := fr(idx); pi3 := fi(idx);
       -- f(plusidx-1) is still in zr1 and zi1
        idx := xr'last-2-plusidx;
        qr3 := br(idx); qi3 := bi(idx);
       -- zr3 := pr3*qr3 - pi3*qi3;
       -- zi3 := pr3*qi3 + pi3*qr3;
        zr3 := zr1*qr3 - zi1*qi3;
        zi3 := zr1*qi3 + zi1*qr3;
        cr(plusidx) := zr3; ci(plusidx) := zi3;
        minidx := plusidx;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          idx := k+1;
          qr1 := xr(idx); qi1 := xi(idx);
          zr1 := pr1*qr1 - pi1*qi1;
          zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          idx := xr'last-k;
          qr2 := xr(idx); qi2 := xi(idx);
          zr2 := pr2*qr2 - pi2*qi2;
          zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          idx := plusidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-plusidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
          idx := minidx-1;
          pr3 := fr(idx); pi3 := fi(idx);
          idx := xr'last-2-minidx;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        idx := plusidx+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-plusidx;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        idx := dim+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        pr3 := xr(1); pi3 := xi(1);
        idx := xr'last-3;
        qr3 := br(idx); qi3 := bi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        pr3 := xr(xr'last); pi3 := xi(xi'last);
        idx := xr'last-3;
        qr3 := fr(idx); qi3 := fi(idx);
        zr3 := pr3*qr3 - pi3*qi3;
        zi3 := pr3*qi3 + pi3*qr3;
        idx := xr'last-2;
        cr(idx) := zr3; ci(idx) := zi3;
      end if;
    else
     -- f(f'first) := x(x'first)*x(x'first+1);
      pr1 := xr(1); pi1 := xi(1);
      idx := xr'first+1;
      qr1 := xr(idx);  qi1 := xi(idx);
      zr1 := pr1*qr1 - pi1*qi1;
      zi1 := pr1*qi1 + pi1*qr1;
      fr(1) := zr1; fi(1) := zi1;
     -- b(b'first) := x(x'last)*x(x'last-1);
      pr2 := xr(xr'last); pi2 := xi(xr'last);
      idx := xi'last-1;
      qr2 := xr(idx); qi2 := xi(idx);
      zr2 := pr2*qr2 - pi2*qi2;
      zi2 := pr2*qi2 + pi2*qr2;
      br(1) := zr2; bi(1) := zi2;
      for k in 2..xr'last-2 loop 
       -- f(k) := f(k-1)*x(k+1);
        pr1 := zr1; pi1 := zi1;
        idx := k+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(k) := zr1; fi(k) := zi1;
       -- b(k) := b(k-1)*x(x'last-k);
        pr2 := zr2; pi2 := zi2;
        idx := xr'last-k;
        qr2 := xr(idx); qi2 := xi(idx);
        zr2 := pr2*qr2 - pi2*qi2;
        zi2 := pr2*qi2 + pi2*qr2;
        br(k) := zr2; bi(k) := zi2;
      end loop;
      if dim > 1 then
        pr1 := zr1; pi1 := zi1;
        idx := dim+1;
        qr1 := xr(idx); qi1 := xi(idx);
        zr1 := pr1*qr1 - pi1*qi1;
        zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
      end if;
      if xr'last > 2 then
        if xr'last = 3 then
         -- c(1) := x(1)*x(3)
          pr3 := xr(1); pi3 := xi(1);
          qr3 := xr(3); qi3 := xi(3);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(1) := zr3; ci(1) := zi3;
        else
         -- c(1) := x(1)*b(x'last-3);
          pr3 := xr(1); pi3 := xi(1);
          idx := xr'last-3;
          qr3 := br(idx); qi3 := bi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          cr(1) := zr3; ci(1) := zi3;
          for k in 2..xr'last-3 loop
           -- c(k) := f(k-1)*b(x'last-2-k);
            idx := k-1;
            pr3 := fr(idx); pi3 := fi(idx);
            idx := xr'last-2-k;
            qr3 := br(idx); qi3 := bi(idx);
            zr3 := pr3*qr3 - pi3*qi3;
            zi3 := pr3*qi3 + pi3*qr3;
            cr(k) := zr3; ci(k) := zi3;
          end loop;
          pr3 := xr(xr'last); pi3 := xi(xi'last);
          idx := xr'last-3;
          qr3 := fr(idx); qi3 := fi(idx);
          zr3 := pr3*qr3 - pi3*qi3;
          zi3 := pr3*qi3 + pi3*qr3;
          idx := xr'last-2;
          cr(idx) := zr3; ci(idx) := zi3;
        end if;
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                fr : in Standard_Floating_Vectors.Link_to_Vector;
                fi : in Standard_Floating_Vectors.Link_to_Vector;
                br : in Standard_Floating_Vectors.Link_to_Vector;
                bi : in Standard_Floating_Vectors.Link_to_Vector;
                cr : in Standard_Floating_Vectors.Link_to_Vector;
                ci : in Standard_Floating_Vectors.Link_to_Vector ) is

    zr1,zi1,pr1,pi1,qr1,qi1 : double_float;
    zr2,zi2,pr2,pi2,qr2,qi2 : double_float;
    zr3,zi3,pr3,pi3,qr3,qi3 : double_float;
    ptr : integer32;
    dim : constant integer32 := idx'last-1;
    firstend,lastend,plusidx,minidx : integer32;

  begin
    if idx'last >= 8 then
      if idx'last mod 2 = 0 then
        lastend := idx'last-4;
        firstend := lastend/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        ptr := idx(1); pr1 := xr(ptr); pi1 := xi(ptr);
        ptr := idx(idx'first+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        ptr := idx(idx'last); pr2 := xr(ptr); pi2 := xi(ptr);
        ptr := idx(idx'last-1); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          ptr := plusidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-plusidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
          plusidx := plusidx+1;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          ptr := minidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-minidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
          minidx := minidx - 1;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(dim+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := br(ptr); qi3 := bi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        ptr := idx(idx'last); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := fr(ptr); qi3 := fi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        ptr := idx'last-2; cr(ptr) := zr3; ci(ptr) := zi3;
      else
        lastend := idx'last-4;
        firstend := (idx'last-3)/2;
       -- f(f'first) := x(x'first)*x(x'first+1);
        ptr := idx(1); pr1 := xr(ptr); pi1 := xi(ptr);
        ptr := idx(idx'first+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(1) := zr1; fi(1) := zi1;
       -- b(b'first) := x(x'last)*x(x'last-1);
        ptr := idx(idx'last); pr2 := xr(ptr); pi2 := xi(ptr);
        ptr := idx(idx'last-1); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..firstend loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        plusidx := firstend+1;
       -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
        ptr := plusidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
        ptr := idx'last-2-plusidx; qr3 := br(ptr); qi3 := bi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        cr(plusidx) := zr3; ci(plusidx) := zi3;
        minidx := plusidx;
        for k in firstend+1..lastend loop
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
         -- c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          ptr := plusidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-plusidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(plusidx) := zr3; ci(plusidx) := zi3;
         -- c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
          ptr := minidx-1; pr3 := fr(ptr); pi3 := fi(ptr);
          ptr := idx'last-2-minidx; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(minidx) := zr3; ci(minidx) := zi3;
        end loop;
        plusidx := lastend+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(plusidx) := f(plusidx-1)*x(plusidx+1);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(plusidx+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(plusidx) := zr1; fi(plusidx) := zi1;
       -- b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
        pr2 := zr2; pi2 := zi2;
        ptr := idx(idx'last-plusidx); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(plusidx) := zr2; bi(plusidx) := zi2;
        plusidx := plusidx+1;
       -- f(f'last) := f(f'last-1)*x(x'last);
        pr1 := zr1; pi1 := zi1;
        ptr := idx(dim+1); qr1 := xr(ptr); qi1 := xi(ptr);
        zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
        fr(dim) := zr1; fi(dim) := zi1;
       -- c(1) := x(1)*b(x'last-3);
        ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := br(ptr); qi3 := bi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        cr(1) := zr3; ci(1) := zi3;
       -- c(x'last-2) := x(x'last)*f(x'last-3);
        ptr := idx(idx'last); pr3 := xr(ptr); pi3 := xi(ptr);
        ptr := idx'last-3; qr3 := fr(ptr); qi3 := fi(ptr);
        zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
        ptr := idx'last-2;
        cr(ptr) := zr3; ci(ptr) := zi3;
      end if;
    else
     -- f(f'first) := x(x'first)*x(x'first+1);
      ptr := idx(1); pr1 := xr(ptr); pi1 := xi(ptr);
      ptr := idx(idx'first+1); qr1 := xr(ptr); qi1 := xi(ptr);
      zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
      fr(1) := zr1; fi(1) := zi1;
      if idx'last > 2 then
       -- b(b'first) := x(x'last)*x(x'last-1);
        ptr := idx(idx'last); pr2 := xr(ptr); pi2 := xi(ptr);
        ptr := idx(idx'last-1); qr2 := xr(ptr); qi2 := xi(ptr);
        zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
        br(1) := zr2; bi(1) := zi2;
        for k in 2..idx'last-2 loop 
         -- f(k) := f(k-1)*x(k+1);
          pr1 := zr1; pi1 := zi1;
          ptr := idx(k+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(k) := zr1; fi(k) := zi1;
         -- b(k) := b(k-1)*x(x'last-k);
          pr2 := zr2; pi2 := zi2;
          ptr := idx(idx'last-k); qr2 := xr(ptr); qi2 := xi(ptr);
          zr2 := pr2*qr2 - pi2*qi2; zi2 := pr2*qi2 + pi2*qr2;
          br(k) := zr2; bi(k) := zi2;
        end loop;
        if dim > 1 then
          pr1 := zr1; pi1 := zi1;
          ptr := idx(dim+1); qr1 := xr(ptr); qi1 := xi(ptr);
          zr1 := pr1*qr1 - pi1*qi1; zi1 := pr1*qi1 + pi1*qr1;
          fr(dim) := zr1; fi(dim) := zi1;
        end if;
        if idx'last = 3 then
         -- c(1) := x(1)*x(3)
          ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
          ptr := idx(3); qr3 := xr(ptr); qi3 := xi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          cr(1) := zr3; ci(1) := zi3;
        else
         -- c(1) := x(1)*b(x'last-3);
          ptr := idx(1); pr3 := xr(ptr); pi3 := xi(ptr);
          ptr := idx'last-3; qr3 := br(ptr); qi3 := bi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
         -- put("Assigning "); put(zr3); put_line(" to cr(1)");
         -- put("Assigning "); put(zi3); put_line(" to ci(1)");
          cr(1) := zr3; ci(1) := zi3;
          for k in 2..idx'last-3 loop
           -- c(k) := f(k-1)*b(x'last-2-k);
            ptr := k-1; pr3 := fr(ptr); pi3 := fi(ptr);
            ptr := idx'last-2-k; qr3 := br(ptr); qi3 := bi(ptr);
            zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
            cr(k) := zr3; ci(k) := zi3;
          end loop;
          ptr := idx(idx'last); pr3 := xr(ptr); pi3 := xi(ptr);
          ptr := idx'last-3; qr3 := fr(ptr); qi3 := fi(ptr);
          zr3 := pr3*qr3 - pi3*qi3; zi3 := pr3*qi3 + pi3*qr3;
          ptr := idx'last-2; cr(ptr) := zr3; ci(ptr) := zi3;
        end if;
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(mxe'range);

  begin
    for k in mxe'range loop
      if mxe(k) > 1 then
        res(k) := new Standard_Floating_Vectors.Vector'(1..mxe(k)-1 => 0.0);
      end if;
    end loop;
    return res;
  end Allocate;

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec ) is

    rlnk : Standard_Floating_Vectors.Link_to_Vector;
    ilnk : Standard_Floating_Vectors.Link_to_Vector;
    zr,zi,yr,yi,xrk,xik : double_float;

  begin
    for k in xr'range loop
      if mxe(k) > 1 then
        rlnk := rpwt(k); ilnk := ipwt(k);
       -- lnk(1) := x(k)*x(k);
        xrk := xr(k); xik := xi(k);
        zr := xrk*xrk - xik*xik;
        zi := 2.0*xrk*xik;
        rlnk(1) := zr; ilnk(1) := zi;
        for i in 2..mxe(k)-1 loop
         -- lnk(i) := lnk(i-1)*x(k);
          yr := zr; yi := zi;
          zr := xrk*yr - xik*yi;
          zi := xrk*yi + xik*yr;
          rlnk(i) := zr; ilnk(i) := zi;
        end loop;
      end if;
    end loop;
  end Power_Table;

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                rcf : in double_float;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                rpf : out double_float ) is

    rpwx : Standard_Floating_Vectors.Link_to_Vector;
    idx,powidx : integer32;

  begin
    idx := fac(fac'first); powidx := xps(idx);
    if powidx = 2 then -- res := cff*x(idx);
      rpf := rcf*xr(idx);
    else
      rpwx := rpwt(idx); -- res := cff*pwx(powidx-2);
      rpf := rcf*rpwx(powidx-2);
    end if;
    for k in fac'first+1..fac'last loop
      idx := fac(k); powidx := xps(idx);
      if powidx = 2 then -- res := res*x(idx);
        rpf := rpf*xr(idx);
      else
        rpwx := rpwt(idx); -- res := res*pwx(powidx-2);
        rpf := rpf*rpwx(powidx-2);
      end if;
    end loop;
  end Multiply_Factor;

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                rcf,icf : in double_float;
                rpwt : in Standard_Floating_VecVecs.VecVec;
                ipwt : in Standard_Floating_VecVecs.VecVec;
                rpf,ipf : out double_float ) is

    rpwx : Standard_Floating_Vectors.Link_to_Vector;
    ipwx : Standard_Floating_Vectors.Link_to_Vector;
    idx,powidx : integer32;
    zr,zi,xrk,xik : double_float;

  begin
    idx := fac(fac'first); powidx := xps(idx);
    if powidx = 2 then
     -- res := cff*x(idx);
      xrk := xr(idx); xik := xi(idx);
      zr := xrk*rcf - xik*icf;
      zi := xrk*icf + xik*rcf;
      rpf := zr; ipf := zi;
    else
      rpwx := rpwt(idx); ipwx := ipwt(idx);
     -- res := cff*pwx(powidx-2);
      xrk := rpwx(powidx-2); xik := ipwx(powidx-2);
      zr := xrk*rcf - xik*icf;
      zi := xrk*icf + xik*rcf;
      rpf := zr; ipf := zi;
    end if;
    for k in fac'first+1..fac'last loop
      idx := fac(k); powidx := xps(idx);
      if powidx = 2 then
       -- res := res*x(idx);
        xrk := xr(idx); xik := xi(idx);
        zr := xrk*rpf - xik*ipf;
        zi := xrk*ipf + xik*rpf;
        rpf := zr; ipf := zi;
      else
        rpwx := rpwt(idx); ipwx := ipwt(idx);
       -- res := res*pwx(powidx-2);
        xrk := rpwx(powidx-2); xik := ipwx(powidx-2);
        zr := xrk*rpf - xik*ipf;
        zi := xrk*ipf + xik*rpf;
        rpf := zr; ipf := zi;
      end if;
    end loop;
  end Multiply_Factor;

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit ) is
  begin
    c.dim := -1;
    c.pdg := -1;
    Standard_Floating_Vectors.Clear(c.rfwd);
    Standard_Floating_Vectors.Clear(c.ifwd);
    Standard_Floating_Vectors.Clear(c.rbck);
    Standard_Floating_Vectors.Clear(c.ibck);
    Standard_Floating_Vectors.Clear(c.rcrs);
    Standard_Floating_Vectors.Clear(c.icrs);
  end Clear;

  procedure Clear ( c : in out Link_to_Circuit ) is

    procedure free is new unchecked_deallocation(Circuit,Link_to_Circuit);

  begin
    if c /= null then
      Clear(c.all);
      free(c);
    end if;
  end Clear;

  procedure Clear ( c : in out Circuits ) is
  begin
    for k in c'range loop
      Clear(c(k));
    end loop;
  end Clear;

  procedure Clear ( s : in out System ) is

    use Standard_Floating_VecVecs;

  begin
    Clear(s.crc);
    Standard_Floating_Vectors.Clear(s.ryd);
    Standard_Floating_Vectors.Clear(s.iyd);
    Standard_Floating_VecVecs.Clear(s.rpwt);
    Standard_Floating_VecVecs.Clear(s.ipwt);
    if s.jrc /= null
     then Standard_Floating_VecVecs.Deep_Clear(s.jrc);
    end if;
    if s.jic /= null
     then Standard_Floating_VecVecs.Deep_Clear(s.jic);
    end if;
    if s.hrp /= null
     then Standard_Floating_VecVecs.Deep_Clear(s.hrp);
    end if;
    if s.hip /= null
     then Standard_Floating_VecVecs.Deep_Clear(s.hip);
    end if;
  end Clear;

  procedure Clear ( s : in out Link_to_System ) is

    procedure free is new unchecked_deallocation(System,Link_to_System);

  begin
    if s /= null then
      Clear(s.all);
      free(s);
    end if;
  end Clear;

  procedure Clear ( s : in out System_Array ) is
  begin
    for k in s'range loop
      Clear(s(k));
    end loop;
  end Clear;

end Standard_Coefficient_Circuits;
