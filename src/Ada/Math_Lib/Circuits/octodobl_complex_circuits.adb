with unchecked_deallocation;
with Octo_Double_Numbers;                 use Octo_Double_Numbers;
with OctoDobl_Complex_Singular_Values;
with Exponent_Indices;
with OctoDobl_Hessian_Updaters;

package body OctoDobl_Complex_Circuits is

-- CONSTRUCTORS :

  function Allocate ( nbr,dim : integer32 ) return Circuit is

    res : Circuit(nbr);
    zero : constant OctoDobl_Complex_Numbers.Complex_Number
         := OctoDobl_Complex_Numbers.Create(integer(0));
    forward : constant OctoDobl_Complex_Vectors.Vector(1..dim-1)
            := (1..dim-1 => zero);
    backward : constant OctoDobl_Complex_Vectors.Vector(1..dim-2)
             := (1..dim-2 => zero);
    cross : constant OctoDobl_Complex_Vectors.Vector(1..dim-2)
          := (1..dim-2 => zero);

  begin
    res.dim := dim;
    res.pdg := -1;
    res.fwd := new OctoDobl_Complex_Vectors.Vector'(forward);
    res.bck := new OctoDobl_Complex_Vectors.Vector'(backward);
    res.crs := new OctoDobl_Complex_Vectors.Vector'(cross);
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
    zero : constant OctoDobl_Complex_Numbers.Complex_Number 
         := OctoDobl_Complex_Numbers.Create(integer(0));

  begin
    res.crc := c;
    res.mxe := Exponent_Maxima(c,dim);
    res.pwt := Allocate(res.mxe);
    res.yd := new OctoDobl_Complex_Vectors.Vector'(0..dim => zero);
    return res;
  end Create;

  function Allocate ( neq,dim : integer32 )
                    return OctoDobl_Complex_VecMats.VecMat is

    res : OctoDobl_Complex_VecMats.VecMat(1..neq);

  begin
    for k in 1..neq loop
      declare
        mat : OctoDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := OctoDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        res(k) := new OctoDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate;

-- SINGULAR VALUE DECOMPOSITIONS :

  procedure Singular_Values
              ( s : in out System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                vh : in OctoDobl_Complex_VecMats.VecMat;
                U : out OctoDobl_Complex_Matrices.Matrix;
                V : out OctoDobl_Complex_Matrices.Matrix;
                e : out OctoDobl_Complex_Vectors.Vector;
                svls : in OctoDobl_Complex_VecVecs.VecVec ) is

    info : integer32;

  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff2(s.crc,x,s.yd,s.pwt,s.fx,s.jm,vh);
    OctoDobl_Complex_Singular_Values.SVD
      (s.jm,s.dim,s.dim,svls(0).all,e,U,V,0,info);
    for k in vh'range loop
      OctoDobl_Complex_Singular_Values.SVD
        (vh(k).all,s.dim,s.dim,svls(k).all,e,U,V,0,info);
    end loop;
  end Singular_Values;

  procedure Singular_Values
              ( s : in Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                vh : in OctoDobl_Complex_VecMats.VecMat;
                U : out OctoDobl_Complex_Matrices.Matrix;
                V : out OctoDobl_Complex_Matrices.Matrix;
                e : out OctoDobl_Complex_Vectors.Vector;
                svls : in OctoDobl_Complex_VecVecs.VecVec ) is
  begin
    if s /= null
     then Singular_Values(s.all,x,vh,U,V,e,svls);
    end if;
  end Singular_Values;

  procedure Singular_Values
              ( c : in Circuit;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                pwt : in OctoDobl_Complex_VecVecs.VecVec;
                A : out OctoDobl_Complex_Matrices.Matrix;
                U : out OctoDobl_Complex_Matrices.Matrix;
                V : out OctoDobl_Complex_Matrices.Matrix;
                e : out OctoDobl_Complex_Vectors.Vector;
                s : out OctoDobl_Complex_Vectors.Vector ) is

    info : integer32;

  begin
    if c.pdg <= 1 then
      s := (s'range => OctoDobl_Complex_Numbers.Create(integer32(0)));
    else
      Speel(c,x,yd,pwt,A);
      OctoDobl_Complex_Singular_Values.SVD(A,c.dim,c.dim,s,e,U,V,0,info);
    end if;
  end Singular_Values;

  procedure Singular_Values
              ( c : in Link_to_Circuit;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                pwt : in OctoDobl_Complex_VecVecs.VecVec;
                A : out OctoDobl_Complex_Matrices.Matrix;
                U : out OctoDobl_Complex_Matrices.Matrix;
                V : out OctoDobl_Complex_Matrices.Matrix;
                e : out OctoDobl_Complex_Vectors.Vector;
                s : out OctoDobl_Complex_Vectors.Vector ) is
  begin
    if c /= null
     then Singular_Values(c.all,x,yd,pwt,A,U,V,e,s);
    end if;
  end Singular_Values;

-- ALGORITHMIC DIFFERENTIATION AND EVALUATION OF CIRCUITS :

  procedure EvalDiff
              ( s : in out System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff(s.crc,x,s.yd,s.pwt,s.fx,s.jm);
  end EvalDiff;

  procedure EvalDiff
              ( s : in Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff(s.crc,x,s.yd,s.pwt,s.fx,s.jm);
  end EvalDiff;

  procedure EvalDiff2
              ( s : in out System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                vh : in OctoDobl_Complex_VecMats.VecMat ) is
  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff2(s.crc,x,s.yd,s.pwt,s.fx,s.jm,vh);
  end EvalDiff2;

  procedure EvalDiff2
              ( s : in Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                vh : in OctoDobl_Complex_VecMats.VecMat ) is
  begin
    Power_Table(s.mxe,x,s.pwt);
    EvalDiff2(s.crc,x,s.yd,s.pwt,s.fx,s.jm,vh);
  end EvalDiff2;

  procedure EvalDiff
              ( c : in Circuits;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                pwt : in OctoDobl_Complex_VecVecs.VecVec;
                fx : out OctoDobl_Complex_Vectors.Vector;
                jm : out OctoDobl_Complex_Matrices.Matrix ) is
  begin
    for i in c'range loop
      Speel(c(i).all,x,yd,pwt);
      fx(i) := yd(0);
      for j in jm'range(2) loop
        jm(i,j) := yd(j);
        yd(j) := OctoDobl_Complex_Numbers.Create(integer(0));
      end loop;
    end loop;
  end EvalDiff;

  procedure EvalDiff2
              ( c : in Circuits;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                pwt : in OctoDobl_Complex_VecVecs.VecVec;
                fx : out OctoDobl_Complex_Vectors.Vector;
                jm : out OctoDobl_Complex_Matrices.Matrix;
                vh : in OctoDobl_Complex_VecMats.VecMat ) is

    mat : OctoDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in c'range loop
      mat := vh(i);
      Speel(c(i).all,x,yd,pwt,mat.all);
      fx(i) := yd(0);
      for j in jm'range(2) loop
        jm(i,j) := yd(j);
        yd(j) := OctoDobl_Complex_Numbers.Create(integer(0));
      end loop;
    end loop;
  end EvalDiff2;

-- ALGORITMIC DIFFERENTIATION AND EVALUATION OF CIRCUIT :

  procedure Indexed_Speel
              ( c : in Circuit;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                h : out OctoDobl_Complex_Matrices.Matrix ) is
  begin
    Indexed_Speel(c.xps,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs,h);
  end Indexed_Speel;

  procedure Indexed_Speel
              ( c : in Circuit;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Indexed_Speel(c.xps,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs);
  end Indexed_Speel;

  procedure Speel ( c : in Circuit;
                    x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in OctoDobl_Complex_VecVecs.VecVec ) is
  begin
    Speel(c.xps,c.idx,c.fac,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs,pwt);
  end Speel;

  procedure Speel ( c : in Circuit;
                    x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in OctoDobl_Complex_VecVecs.VecVec;
                    h : out OctoDobl_Complex_Matrices.Matrix ) is
  begin
    Speel(c.xps,c.idx,c.fac,c.cff,c.cst,x,yd,c.fwd,c.bck,c.crs,pwt,h);
  end Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                cff : in OctoDobl_Complex_Vectors.Vector;
                cst : in OctoDobl_Complex_Numbers.Complex_Number;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                fwd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                bck : in OctoDobl_Complex_Vectors.Link_to_Vector;
                crs : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    kcff : OctoDobl_Complex_Numbers.Complex_Number;

    use OctoDobl_Complex_Numbers,Standard_Integer_Vectors;

  begin
    yd(0) := cst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        if idk'last = 1 then
          idx1 := idk(1); kcff := cff(k);
          yd(0) := yd(0) + kcff*x(idx1);
          yd(idx1) := yd(idx1) + kcff;
        else
          Forward_Backward_Cross(idk.all,x,fwd,bck,crs);
          kcff := cff(k); yd(0) := yd(0) + kcff*fwd(idk'last-1); 
          if idk'last = 2 then
            idx1 := idk(1); idx2 := idk(2);
            yd(idx2) := yd(idx2) + kcff*x(idx1);
            yd(idx1) := yd(idx1) + kcff*x(idx2);
          else -- idk'last > 2
            idx1 := idk(1);
            yd(idx1) := yd(idx1) + kcff*bck(idk'last-2);
            for j in idk'first+1..idk'last-1 loop
              idx1 := idk(j);
              yd(idx1) := yd(idx1) + kcff*crs(j-1);
            end loop;
            idx1 := idk(idk'last);
            yd(idx1) := yd(idx1) + kcff*fwd(idk'last-2);
          end if;
        end if;
      end if;
    end loop;
  end Indexed_Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_Vectors.Vector;
                cff : in OctoDobl_Complex_Numbers.Complex_Number;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                fwd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                bck : in OctoDobl_Complex_Vectors.Link_to_Vector;
                crs : in OctoDobl_Complex_Vectors.Link_to_Vector;
                h : in out OctoDobl_Complex_Matrices.Matrix ) is

    use OctoDobl_Complex_Numbers;

    sz : constant integer32 := idx'last;
    idx1,idx2,idx3,idx4,idx5,idx6 : integer32;
    acc : Complex_Number;

  begin
    Forward_Backward_Cross(idx,x,fwd,bck,crs);
    yd(0) := yd(0) + cff*fwd(sz-1); 
    if sz = 2 then
      idx1 := idx(1); idx2 := idx(2);
      yd(idx2) := yd(idx2) + cff*x(idx1);
      yd(idx1) := yd(idx1) + cff*x(idx2);
      h(idx1,idx2) := h(idx1,idx2) + cff;
    else -- sz > 2
      idx1 := idx(1);
      yd(idx1) := yd(idx1) + cff*bck(sz-2);
      for j in idx'first+1..sz-1 loop
        idx1 := idx(j); yd(idx1) := yd(idx1) + cff*crs(j-1);
      end loop;
      idx1 := idx(sz); yd(idx1) := yd(idx1) + cff*fwd(sz-2);
      if sz = 3 then
        idx1 := idx(1); idx2 := idx(2); idx3 := idx(3);
        h(idx1,idx2) := h(idx1,idx2) + cff*x(idx3);
        h(idx1,idx3) := h(idx1,idx3) + cff*x(idx2);
        h(idx2,idx3) := h(idx2,idx3) + cff*x(idx1);
      else -- sz > 3
        idx5 := idx(sz-1); idx6 := idx(sz); idx3 := sz-3;
       -- last element is copy of fwd(sz-3), multiplied with cff
        h(idx5,idx6) := h(idx5,idx6) + cff*fwd(idx3);
       -- first element is copy of bck(sz-3), multiplied with cff
        idx1 := idx(1); idx2 := idx(2);
        h(idx1,idx2) := h(idx1,idx2) + cff*bck(idx3);
        if sz = 4 then -- special case for all rows
          idx3 := idx(3); idx4 := idx(4);
          acc := cff*x(idx2);
          h(idx1,idx3) := h(idx1,idx3) + acc*x(idx6);
          h(idx1,idx4) := h(idx1,idx4) + acc*x(idx5);
          acc := cff*x(idx1);
          h(idx2,idx3) := h(idx2,idx3) + acc*x(idx6);
          h(idx2,idx4) := h(idx2,idx4) + acc*x(idx5);
        else -- sz > 4
         -- first row is special, starts with x(idx(2)) after diagonal
          idx3 := idx(3); idx4 := sz-4;
          acc := cff*x(idx2);
          h(idx1,idx3) := h(idx1,idx3) + acc*bck(idx4);
          for k in 4..sz-2 loop
            idx4 := idx(k-1); idx5 := idx(k); idx6 := sz-k-1;
            acc := acc*x(idx4);
            h(idx1,idx5) := h(idx1,idx5) + acc*bck(idx6);
          end loop;
          idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
          acc := acc*x(idx4);
          h(idx1,idx5) := h(idx1,idx5) + acc*x(idx6);
          h(idx1,idx6) := h(idx1,idx6) + acc*x(idx5);
         -- second row is special, starts with x(idx(1)) after diagonal
          acc := cff*x(idx1); idx4 := sz-4;
          h(idx2,idx3) := h(idx2,idx3) + acc*bck(idx4);
          for k in 4..sz-2 loop
            idx4 := idx(k-1); idx5 := idx(k); idx6 := sz-k-1;
            acc := acc*x(idx4);
            h(idx2,idx5) := h(idx2,idx5) + acc*bck(idx6);
          end loop;
          idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
          acc := acc*x(idx4);
          h(idx2,idx5) := h(idx2,idx5) + acc*x(idx6);
          h(idx2,idx6) := h(idx2,idx6) + acc*x(idx5);
         -- the row with index sz-2 has a general formula
          idx3 := sz-4; acc := cff*fwd(idx3);
          h(idx4,idx5) := h(idx4,idx5) + acc*x(idx6);
          h(idx4,idx6) := h(idx4,idx6) + acc*x(idx5);
          for rw in 3..sz-3 loop  -- row rw starts with fwd(rw-2)
            idx1 := idx(rw); idx2 := idx(rw+1);
            acc := cff*fwd(rw-2);
            h(idx1,idx2) := h(idx1,idx2) + acc*bck(sz-rw-2);
            for k in rw+2..sz-2 loop
              idx2 := idx(k-1); idx3 := idx(k); idx4 := sz-k-1;
              acc := acc*x(idx2);
              h(idx1,idx3) := h(idx1,idx3) + acc*bck(idx4);
            end loop;
            idx4 := idx(sz-2); idx5 := idx(sz-1); idx6 := idx(sz);
            acc := acc*x(idx4);
            h(idx1,idx5) := h(idx1,idx5) + acc*x(idx6);
            h(idx1,idx6) := h(idx1,idx6) + acc*x(idx5);
          end loop;
        end if;
      end if;
    end if;
  end Indexed_Speel;

  procedure Indexed_Speel
              ( idx : in Standard_Integer_VecVecs.VecVec;
                cff : in OctoDobl_Complex_Vectors.Vector;
                cst : in OctoDobl_Complex_Numbers.Complex_Number;
                x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                fwd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                bck : in OctoDobl_Complex_Vectors.Link_to_Vector;
                crs : in OctoDobl_Complex_Vectors.Link_to_Vector;
                h : out OctoDobl_Complex_Matrices.Matrix ) is

    dim : constant integer32 := x'last;
    idk : Standard_Integer_Vectors.Link_to_Vector;
    idx1 : integer32;
    kcff : OctoDobl_Complex_Numbers.Complex_Number;

    use OctoDobl_Complex_Numbers,Standard_Integer_Vectors;

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        h(i,j) := Create(integer(0));
      end loop;
    end loop;
    yd(0) := cst;
    for k in idx'range loop
      idk := idx(k);
      if idk /= null then
        kcff := cff(k);
        if idk'last = 1 then
          idx1 := idk(1);
          yd(0) := yd(0) + kcff*x(idx1);
          yd(idx1) := yd(idx1) + kcff;
        else
          Indexed_Speel(idk.all,kcff,x,yd,fwd,bck,crs,h);
        end if;
      end if;
    end loop;
    for i in 2..dim loop
      for j in 1..(i-1) loop
        h(i,j) := h(j,i);
      end loop;
    end loop;
  end Indexed_Speel;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    cff : in OctoDobl_Complex_Vectors.Vector;
                    cst : in OctoDobl_Complex_Numbers.Complex_Number;
                    x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    fwd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    bck : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    crs : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in OctoDobl_Complex_VecVecs.VecVec ) is

    idk,fck,xpk : Standard_Integer_Vectors.Link_to_Vector;
    idx1,idx2 : integer32;
    kcff,acc,wrk : OctoDobl_Complex_Numbers.Complex_Number;
    factor : octo_double;

    use OctoDobl_Complex_Numbers,Standard_Integer_Vectors;

  begin
    yd(0) := cst;
    for k in idx'range loop
      idk := idx(k);
      kcff := cff(k);
      if idk = null then -- we have an extra constant coefficient
        yd(0) := yd(0) + kcff;
      else
        if idk'last < idk'first then -- deal with stray 0 exponent
          yd(0) := yd(0) + kcff;
        else
          idx1 := idk(1); fck := fac(k);
          if fck = null then
            if idk'last = 1 then
              yd(0) := yd(0) + kcff*x(idx1);
              yd(idx1) := yd(idx1) + kcff;
            else
              Forward_Backward_Cross(idk.all,x,fwd,bck,crs);
              yd(0) := yd(0) + kcff*fwd(idk'last-1); 
              if idk'last = 2 then
                idx2 := idk(2);
                yd(idx2) := yd(idx2) + kcff*x(idx1);
                yd(idx1) := yd(idx1) + kcff*x(idx2);
              else -- idk'last > 2
                yd(idx1) := yd(idx1) + kcff*bck(idk'last-2);
                for j in idk'first+1..idk'last-1 loop
                  idx2 := idk(j);
                  yd(idx2) := yd(idx2) + kcff*crs(j-1);
                end loop;
                idx2 := idk(idk'last);
                yd(idx2) := yd(idx2) + kcff*fwd(idk'last-2);
              end if;
            end if;
          else
            xpk := xps(k);
            if idk'last = 1 then
              Multiply_Factor(xpk,fck,x,kcff,pwt,acc);
              yd(0) := yd(0) + acc*x(idx1);
              factor := Octo_Double_Numbers.Create(xpk(idx1));
              yd(idx1) := yd(idx1) + factor*acc;
            else
              Forward_Backward_Cross(idk.all,x,fwd,bck,crs);
              Multiply_Factor(xpk,fck,x,kcff,pwt,acc);
              yd(0) := yd(0) + acc*fwd(idk'last-1); 
              if idk'last = 2 then
                idx2 := idk(2);
                wrk := acc*x(idx1);
                factor := Octo_Double_Numbers.Create(xpk(idx2));
                wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
                wrk := acc*x(idx2);
                factor := Octo_Double_Numbers.Create(xpk(idx1));
                wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
              else -- idk'last > 2
                wrk := acc*bck(idk'last-2);
                factor := Octo_Double_Numbers.Create(xpk(idx1));
                wrk := factor*wrk;
                yd(idx1) := yd(idx1) + wrk;
                for j in idk'first+1..idk'last-1 loop
                  wrk := acc*crs(j-1); idx2 := idk(j);
                  factor := Octo_Double_Numbers.Create(xpk(idx2));
                  wrk := factor*wrk;
                  yd(idx2) := yd(idx2) + wrk;
                end loop;
                wrk := acc*fwd(idk'last-2); idx2 := idk(idk'last);
                factor := Octo_Double_Numbers.Create(xpk(idx2));
                wrk := factor*wrk;
                yd(idx2) := yd(idx2) + wrk;
              end if;
            end if;
          end if;
        end if;
      end if;
    end loop;
  end Speel;

  procedure Speel ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                    idx : in Standard_Integer_Vectors.Vector;
                    cff : in OctoDobl_Complex_Numbers.Complex_Number;
                    x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    fwd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    bck : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    crs : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in OctoDobl_Complex_VecVecs.VecVec ) is

    use OctoDobl_Complex_Numbers,Standard_Integer_Vectors;

    idx1,idx2 : integer32;
    acc,wrk : Complex_Number;
    factor : octo_double;

  begin
    idx1 := idx(1);
    if idx'last = 1 then
      Multiply_Factor(xps,fac,x,cff,pwt,acc);
      yd(0) := yd(0) + acc*x(idx1);
      factor := Octo_Double_Numbers.Create(xps(idx1));
      yd(idx1) := yd(idx1) + factor*acc;
    else
      Forward_Backward_Cross(idx,x,fwd,bck,crs);
      Multiply_Factor(xps,fac,x,cff,pwt,acc);
      yd(0) := yd(0) + acc*fwd(idx'last-1); 
      if idx'last = 2 then
        idx2 := idx(2);
        wrk := acc*x(idx1);
        factor := Octo_Double_Numbers.Create(xps(idx2));
        wrk := factor*wrk; yd(idx2) := yd(idx2) + wrk;
        wrk := acc*x(idx2);
        factor := Octo_Double_Numbers.Create(xps(idx1));
        wrk := factor*wrk; yd(idx1) := yd(idx1) + wrk;
      else -- idk'last > 2
        wrk := acc*bck(idx'last-2);
        factor := Octo_Double_Numbers.Create(xps(idx1));
        wrk := factor*wrk;
        yd(idx1) := yd(idx1) + wrk;
        for j in idx'first+1..idx'last-1 loop
          wrk := acc*crs(j-1); idx2 := idx(j);
          factor := Octo_Double_Numbers.Create(xps(idx2));
          wrk := factor*wrk;
          yd(idx2) := yd(idx2) + wrk;
        end loop;
        wrk := acc*fwd(idx'last-2); idx2 := idx(idx'last);
        factor := Octo_Double_Numbers.Create(xps(idx2));
        wrk := factor*wrk;
        yd(idx2) := yd(idx2) + wrk;
      end if;
    end if;
  end Speel;

  procedure Speel ( xps,idx,fac : in Standard_Integer_VecVecs.VecVec;
                    cff : in OctoDobl_Complex_Vectors.Vector;
                    cst : in OctoDobl_Complex_Numbers.Complex_Number;
                    x,yd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    fwd : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    bck : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    crs : in OctoDobl_Complex_Vectors.Link_to_Vector;
                    pwt : in OctoDobl_Complex_VecVecs.VecVec;
                    h : out OctoDobl_Complex_Matrices.Matrix ) is

    use OctoDobl_Complex_Numbers,Standard_Integer_Vectors;

    dim : constant integer32 := x'last;
    idx1 : integer32;
    kcff : Complex_Number;

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        h(i,j) := Create(integer(0));
      end loop;
    end loop;
    yd(0) := cst;
    for k in xps'range loop
      declare
        xpk : constant Standard_Integer_Vectors.Vector := xps(k).all;
        idk : constant Standard_Integer_Vectors.Vector := idx(k).all;
        fck : constant Standard_Integer_Vectors.Link_to_Vector := fac(k);
      begin
        kcff := cff(k);
        if fck = null then
          if idk'last = 1 then
            idx1 := idk(1);
            yd(0) := yd(0) + kcff*x(idx1);
            yd(idx1) := yd(idx1) + kcff;
          else
            Indexed_Speel(idk,kcff,x,yd,fwd,bck,crs,h);
          end if;
        else
          Speel(xps(k),fck,idk,kcff,x,yd,fwd,bck,crs,pwt);
          if idk'last = 1 then -- one variable special case
            OctoDobl_Hessian_Updaters.Speel1(h,kcff,xpk,idk,fck.all,x.all,pwt);
          elsif idk'last = 2 then -- two variable special case
            OctoDobl_Hessian_Updaters.Speel2(h,kcff,xpk,idk,fck.all,x.all,pwt);
          elsif idk'last = 3 then -- three variable special case
            OctoDobl_Hessian_Updaters.Speel3(h,kcff,xpk,idk,fck.all,x.all,pwt);
          elsif idk'last = 4 then -- four variable special case
            OctoDobl_Hessian_Updaters.Speel4(h,kcff,xpk,idk,fck.all,x.all,pwt);
          else
            OctoDobl_Hessian_Updaters.SpeelN
              (h,kcff,xpk,idk,fck.all,x.all,fwd,bck,pwt);
          end if;
        end if;
      end;
    end loop;
    for i in 2..dim loop
      for j in 1..(i-1) loop
        h(i,j) := h(j,i);
      end loop;
    end loop;
  end Speel;

-- AUXILIARY PROCEDURES :

  procedure Forward ( x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                      f : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    use OctoDobl_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
  end Forward;

  procedure Forward_Backward
              ( x,f,b : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    use OctoDobl_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      b(k) := b(k-1)*x(x'last-k);
    end loop;
  end Forward_Backward;

  procedure Fused_Forward_Backward
              ( x,f,b : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    use OctoDobl_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    b(b'first) := x(x'last)*x(x'last-1);
    for k in 2..x'last-2 loop
      f(k) := f(k-1)*x(k+1);
      b(k) := b(k-1)*x(x'last-k);
    end loop;
    if f'last > 1
     then f(f'last) := f(f'last-1)*x(x'last);
    end if;
  end Fused_Forward_Backward;

  procedure Forward_Backward_Cross
              ( x,f,b,c : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    use OctoDobl_Complex_Numbers;

  begin
    f(f'first) := x(x'first)*x(x'first+1);
    for k in 2..x'last-1 loop
      f(k) := f(k-1)*x(k+1);
    end loop;
    if x'last > 2 then
      b(b'first) := x(x'last)*x(x'last-1);
      for k in 2..x'last-2 loop
        b(k) := b(k-1)*x(x'last-k);
      end loop;
      if x'last = 3 then
        c(1) := x(1)*x(3);
      else
        c(1) := x(1)*b(x'last-3);
        for k in 2..x'last-3 loop
          c(k) := f(k-1)*b(x'last-2-k);
        end loop;
        c(x'last-2) := x(x'last)*f(x'last-3);
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Forward_Backward_Cross
              ( idx : in Standard_Integer_Vectors.Vector;
                x,f,b,c : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    use OctoDobl_Complex_Numbers;

  begin
    f(1) := x(idx(1))*x(idx(2));
    for k in 2..idx'last-1 loop
      f(k) := f(k-1)*x(idx(k+1));
    end loop;
    if idx'last > 2 then
      b(b'first) := x(idx(idx'last))*x(idx(idx'last-1));
      for k in 2..idx'last-2 loop
        b(k) := b(k-1)*x(idx(idx'last-k));
      end loop;
      if idx'last = 3 then
        c(1) := x(idx(1))*x(idx(3));
      else
        c(1) := x(idx(1))*b(idx'last-3);
        for k in 2..idx'last-3 loop
          c(k) := f(k-1)*b(idx'last-2-k);
        end loop;
        c(idx'last-2) := x(idx(idx'last))*f(idx'last-3);
      end if;
    end if;
  end Forward_Backward_Cross;

  procedure Fused_Forward_Backward_Cross
              ( x,f,b,c : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    use OctoDobl_Complex_Numbers;

    firstend,lastend,plusidx,minidx : integer32;

  begin
    if x'last >= 8 then
      if x'last mod 2 = 0 then
        lastend := x'last-4;
        firstend := lastend/2;
        f(f'first) := x(x'first)*x(x'first+1);
        b(b'first) := x(x'last)*x(x'last-1);
        for k in 2..firstend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
        end loop;
        minidx := firstend+1; plusidx := minidx+1;
        for k in firstend+1..lastend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
          c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          plusidx := plusidx + 1;
          c(minidx) := f(minidx-1)*b(x'last-2-minidx);
          minidx := minidx - 1;
        end loop;
      else
        lastend := x'last-4;
        firstend := (x'last-3)/2;
        f(f'first) := x(x'first)*x(x'first+1);
        b(b'first) := x(x'last)*x(x'last-1);
        for k in 2..firstend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
        end loop;
        plusidx := firstend+1;
        c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
        minidx := plusidx;
        for k in firstend+1..lastend loop
          f(k) := f(k-1)*x(k+1);
          b(k) := b(k-1)*x(x'last-k);
          plusidx := plusidx + 1;
          c(plusidx) := f(plusidx-1)*b(x'last-2-plusidx);
          minidx := minidx - 1;
          c(minidx) := f(minidx-1)*b(x'last-2-minidx);
        end loop;
      end if;
      plusidx := lastend+1;
      f(plusidx) := f(plusidx-1)*x(plusidx+1);
      b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
      plusidx := plusidx+1;
      f(plusidx) := f(plusidx-1)*x(plusidx+1);
      b(plusidx) := b(plusidx-1)*x(x'last-plusidx);
      f(f'last) := f(f'last-1)*x(x'last);
      c(1) := x(1)*b(x'last-3);
      c(x'last-2) := x(x'last)*f(x'last-3);
    else
      f(f'first) := x(x'first)*x(x'first+1);
      b(b'first) := x(x'last)*x(x'last-1);
      for k in 2..x'last-2 loop
        f(k) := f(k-1)*x(k+1);
        b(k) := b(k-1)*x(x'last-k);
      end loop;
      if f'last > 1
       then f(f'last) := f(f'last-1)*x(x'last);
      end if;
      if x'last > 3 then
        c(1) := x(1)*b(x'last-3);
        for k in 2..x'last-3 loop
          c(k) := f(k-1)*b(x'last-2-k);
        end loop;
        c(x'last-2) := x(x'last)*f(x'last-3);
      elsif x'last = 3 then
        c(1) := x(1)*x(3);
      end if;
    end if;
  end Fused_Forward_Backward_Cross;

  function Allocate
             ( mxe : Standard_Integer_Vectors.Vector )
             return OctoDobl_Complex_VecVecs.VecVec is

    res : OctoDobl_Complex_VecVecs.VecVec(mxe'range);
    zero : constant OctoDobl_Complex_Numbers.Complex_Number
         := OctoDobl_Complex_Numbers.Create(integer(0));

  begin
    for k in mxe'range loop
      if mxe(k) > 1 then
        res(k) := new OctoDobl_Complex_Vectors.Vector'(1..mxe(k)-1 => zero);
      end if;
    end loop;
    return res;
  end Allocate;

  procedure Power_Table
              ( mxe : in Standard_Integer_Vectors.Vector;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                pwt : in OctoDobl_Complex_VecVecs.VecVec ) is

    lnk : OctoDobl_Complex_Vectors.Link_to_Vector;

    use OctoDobl_Complex_Numbers;

  begin
    for k in x'range loop
      if mxe(k) > 1 then
        lnk := pwt(k);
        lnk(1) := x(k)*x(k);
        for i in 2..mxe(k)-1 loop
          lnk(i) := lnk(i-1)*x(k);
        end loop;
      end if;
    end loop;
  end Power_Table;

  procedure Multiply_Factor
              ( xps,fac : in Standard_Integer_Vectors.Link_to_Vector;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                cff : in OctoDobl_Complex_Numbers.Complex_Number;
                pwt : in OctoDobl_Complex_VecVecs.VecVec;
                res : out OctoDobl_Complex_Numbers.Complex_Number ) is

    pwx : OctoDobl_Complex_Vectors.Link_to_Vector;
    idx,powidx : integer32;

    use OctoDobl_Complex_Numbers;

  begin
    idx := fac(fac'first); powidx := xps(idx);
    if powidx = 2 then
      res := cff*x(idx);
    else
      pwx := pwt(idx);
      res := cff*pwx(powidx-2);
    end if;
    for k in fac'first+1..fac'last loop
      idx := fac(k); powidx := xps(idx);
      if powidx = 2 then
        res := res*x(idx);
      else
        pwx := pwt(idx);
        res := res*pwx(powidx-2);
      end if;
    end loop;
  end Multiply_Factor;

-- DESTRUCTORS :

  procedure Clear ( c : in out Circuit ) is
  begin
    c.dim := -1;
    c.pdg := -1;
    OctoDobl_Complex_Vectors.Clear(c.fwd);
    OctoDobl_Complex_Vectors.Clear(c.bck);
    OctoDobl_Complex_Vectors.Clear(c.crs);
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
  begin
    Clear(s.crc);
    OctoDobl_Complex_Vectors.Clear(s.yd);
    OctoDobl_Complex_VecVecs.Clear(s.pwt);
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

end OctoDobl_Complex_Circuits;
