-- with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;

package body Standard_Inlined_Singular_Values is

  function Min0 ( a,b : integer32 ) return integer32 is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Min0;

  procedure SVD ( xrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  xiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  n,p : in integer32;
                  s : in Standard_Floating_Vectors.Link_to_Vector;
                  e : in Standard_Floating_Vectors.Link_to_Vector;
                  urv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  uiv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  vrv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  viv : in Standard_Floating_VecVecs.Link_to_VecVec;
                  job : in integer32; info : out integer32;
                  rwrk : in Standard_Floating_Vectors.Link_to_Vector;
                  iwrk : in Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    null;
  end SVD;

  procedure SVD ( x : in out Standard_Complex_Matrices.Matrix;
                  n,p : in integer32;
                  s,e : out Standard_Complex_Vectors.Vector;
                  u : out Standard_Complex_Matrices.Matrix;
                  v : out Standard_Complex_Matrices.Matrix;
                  job : in integer32; info : out integer32 ) is

    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    urv,uiv : Standard_Floating_VecVecs.Link_to_VecVec;
    vrv,viv : Standard_Floating_VecVecs.Link_to_VecVec;
    mm : constant integer32 := Min0(n+1,p);
    sf,ef,wr,wi : Standard_Floating_Vectors.Link_to_Vector;

  begin
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    urv := Standard_Vector_Splitters.Allocate(n,n,1,1);
    uiv := Standard_Vector_Splitters.Allocate(n,n,1,1);
    vrv := Standard_Vector_Splitters.Allocate(p,p,1,1);
    viv := Standard_Vector_Splitters.Allocate(p,p,1,1);
    sf := new Standard_Floating_Vectors.Vector'(1..mm => 0.0);
    ef := new Standard_Floating_Vectors.Vector'(1..p => 0.0);
    wr := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    wi := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    Standard_Matrix_Splitters.Complex_Parts(x,xrv,xiv);
    SVD(xrv,xiv,n,p,sf,ef,urv,uiv,vrv,viv,job,info,wr,wi);
    for k in s'range loop
      s(k) := Standard_Complex_Numbers.Create(sf(k),0.0);
    end loop;
    for k in e'range loop
      e(k) := Standard_Complex_Numbers.Create(ef(k),0.0);
    end loop;
    Standard_Matrix_Splitters.Complex_Merge(xrv,xiv,x);
    Standard_Matrix_Splitters.Complex_Merge(urv,uiv,u);
    Standard_Matrix_Splitters.Complex_Merge(vrv,viv,v);
    Standard_Floating_Vectors.Clear(sf);
    Standard_Floating_Vectors.Clear(ef);
    Standard_Floating_Vectors.Clear(wr);
    Standard_Floating_Vectors.Clear(wi);
    Standard_Floating_VecVecs.Deep_Clear(xrv);
    Standard_Floating_VecVecs.Deep_Clear(xiv);
    Standard_Floating_VecVecs.Deep_Clear(urv);
    Standard_Floating_VecVecs.Deep_Clear(uiv);
    Standard_Floating_VecVecs.Deep_Clear(vrv);
    Standard_Floating_VecVecs.Deep_Clear(viv);
  end SVD;

end Standard_Inlined_Singular_Values;
