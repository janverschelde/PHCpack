with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;

package body QuadDobl_Random_VecVecs is

  function Random_VecVec ( n,m : natural32 )
                         return Quad_Double_VecVecs.VecVec is

    res : Quad_Double_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Quad_Double_Vectors.Vector'
                      (Random_Vector(1,integer32(n)));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,g : natural32 )
                         return Quad_Double_VecVecs.VecVec is

    res : Quad_Double_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Quad_Double_Vectors.Vector'
                      (Random_Vector(1,integer32(n),g));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m : natural32 )
                         return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'
                      (Random_Vector(1,integer32(n)));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,g : natural32 )
                         return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new QuadDobl_Complex_Vectors.Vector'
                      (Random_Vector(1,integer32(n),g));
    end loop;
    return res;
  end Random_VecVec;

end QuadDobl_Random_VecVecs;
