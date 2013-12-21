with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;

package body DoblDobl_Random_VecVecs is

  function Random_VecVec ( n,m : natural32 )
                         return Double_Double_VecVecs.VecVec is

    res : Double_Double_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Double_Double_Vectors.Vector'
                      (Random_Vector(1,integer32(n)));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,g : natural32 )
                         return Double_Double_VecVecs.VecVec is

    res : Double_Double_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Double_Double_Vectors.Vector'
                      (Random_Vector(1,integer32(n),g));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m : natural32 )
                         return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'
                      (Random_Vector(1,integer32(n)));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,g : natural32 )
                         return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new DoblDobl_Complex_Vectors.Vector'
                       (Random_Vector(1,integer32(n),g));
    end loop;
    return res;
  end Random_VecVec;

end DoblDobl_Random_VecVecs;
