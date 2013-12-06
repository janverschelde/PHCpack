with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;

package body Standard_Random_VecVecs is

  function Random_VecVec ( n,m : natural32; low,upp : integer32 )
                         return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Standard_Integer_Vectors.Vector'
                      (Random_Vector(1,integer32(n),low,upp));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m : natural32 )
                         return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Standard_Floating_Vectors.Vector'
                      (Random_Vector(1,integer32(n)));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,g : natural32 )
                         return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Standard_Floating_Vectors.Vector'
                      (Random_Vector(1,integer32(n),g));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m : natural32 )
                         return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Standard_Complex_Vectors.Vector'
                      (Random_Vector(1,integer32(n)));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,g : natural32 )
                         return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Standard_Complex_Vectors.Vector'
                      (Random_Vector(1,integer32(n),g));
    end loop;
    return res;
  end Random_VecVec;

end Standard_Random_VecVecs;
