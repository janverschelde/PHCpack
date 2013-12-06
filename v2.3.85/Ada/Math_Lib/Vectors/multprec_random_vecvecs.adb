with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Integer_Vectors;
with Multprec_Floating_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;

package body Multprec_Random_VecVecs is

  function Random_VecVec ( n,m,sz : natural32 )
                         return Multprec_Integer_VecVecs.VecVec is

    res : Multprec_Integer_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Multprec_Integer_Vectors.Vector'
                      (Random_Vector(1,integer32(n),sz));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,sz : natural32 )
                         return Multprec_Floating_VecVecs.VecVec is

    res : Multprec_Floating_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Multprec_Floating_Vectors.Vector'
                      (Random_Vector(1,integer32(n),sz));
    end loop;
    return res;
  end Random_VecVec;

  function Random_VecVec ( n,m,sz : natural32 )
                         return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(1..integer32(m));

  begin
    for i in 1..integer32(m) loop
      res(i) := new Multprec_Complex_Vectors.Vector'
                      (Random_Vector(1,integer32(n),sz));
    end loop;
    return res;
  end Random_VecVec;

end Multprec_Random_VecVecs;
