with unchecked_deallocation;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with Newton_Convolutions;

package body QuadDobl_Predictor_Convolutions is

  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32 ) return Predictor is

    res : Predictor(sol'last,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(integer(0));
    knum : constant QuadDobl_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant QuadDobl_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new QuadDobl_Complex_Vectors.Vector'(1..neq => zero);
    for k in sol'range loop
      res.numcff(k) := new QuadDobl_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new QuadDobl_Complex_Vectors.Vector'(kden);
    end loop;
    return res;
  end Create;

  function Create ( sol : QuadDobl_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_Predictor is

    prd : constant Predictor(sol'last,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_Predictor := new Predictor'(prd);

  begin
    return res;
  end Create;

  procedure Clear ( p : in out Link_to_Predictor ) is

    procedure free is new unchecked_deallocation(Predictor,Link_to_Predictor);

  begin
    if p /= null then
      QuadDobl_Complex_VecVecs.Clear(p.sol);
      QuadDobl_Complex_VecVecs.Clear(p.numcff);
      QuadDobl_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

end QuadDobl_Predictor_Convolutions;
