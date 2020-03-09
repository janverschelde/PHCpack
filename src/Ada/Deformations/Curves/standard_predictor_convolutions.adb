with unchecked_deallocation;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Speelpenning_Convolutions;
with Newton_Convolutions;

package body Standard_Predictor_Convolutions is

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) return LU_Predictor is

    res : LU_Predictor(sol'last,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(0.0);
    knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new Standard_Complex_Vectors.Vector'(1..neq => zero);
    for k in sol'range loop
      res.numcff(k) := new Standard_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new Standard_Complex_Vectors.Vector'(kden);
    end loop;
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return SVD_Predictor is

    dim : constant integer32 := sol'last;
    res : SVD_Predictor(neq,dim,dim+1,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(0.0);
    knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

    use Standard_Speelpenning_Convolutions;

  begin
    res.sol := Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new Standard_Complex_Vectors.Vector'(1..neq => zero);
    res.ewrk := new Standard_Complex_Vectors.Vector'(1..dim => zero);
    res.dx := Allocate_Coefficients(dim,deg);
    res.xd := Linearized_Allocation(dim,deg);
    for k in sol'range loop
      res.numcff(k) := new Standard_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new Standard_Complex_Vectors.Vector'(kden);
    end loop;
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_LU_Predictor is

    prd : constant LU_Predictor(sol'last,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_LU_Predictor := new LU_Predictor'(prd);

  begin
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_SVD_Predictor is

    dim : constant integer32 := sol'last;
    prd : constant SVD_Predictor(neq,dim,dim+1,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_SVD_Predictor := new SVD_Predictor'(prd);

  begin
    return res;
  end Create;

  procedure Clear ( p : in out Link_to_LU_Predictor ) is

    procedure free is
      new unchecked_deallocation(LU_Predictor,Link_to_LU_Predictor);

  begin
    if p /= null then
      Standard_Complex_VecVecs.Clear(p.sol);
      Standard_Complex_Vectors.Clear(p.wrk);
      Standard_Complex_VecVecs.Clear(p.numcff);
      Standard_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Link_to_SVD_Predictor ) is

    procedure free is
      new unchecked_deallocation(SVD_Predictor,Link_to_SVD_Predictor);

  begin
    if p /= null then
      Standard_Complex_VecVecs.Clear(p.sol);
      Standard_Complex_Vectors.Clear(p.wrk);
      Standard_Complex_Vectors.Clear(p.ewrk);
      Standard_Complex_VecVecs.Clear(p.dx);
      Standard_Complex_VecVecs.Clear(p.xd);
      Standard_Complex_VecVecs.Clear(p.numcff);
      Standard_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

end Standard_Predictor_Convolutions;
