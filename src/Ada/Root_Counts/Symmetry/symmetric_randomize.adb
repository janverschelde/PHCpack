with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Permutations,Permute_Operations;    use Permutations,Permute_Operations;

function Symmetric_Randomize ( p : Laur_Sys; v,w : List_of_Permutations )
                             return Laur_Sys is

  res : Laur_Sys(p'range);

  procedure Permute_and_Randomize ( t : in Term; index : integer32 ) is

    tmpv,tmpw : List_of_Permutations;

  begin
    tmpv := v;  tmpw := w;
    while not Is_Null(tmpv) loop
      declare
        permt : Term := Permutation(Head_Of(tmpv).all)*t;
        indw : constant integer32 := Head_Of(tmpw)(index);
      begin
        if Coeff(res(indw),permt.dg) = Create(0.0)
         then Add(res(indw),permt);
        end if;
        Clear(permt);
      end;
      tmpv := Tail_Of(tmpv); 
      tmpw := Tail_Of(tmpw);
    end loop;
  end Permute_and_Randomize;
  
  procedure Symmetric_Randomize_Terms
              ( index : in integer32; py : in Poly ) is

    tpy : Term;

    procedure Pick_Term ( t : in Term; cont : out boolean ) is
    begin
      if Coeff(res(index),t.dg) = Create(0.0) then
        Copy(t,tpy);
        tpy.cf := Random1;
        cont := false;
      else
        cont := true;
      end if;
    end Pick_Term;
    procedure Pick_A_Term is new Visiting_Iterator(Pick_Term);

  begin
    tpy.cf := Create(0.0);
    Pick_A_Term(py);
    if tpy.cf /= Create(0.0) then
      Permute_and_Randomize(tpy,index);
      Clear(tpy);
    end if;
  end Symmetric_Randomize_Terms;

begin
  res := (res'range => Null_Poly);
  for k in res'range loop
    while Number_of_Terms(res(k)) < Number_of_Terms(p(k)) loop
      Symmetric_Randomize_Terms(k,p(k));
    end loop;
  end loop;
  return res;
end Symmetric_Randomize;
