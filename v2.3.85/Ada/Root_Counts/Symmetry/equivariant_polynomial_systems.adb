with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Permutations,Permute_Operations;    use Permutations,Permute_Operations;

package body Equivariant_Polynomial_Systems is

  procedure Act ( v : in List_of_Permutations; s : in Poly_Sys;
                  w : in out List_of_Permutations;
                  fail,inva,equi : out boolean ) is
 
    min_s : Poly_Sys(s'range);
    last_w : List_of_Permutations;
    wrkinva,wrkequi : boolean;

    procedure Process ( p : in Permutation; cont : in out boolean ) is

      ps : Poly; -- the permuted polynomial
      pp : Permutation(p'range);

    begin
      for i in s'range loop
        ps := p*s(i);
        pp(i) := p'last+1;
        for j in s'range loop
          if Equal(ps,s(j))
           then pp(i) := j;
           elsif Equal(ps,min_s(j))
               then pp(i) := -j;
          end if;
        end loop;
        if pp(i) = p'last+1
         then fail := true;
        end if;
        Clear(ps);
      end loop;
      if wrkinva then
        for j in pp'range loop
          wrkinva := (pp(j) = j);
          exit when not wrkinva;
        end loop;
      end if;
      if wrkequi
       then wrkequi := Equal(pp,p);
      end if;
      Append(w,last_w,pp);
      cont := true;
    end Process;
    procedure Act_of_Permutations is new Iterator(Process);

  begin
    min_s := -s;
    fail := false;
    wrkinva := true; wrkequi := true;
    Act_of_Permutations(v);
    inva := wrkinva; equi := wrkequi;
    Clear(min_s);
  end Act;

  function Symmetric ( s : Poly_Sys; v,w : List_of_Permutations )
                     return boolean is

    lw,plw,pw : List_of_Permutations;
    fail,inva,equi : boolean;

  begin
    Act(v,s,lw,fail,inva,equi);
    pw := w;
    plw := lw;
    while not Is_Null(plw) loop
      if not Equal(Permutation(Head_Of(plw).all),Permutation(Head_Of(pw).all))
       then Clear(lw);
	    return false;
       else plw := Tail_Of(plw);
	    pw := Tail_Of(pw);
      end if;
    end loop;
    Clear(lw);
    return true;
  end Symmetric;

end Equivariant_Polynomial_Systems;
