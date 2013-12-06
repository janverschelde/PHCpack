with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Laur_Poly_Convertors;      use Standard_Laur_Poly_Convertors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Standard_Simpomial_Solvers;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;
with Symmetric_BKK_Bound_Solvers;        use Symmetric_BKK_Bound_Solvers;
with Orbits_of_Solutions;                use Orbits_of_Solutions;

package body Symmetric_Polyhedral_Continuation is

  function Select_Subsystem ( p : Laur_Sys; mix : Vector; mic : Mixed_Cell )
                            return Laur_Sys is

    res : Laur_Sys(p'range);
    cnt : integer32 := 0;

  begin
    for k in mix'range loop
      for l in 1..mix(k) loop
        cnt := cnt + 1;
        res(cnt) := Select_Terms(p(cnt),mic.pts(k));
      end loop;
    end loop;
    return res;
  end Select_Subsystem;

  function Symmetric_Mixed_Solve
                ( file : file_type; grp : List_of_Permutations; sign : boolean;
                  p : Laur_Sys; mixsub : Mixed_Subdivision;
                  n : integer32; mix : Vector ) return Solution_List is

    sols,sols_last : Solution_List;
    cnt : integer32 := 0;
    tmp : Mixed_Subdivision := mixsub;
    tol_zero : constant double_float := 1.0E-12;

    procedure Solve_Subsystem ( mic : in Mixed_Cell ) is
  
      q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
      sq : Laur_Sys(q'range);
      qsols : Solution_List;
      fail,zero_y : boolean;
      eps : constant double_float := 10.0**(-10);

    begin
      new_line(file);
      put(file,"*** CONSIDERING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Reduce(n+1,q); sq := Shift(q);
      declare
        pq : Poly_Sys(q'range) := Laurent_to_Polynomial_System(sq);
      begin
        Standard_Simpomial_Solvers.Solve(sq,tol_zero,qsols,fail,zero_y);
        if not fail then
          put_line(file,"It is a simplex system.");
        else
          put_line(file,"No simplex system.");
          if mic.sub = null then
            put_line(file,"Calling the black box solver.");
            qsols := Symmetric_BKK_Solve(file,pq,grp,sign);
          else
            put_line(file,"Using the refinement of the cell.");
            declare
              sup : Array_of_Lists(q'range);
              cnt : integer32 := sup'first;
              lif : Array_of_Lists(mix'range);
              lifq : Laur_Sys(q'range);
            begin
              for i in mic.pts'range loop
                sup(cnt) := Reduce(mic.pts(i),q'last+1);
                for j in 1..(mix(i)-1) loop
                  Copy(sup(cnt),sup(cnt+j));
                end loop;
                cnt := cnt + mix(i);
              end loop;
              lif := Induced_Lifting(n,mix,sup,mic.sub.all);
              lifq := Perform_Lifting(n,mix,lif,q);
              qsols := Symmetric_Mixed_Solve
                         (file,grp,sign,lifq,mic.sub.all,n,mix);
              Deep_Clear(sup); Deep_Clear(lif); Clear(lifq);
            end;
          end if;
          Set_Continuation_Parameter(qsols,Create(0.0));
        end if;
        put(file,Length_Of(qsols),1);
        put_line(file," solutions found.");
        if not Is_Null(qsols) then
          Analyze(grp,sign,eps,qsols);
          put(file,Length_Of(qsols),1);
          put_line(file," generating solutions found.");
          Mixed_Continuation(file,p,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
        end if;
        Clear(pq); Clear(sq);
      end;
      Clear(q); -- Shallow_Clear(qsols);
    end Solve_Subsystem;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      Solve_Subsystem(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return sols;
  end Symmetric_Mixed_Solve;

  function Symmetric_Mixed_Solve
                ( file : file_type; sign : boolean; p : Laur_Sys;
                  mixsub : Mixed_Subdivision; n : integer32;
                  mix : Vector ) return Solution_List is

    sols,sols_last : Solution_List;
    cnt : integer32 := 0;
    tmp : Mixed_Subdivision := mixsub;
    tol_zero : constant double_float := 1.0E-12;

    procedure Solve_Subsystem ( mic : in Mixed_Cell ) is
  
      q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
      sq : Laur_Sys(q'range);
      qsols,genqsols : Solution_List;
      fail,zero_y : boolean;
      eps : constant double_float := 10.0**(-10);

    begin
      new_line(file);
      put(file,"*** CONSIDERING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Reduce(n+1,q); sq := Shift(q);
      declare
        pq : Poly_Sys(q'range) := Laurent_to_Polynomial_System(sq);
      begin
        Standard_Simpomial_Solvers.Solve(sq,tol_zero,qsols,fail,zero_y);
        if not fail then
          put_line(file,"It is a simplex system.");
        else
          put_line(file,"No simplex system.");
          if mic.sub = null then
            put_line(file,"Calling the black box solver.");
            qsols := Symmetric_BKK_Solve(file,pq,sign);
          else
            put_line(file,"Using the refinement of the cell.");
            declare
              sup : Array_of_Lists(q'range);
              cnt : integer32 := sup'first;
              lif : Array_of_Lists(mix'range);
              lifq : Laur_Sys(q'range);
            begin
              for i in mic.pts'range loop
                sup(cnt) := Reduce(mic.pts(i),q'last+1);
                for j in 1..(mix(i)-1) loop
                  Copy(sup(cnt),sup(cnt+j));
                end loop;
                cnt := cnt + mix(i);
              end loop;
              lif := Induced_Lifting(n,mix,sup,mic.sub.all);
              lifq := Perform_Lifting(n,mix,lif,q);
              qsols := Symmetric_Mixed_Solve(file,sign,lifq,mic.sub.all,n,mix);
              Deep_Clear(sup); Deep_Clear(lif); Clear(lifq);
            end;
          end if;
          Set_Continuation_Parameter(qsols,Create(0.0));
        end if;
        put(file,Length_Of(qsols),1);
        put_line(file," solutions found.");
        if not Is_Null(qsols) then
          genqsols := Generating(qsols,sign,eps);
          put(file,Length_Of(genqsols),1);
          put_line(file," generating solutions found.");
          Mixed_Continuation(file,p,mic.nor.all,genqsols);
          Concat(sols,sols_last,genqsols);
        end if;
        Clear(pq); Clear(sq);
      end;
      Clear(q); -- Shallow_Clear(genqsols);
    end Solve_Subsystem;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      Solve_Subsystem(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return sols;
  end Symmetric_Mixed_Solve;

end Symmetric_Polyhedral_Continuation;
