with Dictionaries;

package body Linear_Programming is

  procedure Primal_Simplex 
                 ( dic : in out Matrix; eps : in double_float;
                   in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                   nit : in out natural32; unbounded : out boolean ) is

    unbound : boolean;

  begin
    while not Dictionaries.Primal_Optimal(dic,eps) loop
      Dictionaries.Primal_Update(dic,in_bas,out_bas,eps,unbound);
      nit := nit + 1;
      exit when unbound;
    end loop;
    unbounded := unbound;
  end Primal_Simplex;

  procedure Generic_Primal_Simplex
                 ( dic : in out Matrix; eps : in double_float;
                   in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                   nit : in out natural32; unbounded : out boolean ) is

    unbound : boolean;

  begin
    while not Dictionaries.Primal_Optimal(dic,eps) loop
      Report(dic,in_bas,out_bas);
      Dictionaries.Primal_Update(dic,in_bas,out_bas,eps,unbound);
      nit := nit + 1;
      exit when unbound;
    end loop;
    Report(dic,in_bas,out_bas);
    unbounded := unbound;
  end Generic_Primal_Simplex;

  procedure Dual_Simplex
                 ( dic : in out Matrix; eps : in double_float;
                   in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                   nit : in out natural32; feasible : out boolean ) is

    feasi : boolean;

  begin
    while not Dictionaries.Dual_Optimal(dic,eps) loop
      Dictionaries.Dual_Update(dic,in_bas,out_bas,eps,feasi);
      nit := nit + 1;
      exit when not feasi;
    end loop;
    feasible := feasi;
  end Dual_Simplex;

  procedure Generic_Dual_Simplex
                 ( dic : in out Matrix; eps : in double_float;
                   in_bas,out_bas : in out Standard_Integer_Vectors.Vector;
                   nit : in out natural32; feasible : out boolean ) is

    feasi : boolean;

  begin
    while not Dictionaries.Dual_Optimal(dic,eps) loop
      Report(dic,in_bas,out_bas);
      Dictionaries.Dual_Update(dic,in_bas,out_bas,eps,feasi);
      nit := nit + 1;
      exit when not feasi;
    end loop;
    Report(dic,in_bas,out_bas);
    feasible := feasi;
  end Generic_Dual_Simplex;

end Linear_Programming;
