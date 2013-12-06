with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;

package body Multprec_Residual_Evaluations is

  function Residual ( p_eval : Eval_Poly_Sys; zero : Vector )
                    return Floating_Number is

    res : Floating_Number;
    eva : Vector(p_eval'range) := Eval(p_eval,zero);

  begin
    res := Max_Norm(eva);
    Clear(eva);
    return res;
  end Residual;

  procedure Residual ( file : in file_type;
                       p_eval : in Eval_Poly_Sys; sol : in Solution ) is

    res : Floating_Number := Residual(p_eval,sol.v);

  begin
    put(file,res);
    Clear(res);
  end Residual;

  procedure Residuals ( file : in file_type;
                        p_eval : in Eval_Poly_Sys; sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      put(file,"residual "); put(file,i,1); put(file," : ");
      Residual(file,p_eval,ls.all);
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Residuals;

end Multprec_Residual_Evaluations;
