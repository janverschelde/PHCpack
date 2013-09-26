with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Polyhedral_Coefficient_Predictors is

  procedure Predictor ( s,ds : in out double_float ) is
  begin
    s := s + ds;
    if s > 0.0
     then ds := s - ds;  -- old value for s
          s := 0.0;
          ds := abs(ds); -- do not overshoot again
    end if;
  end Predictor;

  procedure Step_Control
              ( fail : in boolean; s,ds : in out double_float;
                max_ds : in double_float; cnt : in out natural32;
                x,backup : in out Standard_Complex_Vectors.Vector ) is

    expansion : constant double_float := 1.5;
    reduction : constant double_float := 2.0;
    threshold : constant natural32 := 5;

  begin
    if fail then
      s := s - ds;
      ds := ds/reduction;
      x := backup;
      cnt := 0;
    else
      if (ds < max_ds) and (cnt > threshold)
       then ds := expansion*ds;
      end if;
      backup := x;
      cnt := cnt + 1;
    end if;
  end Step_Control;

  procedure Step_Control
              ( fail : in boolean; s,ds : in out double_float;
                max_ds : in double_float; cnt : in out natural32;
                x,backup,px : in out Standard_Complex_Vectors.Vector ) is

    expansion : constant double_float := 1.5;
    reduction : constant double_float := 2.0;
    threshold : constant natural32 := 5;

  begin
    if fail then
      s := s - ds;
      ds := ds/reduction;
      x := backup;
      cnt := 0;
    else
      if (ds < max_ds) and (cnt > threshold)
       then ds := expansion*ds;
      end if;
      px := backup;
      backup := x;
      cnt := cnt + 1;
    end if;
  end Step_Control;

  procedure Secant_Predictor 
              ( x : in out Standard_Complex_Vectors.Vector;
                px : in Standard_Complex_Vectors.Vector;
                h : in double_float ) is
  begin
    if h > 0.0 then
      for i in x'range loop
        x(i) := x(i) + h*(x(i) - px(i));
      end loop;
    end if;
  end Secant_Predictor;

end Polyhedral_Coefficient_Predictors;
