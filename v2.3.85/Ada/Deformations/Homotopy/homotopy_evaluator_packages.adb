with text_io;                            use text_io;
with Standard_Evaluator_Packages;        use Standard_Evaluator_Packages;

package body Homotopy_Evaluator_Packages is

  procedure Create_Homotopy_Constants ( file : in file_type ) is

  -- DESCRIPTION :
  --   Writes the code to initialize the homotopy constants.

  begin
    put_line(file,
      "  procedure Homotopy_Constants ( a : in Complex_Number; "
                                     & "k : in positive ) is");
    put_line(file,"  begin");
    put_line(file,"    aa := a;");
    put_line(file,"    kk := k;");
    put_line(file,"  end Homotopy_Constants;");
  end Create_Homotopy_Constants;

  procedure Create_Inline_Homotopy_Evaluator ( file : in file_type ) is

  -- DESCRIPTION :
  --   Writes the code to evaluate the homotopy.

  begin
    put_line(file,
      "  function Eval_Homotopy ( x : Vector; t : Complex_Number ) "
           & "return Vector is");
    new_line(file);
    put_line(file,"    y : Vector(x'range);                                  ");
    put_line(file,"    eval_target : Vector(x'range) := Eval_Target_Sys(x);  ");
    put_line(file,"    eval_astart : Vector(x'range) := aa*Eval_Start_Sys(x);");
    put_line(file,"    tpk : constant Complex_Number := t**kk;               ");
    put_line(file,"    mtk : constant Complex_Number := (Create(1.0)-t)**kk; ");
    new_line(file);
    put_line(file,"  begin");
    put_line(file,"    for i in y'range loop");
    put_line(file,"      y(i) := mtk*eval_astart(i) + tpk*eval_target(i);");
    put_line(file,"    end loop;");
    put_line(file,"    return y;");
    put_line(file,"  end Eval_Homotopy;");
  end Create_Inline_Homotopy_Evaluator;

  procedure Create_Inline_Homotopy_Differentiator1 ( file : in file_type ) is

  -- DESCRIPTION :
  --   Writes the code to differentiate the homotopy w.r.t. the variables.

  begin
    put_line(file,
      "  function Diff_Homotopy ( x : Vector; t : Complex_Number ) "
           & "return Matrix is");
    new_line(file);
    put_line(file,"    y : Matrix(x'range,x'range);                          ");
    put_line(file,"    eval_target : Matrix(x'range,x'range)"
                               & " := Eval_Target_Jaco(x);  ");
    put_line(file,"    eval_astart : Matrix(x'range,x'range)"
                               & " := Eval_Start_Jaco(x);");
    put_line(file,"    tpk : constant Complex_Number := t**kk;               ");
    put_line(file,"    mtk : constant Complex_Number"
                       & " := aa*(Create(1.0)-t)**kk; ");
    new_line(file);
    put_line(file,"  begin");
    put_line(file,"    for i in y'range(1) loop");
    put_line(file,"      for j in y'range(2) loop");
    put_line(file,"        y(i,j) := mtk*eval_astart(i,j) "
                                & "+ tpk*eval_target(i,j);");
    put_line(file,"      end loop;");
    put_line(file,"    end loop;");
    put_line(file,"    return y;");
    put_line(file,"  end Diff_Homotopy;");
  end Create_Inline_Homotopy_Differentiator1;

  procedure Create_Inline_Homotopy_Differentiator2 ( file : in file_type ) is

  -- DESCRIPTION :
  --   Writes the code to differentiate the homotopy w.r.t. t.

  begin
    put_line(file,
      "  function Diff_Homotopy ( x : Vector; t : Complex_Number ) "
           & "return Vector is");
    new_line(file);
    put_line(file,"    y : Vector(x'range);");
    new_line(file);
    put_line(file,"  begin");
    put_line(file,"    return y;");
    put_line(file,"  end Diff_Homotopy;");
  end Create_Inline_Homotopy_Differentiator2;

  procedure Create_Package_Specification
                ( file : in file_type; packname : in String ) is

  -- DESCRIPTION :
  --   Writes the specification of the homotopy evaluator package.

  begin
    put_line(file,"with Standard_Complex_Numbers;           "
                 & "use Standard_Complex_Numbers;");
    put_line(file,"with Standard_Complex_Vectors;           "
                 & "use Standard_Complex_Vectors;");
    put_line(file,"with Standard_Complex_Matrices;          "
                 & "use Standard_Complex_Matrices;");
    new_line(file);
    put_line(file,"package " & packname & " is");
    new_line(file);
    put_line(file,
      "  procedure Homotopy_Constants ( a : in Complex_Number; "
                                     & "k : in positive );");
    new_line(file);
    put_line(file,
      "  function Eval_Homotopy ( x : Vector; t : Complex_Number ) "
           & "return Vector;");
    put_line(file,
      "  function Diff_Homotopy ( x : Vector; t : Complex_Number ) "
           & "return Matrix;");
    put_line(file,
      "  function Diff_Homotopy ( x : Vector; t : Complex_Number ) "
           & "return Vector;");
    new_line(file);
    put_line(file,"end " & packname & ";");
  end Create_Package_Specification;

  procedure Create_Package_Implementation
                 ( file : in file_type; packname : in String;
                   p,q : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Writes the implementation of the homotopy evaluator package.

  begin
    put_line(file,"with Standard_Floating_Numbers;          "
                 & "use Standard_Floating_Numbers;");
    new_line(file);
    put_line(file,"package body " & packname & " is");
    new_line(file);
    put_line(file,"  aa : Complex_Number;");
    put_line(file,"  kk : positive;");
    new_line(file);
    Create_Inline_System_Evaluator(file,"Eval_Target_Sys",p);
    new_line(file);
    Create_Inline_System_Evaluator(file,"Eval_Start_Sys",q);
    new_line(file);
    Create_Inline_Jacobian_Evaluator(file,"Eval_Target_Jaco",p);
    new_line(file);
    Create_Inline_Jacobian_Evaluator(file,"Eval_Start_Jaco",q);
    new_line(file);
    Create_Homotopy_Constants(file);
    new_line(file);
    Create_Inline_Homotopy_Evaluator(file);
    new_line(file);
    Create_Inline_Homotopy_Differentiator1(file);
    new_line(file);
    Create_Inline_Homotopy_Differentiator2(file);
    new_line(file);
    put_line(file,"end " & packname & ";");
  end Create_Package_Implementation;

  procedure Create ( packname : in String; p,q : in Poly_Sys ) is

    specfile,bodyfile : file_type;

  begin
    Replace_Symbols;
    Create(specfile,out_file,packname & ".ads");
    Create_Package_Specification(specfile,packname);
    Close(specfile);
    Create(bodyfile,out_file,packname & ".adb");
    Create_Package_Implementation(bodyfile,packname,p,q);
    Close(bodyfile);
  end Create;

  procedure Create ( p,q : in Poly_Sys ) is

    packname : String := Read_Package_Name;

  begin
    Create(packname,p,q);
  end Create;

end Homotopy_Evaluator_Packages;
