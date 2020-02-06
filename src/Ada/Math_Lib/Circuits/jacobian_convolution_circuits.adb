with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;

package body Jacobian_Convolution_Circuits is

  function Jacobian ( c : Standard_Speelpenning_Convolutions.Circuits;
                      x : Standard_Complex_Vectors.Vector ) 
                    return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(c'range,x'range);

  begin
    for i in c'range loop
      for j in x'range loop
        res(i,j) := Standard_Speelpenning_Convolutions.Diff(c(i),x,j);
      end loop;
    end loop;
    return res;
  end Jacobian;

  function Jacobian ( c : DoblDobl_Speelpenning_Convolutions.Circuits;
                      x : DoblDobl_Complex_Vectors.Vector ) 
                    return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(c'range,x'range);

  begin
    for i in c'range loop
      for j in x'range loop
        res(i,j) := DoblDobl_Speelpenning_Convolutions.Diff(c(i),x,j);
      end loop;
    end loop;
    return res;
  end Jacobian;

  function Jacobian ( c : QuadDobl_Speelpenning_Convolutions.Circuits;
                      x : QuadDobl_Complex_Vectors.Vector ) 
                    return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(c'range,x'range);

  begin
    for i in c'range loop
      for j in x'range loop
        res(i,j) := QuadDobl_Speelpenning_Convolutions.Diff(c(i),x,j);
      end loop;
    end loop;
    return res;
  end Jacobian;

   procedure Singular_Values
               ( c : in Standard_Speelpenning_Convolutions.Circuits;
                 x : in Standard_Complex_Vectors.Vector;
                 A : out Standard_Complex_Matrices.Matrix;
                 U : out Standard_Complex_Matrices.Matrix;
                 V : out Standard_Complex_Matrices.Matrix;
                 e : out Standard_Complex_Vectors.Vector;
                 s : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 11;
    info : integer32;

  begin
    A := Jacobian(c,x);
    Standard_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
  end Singular_Values;

   procedure Singular_Values
               ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                 x : in DoblDobl_Complex_Vectors.Vector;
                 A : out DoblDobl_Complex_Matrices.Matrix;
                 U : out DoblDobl_Complex_Matrices.Matrix;
                 V : out DoblDobl_Complex_Matrices.Matrix;
                 e : out DoblDobl_Complex_Vectors.Vector;
                 s : out DoblDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 11;
    info : integer32;

  begin
    A := Jacobian(c,x);
    DoblDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
  end Singular_Values;

   procedure Singular_Values
               ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                 x : in QuadDobl_Complex_Vectors.Vector;
                 A : out QuadDobl_Complex_Matrices.Matrix;
                 U : out QuadDobl_Complex_Matrices.Matrix;
                 V : out QuadDobl_Complex_Matrices.Matrix;
                 e : out QuadDobl_Complex_Vectors.Vector;
                 s : out QuadDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 11;
    info : integer32;

  begin
    A := Jacobian(c,x);
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,job,info);
  end Singular_Values;

end Jacobian_Convolution_Circuits;
