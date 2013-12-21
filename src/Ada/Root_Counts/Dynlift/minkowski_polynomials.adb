with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Cayley_Embedding;                   use Cayley_Embedding;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;

package body Minkowski_Polynomials is

  function Minkowski_Polynomial ( n,r : natural32 ) return Poly is

    res : Poly := Null_Poly;
    acc : Degrees
        := new Standard_Natural_Vectors.Vector'(1..integer32(r) => 0);

    procedure Generate_Monomials
                ( k,sum : in natural32; deg : in out Degrees ) is

    -- DESCRIPTION :
    --   Generates all exponent vectors whose sum equals n.

      t : Term;

    begin
      if k = r then
        t.cf := Create(1.0);
        deg(integer32(r)) := n-sum;
        t.dg := deg;
        Add(res,t);
      else
        for i in 0..(n-sum) loop
          deg(integer32(k)) := i;
          Generate_Monomials(k+1,sum+i,deg);
        end loop;
      end if;
    end Generate_Monomials;

  begin
    Generate_Monomials(1,0,acc);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(acc));
    return res;
  end Minkowski_Polynomial;

  function Convert ( dg : Degrees ) return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the degrees vector to a vector with integer numbers.

    res : Standard_Integer_Vectors.Vector(dg'range);

  begin
    for i in res'range loop
      res(i) := integer32(dg(i));
    end loop;
    return res;
  end Convert;

  procedure Minkowski_Polynomial
                  ( p : in out Poly; t : in Triangulation; n : in natural32;
                    mix : in Vector; mixsub : out Mixed_Subdivision ) is

    procedure Coefficient_Volume
                  ( submix : in Vector; sub : in Mixed_Subdivision;
                    vol : out natural32 ) is
    begin
      vol := Mixed_Volume(integer32(n),submix,sub);
    end Coefficient_Volume;
    procedure Coefficient_Volumes is
      new Minkowski_Polynomial_Subdivisions(Coefficient_Volume);

  begin
    Coefficient_Volumes(p,t,n,mix,mixsub);
  end Minkowski_Polynomial;

  procedure Minkowski_Polynomial_Subdivisions
                 ( p : in out Poly; t : in Triangulation; n : in natural32;
                   mix : in Vector; mixsub : out Mixed_Subdivision ) is

    procedure Coefficient_Volume ( tt : in out Term; cont : out boolean ) is

      wrkmix : Vector(mix'range) := Convert(tt.dg);
      wrksub : Mixed_Subdivision := Extract_Mixed_Cells(integer32(n),wrkmix,t);
      vol : natural32;

    begin
      Deflate(integer32(n),wrksub);
      Process(wrkmix,wrksub,vol);
      tt.cf := Create(double_float(vol));
      if wrkmix = mix
       then mixsub := wrksub;
       else Deep_Clear(wrksub);
      end if;
      cont := true;
    end Coefficient_Volume;
    procedure Coefficient_Volumes is new Changing_Iterator(Coefficient_Volume);

  begin
    Coefficient_Volumes(p);
  end Minkowski_Polynomial_Subdivisions;

end Minkowski_Polynomials;
