with Bracket_Monomials;                use Bracket_Monomials;
with DoblDobl_Complex_Numbers_cv;      use DoblDobl_Complex_Numbers_cv; 
with QuadDobl_Complex_Numbers_cv;      use QuadDobl_Complex_Numbers_cv; 

package body Bracket_Polynomial_Convertors is

  function Convert
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial )
             return  DoblDobl_Bracket_Polynomials.Bracket_Polynomial is

    res : DoblDobl_Bracket_Polynomials.Bracket_Polynomial
        := DoblDobl_Bracket_Polynomials.Null_Bracket_Poly;

    procedure Convert_Term
                ( t : in Standard_Bracket_Polynomials.Bracket_Term;
                  c : out boolean ) is

      dd_t : DoblDobl_Bracket_Polynomials.Bracket_Term;

    begin
      dd_t.coeff := Standard_to_DoblDobl_Complex(t.coeff);
      Copy_Append(t.monom,dd_t.monom);
      DoblDobl_Bracket_Polynomials.Frontal_Construct(res,dd_t);
      c := true;
    end Convert_Term;
    procedure Convert_Terms is new
      Standard_Bracket_Polynomials.Enumerate_Terms(Convert_Term);

  begin
    Convert_Terms(p);
    return res;
  end Convert;

  function Convert
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial )
             return QuadDobl_Bracket_Polynomials.Bracket_Polynomial is

    res : QuadDobl_Bracket_Polynomials.Bracket_Polynomial
        := QuadDobl_Bracket_Polynomials.Null_Bracket_Poly;

    procedure Convert_Term
                ( t : in Standard_Bracket_Polynomials.Bracket_Term;
                  c : out boolean ) is

      dd_t : QuadDobl_Bracket_Polynomials.Bracket_Term;

    begin
      dd_t.coeff := Standard_to_QuadDobl_Complex(t.coeff);
      Copy_Append(t.monom,dd_t.monom);
      QuadDobl_Bracket_Polynomials.Frontal_Construct(res,dd_t);
      c := true;
    end Convert_Term;
    procedure Convert_Terms is new
      Standard_Bracket_Polynomials.Enumerate_Terms(Convert_Term);

  begin
    Convert_Terms(p);
    return res;
  end Convert;

end Bracket_Polynomial_Convertors;
