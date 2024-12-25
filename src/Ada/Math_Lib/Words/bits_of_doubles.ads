with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package Bits_of_Doubles is

-- DESCRIPTION :
--   Provides functions to manipulate the bits of the fraction
--   of a double precision floating point number.

  procedure expand_52bits
              ( bits : out Standard_Natural_Vectors.Vector;
                nbr : in integer64 );

  -- DESCRIPTION :
  --   Given in bits space for 52 bits, in the vector bits of range 0..51,
  --   fills the bits with the binary expansin of the number nbr.

  procedure write_52bits
              ( bits : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes 52 bits in groups of 4.

  procedure Fraction_Exponent
              ( x : in double_float;
                f : out integer64; e : out integer32 );

  -- DESCRIPTION :
  --   Given a double float number x, returns the fraction f
  --   as a 64-bit integer from the fraction of x, multiplied by 2^52,
  --   and returns the exponent e as a 32-bit integer.

  procedure write_52bits_expo ( x : in double_float );

  -- DESCRIPTION :
  --   Shows the bits sequence of the fraction of x, in groups of 4,
  --   and the decimal representation of the exponent of x.

  function Bit_Equal ( x,y : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x and y have the same sign, same exponent,
  --   and the same bits in their fraction.

  procedure write_fraction_bits ( nbr : in double_float );

  -- DESCRIPTION :
  --   Writes the bits of the fraction of nbr in binary,
  --   without the use of an auxiliary vector of natural numbers.

  function value_52bits
             ( bits : Standard_Natural_Vectors.Vector ) return integer64;

  -- DESCRIPTION :
  --   Given 52 bits in the vector bits, or range 0..51,
  --   returns its value as a 64-bit integer.

  function chop_last_bits
             ( nbr : double_float; lastbits : natural32 )
             return double_float;

  -- DESCRIPTION :
  --   Returns the number nbr with the last bits, equal as lastbits,
  --   chopped off, that is: set to zero.

  -- REQUIRED : lastbits < 52.

  procedure chop_last_bits
             ( nbr : in out double_float; lastbits : in natural32;
               headbits : out Standard_Natural_Vectors.Vector;
               tailbits : out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the last bits of the number nbr to zero.

  -- REQUIRED : lastbits < 52.

  -- ON ENTRY :
  --   nbr     a 64-bit double floating-point number;
  --   lastbits is the number of last bits to set to zero.

  -- ON RETURN :
  --   nbr     the number with the last bits removed;
  --   headbits are the leading bits which remain unchanged;
  --   tailbits are the last bits of the original number nbr.

  procedure insert_first_bits
              ( bits : in out Standard_Natural_Vectors.Vector;
                firstbits : in natural32;
                headbits : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Inserts the first bits, as many as the value of firstbits,
  --   of the headbits, into the bits.

  -- REQUIRED : firstbits < 52.

  -- ON ENTRY :
  --   bits     a sequence of 52 bits, of range 0..51;
  --   firstbits is the number of bits to be inserted;
  --   headbits contains the first bits to be inserted.

  -- ON RETURN :
  --   bits     contains the first bits of headbits,
  --            with the original bits shifted to the end.

  procedure insert_first_bits
              ( nbr : in out double_float;
                firstbits : in natural32;
                headbits : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Inserts the first bits, as many as the value of firstbits,
  --   of the headbits, into the fraction of the number nbr.

  -- REQUIRED : firstbits < 52.

  -- ON ENTRY :
  --   nbr      a 64-bit float;
  --   firstbits is the number of bits to be inserted;
  --   headbits contains the first bits to be inserted.

  -- ON RETURN :
  --   nbr      its fraction contains the first bits of headbits,
  --            with the original bits shifted to the end.

  function insert_first_bits
             ( nbr : double_float;
               firstbits : natural32;
               headbits : Standard_Natural_Vectors.Vector )
             return double_float;

  -- DESCRIPTION :
  --   Inserts the first bits, as many as the value of firstbits,
  --   of the headbits, into the fraction of the number nbr.
  --   Returns the number with the original bits of the fraction
  --   of nbr shifted to the end and the first bits inserted
  --   from headbits.

  -- REQUIRED : firstbits < 52.

  -- ON ENTRY :
  --   nbr      a 64-bit float;
  --   firstbits is the number of bits to be inserted;
  --   headbits contains the first bits to be inserted.

  procedure Mod_Split ( x : in double_float;
                        xhi,xlo : out double_float );

  -- DESCRIPTION :
  --   Splits the fraction of the double x in two equal halves,
  --   using modular arithmetic, without vectors of natural numbers.

  procedure Vec_Split ( x : in double_float;
                        xhi,xlo : out double_float );

  -- DESCRIPTION :
  --   Splits the fraction of the double x in two equal halves,
  --   using vectors of natural numbers.

  procedure Split ( x : in double_float;
                    x0,x1,x2,x3 : out double_float );

  -- DESCRIPTION :
  --   Splits the 52 bits in the fraction of x in four equal parts,
  --   returned in x0, x1, x2, x3, with x0 > x1 > x2 > x3.
  --   On return: Bit_Equal(x,x0+x1+x2+x3) is true
  --   and x - (x0 + x1 + x2 + x3) is exactly zero.

  function Different_Sign ( x,y : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if x and y have different signs, false otherwise.
  --   This function should be used to avoid underflow of x*y
  --   which causes the test x*y < 0.0 to fail.

  procedure Sign_Balance ( hi,lo : in out double_float;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   Given hi*lo < 0.0, balances the sign by redistributing
  --   the bits from hi to lo.
  --   If verbose, prints results of intermediate computations.

  procedure Sign_Balance ( x : in out double_double;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   If the high and the low part of x have different signs,
  --   then the bits of x are redistributed so the double double
  --   representation of x is sign balanced.
  --   If verbose, prints results of intermediate computations.

  procedure Sign_Balance ( hihi,lohi,hilo,lolo : in out double_float;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   Balances the doubles of a quad double so all parts have the same sign.

  procedure Sign_Balance
              ( hihihi,lohihi,hilohi,lolohi : in out double_float;
                hihilo,lohilo,hilolo,lololo : in out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Balances the doubles of an octo double so all parts have the same sign.

  procedure Sign_Balance
              ( hihihihi,lohihihi,hilohihi,lolohihi : in out double_float;
                hihilohi,lohilohi,hilolohi,lololohi : in out double_float;
                hihihilo,lohihilo,hilohilo,lolohilo : in out double_float;
                hihilolo,lohilolo,hilololo,lolololo : in out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Balances the doubles of a hexa double so all parts have the same sign.

  procedure Sign_Balance ( x : in out quad_double;
                           verbose : in boolean := true );
  procedure Sign_Balance ( x : in out octo_double;
                           verbose : in boolean := true );
  procedure Sign_Balance ( x : in out hexa_double;
                           verbose : in boolean := true );

  -- DESCRIPTION :
  --   Redistributes the bits so all parts of have the same sign.

  function Is_Sign_Balanced ( x : double_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true of both the high and low part of x have the same sign.

  function Is_Sign_Balanced ( x : quad_double ) return boolean;
  function Is_Sign_Balanced ( x : octo_double ) return boolean;
  function Is_Sign_Balanced ( x : hexa_double ) return boolean;

  -- DESCRIPTION :
  --   Returns true of all parts of x have the same sign.

end Bits_of_Doubles;
