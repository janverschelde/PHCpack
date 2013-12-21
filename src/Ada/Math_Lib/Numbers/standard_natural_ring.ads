with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring;

package Standard_Natural_Ring is
  new Abstract_Ring(natural32,0,1,
                    Standard_Natural_Numbers.Create,
                    Standard_Natural_Numbers.Equal,
                    Standard_Natural_Numbers.Copy,
                    "+",
                    "+",
                    "-",
                    "-",
                    "*",
                    Standard_Natural_Numbers.Add,
                    Standard_Natural_Numbers.Sub,
                    Standard_Natural_Numbers.Min,
                    Standard_Natural_Numbers.Mul,
                    Standard_Natural_Numbers.Clear);

-- DESCRIPTION :
--   Defines the ring of standard natural numbers.
