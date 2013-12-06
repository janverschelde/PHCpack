with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Equation_by_Equation_Solvers is 

-- DESCRIPTION :
--   This package offers equation-by-equation solvers, 
--   for polynomial systems given as polynomials (prefix P_),
--   or given as polynomial functions (prefix G_).

  procedure P_Solve_Equation_by_Equation
              ( file : in file_type;
                filename : in string; report,step : in boolean;
                ne,nv,k : in integer32; p : in Poly_Sys;
                witset : in out Array_of_Solution_Lists;
                planes : in out VecMat );

  -- DESCRIPTION :
  --   Applies the equation-by-equation solver to the system p.

  -- REQUIRED : witset'range = planes'range = 1..max(ne,nv).

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   filename is used to create separate output files;
  --   report   if false, then path trackers will remain silent,
  --            otherwise diagnostics along paths are reported;
  --   step     flag to ask for stepping stones;
  --   ne       number of equations of the system p, p'range = 1..ne;
  --   nv       number of variables of the system p;
  --   k        index of the current equation, if k > 0, a witness
  --            set is expected for the first k equations of p;
  --   p        system of ne equations in nv unknowns;
  --   witset   witness sets for the first k equations of p;
  --   planes   corresponding planes defining the witness sets
  --            for the first k equations of p.

  -- ON RETURN :
  --   witset   witness sets, the i-th entry is the i-dimensional
  --            solution set of the system p;
  --   planes   planes(i) cuts out the points in witset(i).

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;
     -- returns f[k](x), the value of the k-th polynomial at x

    with function jf ( k : integer32; x : Vector ) return Vector;
     -- returns jf[k](x), the k-th row of the Jacobian matrix at x

  procedure G_Solve_Equation_by_Equation
              ( file : in file_type; filename : in string;
                report,step : in boolean;
                ne,nv,k : in integer32; d : in Standard_Natural_Vectors.Vector;
                witset : out Array_of_Solution_Lists; planes : out VecMat );

  -- DESCRIPTION :
  --   Applies the equation-by-equation solver to the system,
  --   evaluated by f and with Jacobian matrix in jf.

  -- REQUIRED : witset'range = planes'range = 1..max(ne,nv).

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   filename is used to create separate output files;
  --   report   if false, then path trackers will remain silent,
  --            otherwise diagnostics along paths are reported;
  --   step     flag to ask for stepping stones;
  --   ne       number of equations of the system;
  --   nv       number of variables of the system;
  --   k        index of the current equation, if k > 0, a witness
  --            set is expected for the first k equations of p;
  --   d        degrees of the polynomials in the system.

  -- ON RETURN :
  --   witset   witness sets, the i-th entry is the i-dimensional
  --            solution set of the system p;
  --   planes   planes(i) cuts out the points in witset(i).

  generic

    with function f ( k : integer32; x : Vector ) return Complex_Number;
     -- returns f[k](x), the value of the k-th polynomial at x

    with function jf ( k : integer32; x : Vector ) return Vector;
     -- returns jf[k](x), the k-th row of the Jacobian matrix at x

    with function Q ( x : Vector ) return Vector;
     -- evaluates at equations which discriminate solutions

  procedure GQ_Solve_Equation_by_Equation
              ( file : in file_type; filename : in string;
                report,step : in boolean;
                ne,nv,k : in integer32; d : in Standard_Natural_Vectors.Vector;
                witset : out Array_of_Solution_Lists; planes : out VecMat );

  -- DESCRIPTION :
  --   Applies the equation-by-equation solver to the system,
  --   evaluated by f and with Jacobian matrix in jf.
  --   Solutions which satisfy any of the Q equations are discarded.

  -- REQUIRED : witset'range = planes'range = 1..max(ne,nv).

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   filename is used to create separate output files;
  --   report   if false, then path trackers will remain silent,
  --            otherwise diagnostics along paths are reported;
  --   step     flag to ask for stepping stones;
  --   ne       number of equations of the system;
  --   nv       number of variables of the system;
  --   k        index of the current equation, if k > 0, a witness
  --            set is expected for the first k equations of p;
  --   d        degrees of the polynomials in the system.

  -- ON RETURN :
  --   witset   witness sets, the i-th entry is the i-dimensional
  --            solution set of the system p;
  --   planes   planes(i) cuts out the points in witset(i).

end Equation_by_Equation_Solvers;
