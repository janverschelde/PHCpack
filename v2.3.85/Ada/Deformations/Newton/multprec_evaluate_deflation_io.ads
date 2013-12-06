with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Multprec_Evaluate_Deflation;       use Multprec_Evaluate_Deflation;

package Multprec_Evaluate_Deflation_io is

-- DESCRIPTION :
--   This package provides output facilities for the structures to
--   evaluate deflation matrices efficiently.

  procedure Write_Derivative_Operator
               ( d : in Standard_Natural_Vectors.Vector );
  procedure Write_Derivative_Operator
               ( file : in file_type;
                 d : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the derivative operator, d(i) means the
  --   i-th variable has to be derived d(i) times.

  procedure Write_Derivative_Operator
               ( d : in Standard_Natural_Vectors.Vector;
                 k : in natural32 );
  procedure Write_Derivative_Operator
               ( file : in file_type;
                 d : in Standard_Natural_Vectors.Vector;
                 k : in natural32 );

  -- DESCRIPTION :
  --   Writes the application of the derivative operator
  --   in v and d to the deflation matrix A(k).

  procedure Write_Derivative_Operator
               ( d : in Standard_Natural_Vectors.Vector;
                 k,l : in natural32 );
  procedure Write_Derivative_Operator
               ( file : in file_type;
                 d : in Standard_Natural_Vectors.Vector;
                 k,l : in natural32 );

  -- DESCRIPTION :
  --   Writes the derivative operator d applied to A(k),
  --   with the indentation of l double spaces.

  procedure Write_Spaces ( l : in natural32 );
  procedure Write_Spaces ( file : in file_type; l : in natural32 );

  -- DESCRIPTION :
  --   Writes l double spaces.

  procedure Write_Zero ( l : in natural32 );
  procedure Write_Zero ( file : in file_type; l : in natural32 );

  -- DESCRIPTION :
  --   Writes zero after l double spaces.

  procedure Write ( evt : in Eval_Tree );
  procedure Write ( file : in file_type; evt : in Eval_Tree );
  procedure Write ( evt : in Eval_Tree;
                    nv,nq,R1 : in Standard_Natural_Vectors.Vector );
  procedure Write ( file : in file_type; evt : in Eval_Tree;
                    nv,nq,R1 : in Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the evaluation tree to standard output.
  --   With the (nv,nq,R1), the dimensions of the matrices are counted
  --   and also written to output.

end Multprec_Evaluate_Deflation_io;
