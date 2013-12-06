with text_io;                             use text_io;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Sample_Point_Lists;                  use Sample_Point_Lists;

package Sample_Point_Lists_io is

-- DESCRIPTION :
--   This package provides input/output routines for lists of samples.
--   The input format is the same as the output format.

  procedure get ( samples,samples_last : in out Standard_Sample_List );
  procedure get ( len,n,k : in natural32;
                  samples,samples_last : in out Standard_Sample_List );
  procedure get ( file : in file_type;
                  samples,samples_last : in out Standard_Sample_List );
  procedure get ( file : in file_type; len,n,k : in natural32;
                  samples,samples_last : in out Standard_Sample_List );
  procedure get ( samples,samples_last : in out Multprec_Sample_List );
  procedure get ( len,n,k : in natural32;
                  samples,samples_last : in out Multprec_Sample_List );
  procedure get ( file : in file_type;
                  samples,samples_last : in out Multprec_Sample_List );
  procedure get ( file : in file_type; len,n,k : in natural32;
                  samples,samples_last : in out Multprec_Sample_List );

  -- DESCRIPTION :
  --   If the parameters len, n, and k are not provided on input, then
  --   reads the length of the list, the dimension of every solution,
  --   and the number of hyperplane sections for every sample first
  --   before appending the samples to the list.

  procedure put ( samples : in Standard_Sample_List );
  procedure put ( len,n,k : in natural32; samples : in Standard_Sample_List );
  procedure put ( file : in file_type; samples : in Standard_Sample_List );
  procedure put ( file : in file_type;
                  len,n,k : in natural32; samples : in Standard_Sample_List );
  procedure put ( samples : in Multprec_Sample_List );
  procedure put ( len,n,k : in natural32; samples : in Multprec_Sample_List );
  procedure put ( file : in file_type; samples : in Multprec_Sample_List );
  procedure put ( file : in file_type;
                  len,n,k : in natural32; samples : in Multprec_Sample_List );

  -- DESCRIPTION :
  --   The parameters len,n,k are written first on file or standard output
  --   before the list of samples.

  procedure Write_Summaries
                ( file : in file_type; samples : in Standard_Sample_List );

  -- DESCRIPTION :
  --   Writes the summary (err,rco,res) for every solution on file.

end Sample_Point_Lists_io;
