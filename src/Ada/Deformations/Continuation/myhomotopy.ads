with Generic_Inline_Homotopy;
with eval_homcyc7;

package MyHomotopy is
  new Generic_Inline_Homotopy(eval_homcyc7.Homotopy_Constants,
                              eval_homcyc7.Eval_Homotopy,
                              eval_homcyc7.Diff_Homotopy,
                              eval_homcyc7.Diff_Homotopy);
