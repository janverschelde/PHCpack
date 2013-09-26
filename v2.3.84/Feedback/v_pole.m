function [diff_a] = v_pole(given, compute, dim)
% Compute the difference between given and the computed eigenvalues
for j=1:dim
   for k=1:dim
       diff=abs(given(j)-compute(k));
       if(diff<1e-5)
            diff_a(k)=diff/abs(given(j));
       end
    end
end

    