function [v] = gs(A)
[Arow,Acolumn] = size(A); 
v(:,1) = A(:,1)/norm(A(:,1));
for k = 2:Acolumn
    for i = k:Acolumn
        A(:,i) = A(:,i)-A(:,i)'*v(:,k-1)*v(:,k-1);
    end
    v(:,k) = A(:,k)/norm(A(:,k));
end