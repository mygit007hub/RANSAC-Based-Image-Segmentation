function B = unif_quan( inp , q )

A=inp;
[s1 s2]=size(A);
B=zeros(size(A));
for i=1:s1
    for j=1:s2
        B(i,j)= floor( (A(i,j) + q/2 )/q) * q ;
    end
end


end

