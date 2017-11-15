function basis_pq = shervin_2dctillus( p, q, N ) %%%dct bases corresponds to D(p,q)

A=zeros(N,N);

for m=0:N-1
    for n=0:N-1
        A(m+1,n+1)= cos(pi*(2*m+1)*p/(2*N)) * cos(pi*(2*n+1)*q/(2*N));
    end
end       

if (p==0)
    alpha_p = sqrt(1/N);
else
    alpha_p = sqrt(2/N);
end

if (q==0)
    alpha_q = sqrt(1/N);
else
    alpha_q = sqrt(2/N);
end

basis_pq = alpha_q*alpha_p*A;

end