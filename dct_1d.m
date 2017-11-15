function dct_func = dct_1d( input_vector )

A= input_vector;
N=length(A);

B=zeros(size(A));
W=zeros(size(A));
W(1)=sqrt(1/N);
W(2:N)=sqrt(2/N);

for i=1:N
    sum=0;
    for j=1:N
        sum=sum+A(j)*cos(  pi*(2*j-1)*(i-1)/(2*N)  );
    end
    B(i)=W(i)*sum;
end

dct_func=B;

end
