function dct2_func = dct_2d( input_image )

A= input_image;
[s1 s2]=size(A);
B=zeros(size(A));
A_mid=zeros(size(A));
dct2_func=zeros(size(A));

for i=1:s1
    A_mid(i,:)=dct_1d(A(i,:));
end
for j=1:s2
    B(:,j)=dct_1d(A(:,j));
end

dct2_func=B;

end