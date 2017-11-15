function [ isolated_remove ] = isolated_remove( input_image, neighborhood_kind )

X=input_image;
[s1,s2]=size(X);

Y= zeros(s1+2,s2+2); %%%% adding zeros columns and rows in the first and last one
Y(2:s1+1,2:s2+1)=X;
Z= Y; %%% after removing isolated points in 4-neighborhood

if (neighborhood_kind==4)
    se= [ 0 1 0; 1 1 1; 0 1 0]; %%% structural element for morphological operations
else
    se= ones(3,3);
end
neighbors=zeros(3,3);

for i=2:s1+1
    for j=2:s2+1
        neighbors= se.*Y(i-1:i+1,j-1:j+1);
        if ( nnz(neighbors)==1 )
            Z(i,j)=0;
        end
    end
end

isolated_remove= zeros(s1,s2);
isolated_remove= Z(2:s1+1,2:s2+1);

end