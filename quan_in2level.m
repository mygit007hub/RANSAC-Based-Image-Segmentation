function [ background_map, true_map, num_col ] = quan_in2level( input_block, quantization_step, max_num_col )

%%%%% this function tries to do quantization and find out how many dominant
%%%%% colors is available in a block and then find the background and foreground
%%%% true_map denotes a map which is 1 if we have up to 4 different colors,  else zero

X= input_block;
q= quantization_step;
[s1,s2]= size(X);
true_map= zeros(size(X));
background_map= zeros(size(X));

B= floor( (X + q/2 )/q) * q ;

A=unique(B);
num_col= length(A(:));
m=zeros( 1, length(A(:)) );

if (  (1<= length(A(:)) ) &&  ( length(A(:))<= max_num_col) )
    for i=1:length(A(:))
    m(i) = sum( B(:)== A(i) );
    end  
    [c, ind]= max(m);
    background_map( B==A(ind) )=1;
    true_map= 1+ true_map;
end


end

