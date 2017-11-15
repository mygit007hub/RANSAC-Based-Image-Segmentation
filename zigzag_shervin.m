function [ matrix_zigzag_order ] = zigzag_shervin( block_size )
 
%%%% this program generates a 2D matrix of size N*N which contains the
%%%% zigzag order of elements

N=block_size;
X=zeros(N);

for k=2:N+1
        number_previous_elements = (k-2)*(k-1)/2;
        temp = (k-2)*(k-1)/2+1;
    if ( rem(k,2) == 1 )
        for i=1:k-1
        j=k-i;
        X(i,j) = temp;
        temp=temp+1;
    end

    else
        for j=1:k-1
        i=k-j;
        X(i,j) = temp;
        temp=temp+1;
        end
    end
end

for k=N+2:2*N   %%%% 2*N+1-k elements in each diagonal
        num_dia=2*N+1-k;
        number_previous_elements = N^2 - (2*N-k+2)*(2*N-k+1)/2;
        temp = number_previous_elements+1;
    if ( rem(k,2) == 0 )
        for i=N: -1: N-num_dia+1
        j=k-i;
        X(i,j) = temp;
        temp=temp+1;
    end

    else
        for j=N: -1: N-num_dia+1
        i=k-j;
        X(i,j) = temp;
        temp=temp+1;
        end
    end
end
           
matrix_zigzag_order=X;

end

