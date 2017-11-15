
function [ foreground_map, recons_back, estimated_inlier_ratio, r_index ] = leastsquare_linearchange_2DCT( input_block , polynomial_basis , maximum_distortion )

%%%% input_block is a 2D block
%%%% polynomial_basis contains a matrix which its columns correspond to different basis

X= input_block;
[s1,s2]=size(X);
% Y= dct_2d(X);
Y= dct2(X);
Z= zeros(s1,s2);
Z(1:4,1)= Y(1:4,1);
Z(1:3,2)= Y(1:3,2);
Z(1:2,3)= Y(1:2,3);
Z(1,4)= Y(1,4);

X_rec= idct2(Z);
Err= X_rec - X;
Err_1d= reshape(Err, s1*s2,1);
Err_1d= abs(Err_1d);
[r,c,v] = find(Err_1d < maximum_distortion);
estimated_inlier_ratio= length(r(:))/(s1*s2);    %%% here we get initial background pixels using minimum sample set 
r_index=r;

recons_back = X_rec;
foreground_map_1d= ones(s1*s2,1);
foreground_map_1d(r)=0;
foreground_map = reshape(foreground_map_1d, s1, s2);


end

