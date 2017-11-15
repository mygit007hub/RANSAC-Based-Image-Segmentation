%%%%%%%%%%%%%%%%% This m-file implements a modified version RANSAC based image foreground segmentation
%%%%%%%%%%%%%%%%% This work has been done on July 31, 2014
clc
close all
clear all



[AA(:,:,1), AA(:,:,2), AA(:,:,3)] = yuv_read('MissionControlClip3_380.yuv', 8, 1920, 1080, 0, 'yv24'); %%%11th watch image in paper with no quad-tree decomposition
 AA= AA(3*64:5*64-1, 27*64:29*64-1,:);  %%%%  N=64, M=12, th= 27

% [AA(:,:,1), AA(:,:,2), AA(:,:,3)] = yuv_read('SlideShow_1280x720_110.yuv', 8, 1280, 720, 0, 'yv24'); %%% 12 linear change, range dis, no quad
% AA= AA(7*64:10*64-1, 10*64:13*64-1,:); %%%% N=64, M=12, th= 24



input_image= double(AA(:,:,1));
U_com= double(AA(:,:,2));
V_com= double(AA(:,:,3));
[sh sv] = size(input_image);
N= 64;
Background_filled= zeros(N,N,3);
Foreground_filled= zeros(N,N,3);
Original_image =   zeros(N,N,3);
Original_image = double(AA);



zig_zag= zigzag_shervin( N );
DCT_1d = zeros(N*N);  %%% each dct basis in one column
for k=0:N*N-1
    q= floor(k/N);
    p= k - q*N;
    DCT_Bases(:,:,k+1)= shervin_2dctillus( p, q, N );
    DCT_1d(:,zig_zag(p+1,q+1))= reshape(DCT_Bases(:,:,k+1),N*N,1);
    ind_check(k+1)=zig_zag(p+1,q+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANSAC using first K DCT basis
type_block=zeros(size(input_image)); %%% image which shows each blocks is processed by which algorithm, variance, least sqaure or ransac
M= 12; %%%% number of used dct basis
P= DCT_1d(:,1:M);
Foreground_map=zeros(size(input_image));
inlier_ratio= 0.9;  %%%% adaptive inlier ratio initialization
max_std_block= 255/2;
number_point_step=size(P,2);
max_dis= 27;  %%%% maximum_distortion
max_std_background= 5;
max_iter= 200;
estimated_inlier_ratio = zeros(sh/N , sv/N);
no_inner_loop = zeros(sh/N , sv/N);
avg_err = zeros(sh/N , sv/N);
temppp=0;
maximum_ls_err=5;
small_inlier = zeros(size(input_image));
aa_chose = zeros(sh/N , sv/N);
salam_check = zeros(sh/N , sv/N);
ls_min_corr_ratio= 0.99; %%% minimum number of pixels which should be predicted correctly with least square

T1 = 20; %%%% threshold for the number of colors in a block
T2 = 3;  %%%% threshold for pixel intensity standard variation

cur_for_ind= zeros(N,N);
cur_back_ind= zeros(N,N);

tic;

const_back=0;
least_squa=0;
text_back= 0;
RANSAC= 0;
check_ifsame= zeros(sh/N, sv/N);

for i = 1:sh/N
    for j=1:sv/N
        
        X= input_image((i-1)*N+1:i*N,(j-1)*N+1:j*N);
        
        [ foreground_map, recons_back, estimated_inr, r_index ] = leastsquare_linearchange_2DCT( X , P , maximum_ls_err );

        if       (   std(X(:)) <= max_std_background  )  %%%% smooth background category, first category
            Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 0*Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N);
            type_block((i-1)*N+1:i*N,(j-1)*N+1:j*N)=0;  %%%% black for smooth background 
            estimated_inlier_ratio(i,j)=1;
            const_back=const_back+1;
            x1d_1=reshape(X,N*N,1);
            poly_weight_1= (P'*P)\(P'*x1d_1);
            back_1= P*poly_weight_1;
            Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)= reshape(back_1,N,N);
            
        elseif    ( ( std(X(:)) ~= 0 )  &&  ( estimated_inr >= ls_min_corr_ratio) )        %%%% least square fitting,second category
            Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)= zeros(N,N);
            Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)= recons_back;
            type_block((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 60;  %%%% dark gray for lleast-sqaure
            least_squa=least_squa+1;
            estimated_inlier_ratio(i,j)=estimated_inr;
            
        elseif    ( ( length(unique(X)) < T1 ) || ( std(X(:)) == 0 ) )                    %%%% third category, text/graphics
            quantization_step=1;
            [ Background_map((i-1)*N+1:i*N,(j-1)*N+1:j*N), true_map((i-1)*N+1:i*N,(j-1)*N+1:j*N), num_col(i,j) ] = quan_in2level( X, quantization_step, T1 );
            aa=Background_map((i-1)*N+1:i*N,(j-1)*N+1:j*N);
            check_ratio= sum(sum(aa))/(N*N);
            if ( check_ratio>=0.4)
            Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 1-Background_map((i-1)*N+1:i*N,(j-1)*N+1:j*N);
            estimated_inlier_ratio(i,j)= sum(sum(Background_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)))/(N*N);
            Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)= Background_map((i-1)*N+1:i*N,(j-1)*N+1:j*N).*X;
            type_block((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 120;  %%%% bright gray for least-sqaure
            text_back=text_back+1;
            aa_3= reshape( aa, N*N, 1);
            [r_3, c_3, v_3]= find(aa_3 == 1);
            x1d_3 = reshape(X, N*N,1);
            x1d_3_r = x1d_3(r_3,:);
            P_3= P(r_3,:);
            alpha_3 = (P_3'*P_3)\(P_3'*x1d_3_r);
            back_1d_3 = P*alpha_3;
            Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)= reshape(back_1d_3,N,N);
            else
            [ poly_weight, estimated_inlier_ratio(i,j), avg_err(i,j), foreground_map, background_rec, no_inner_loop(i,j)] = RANSAC_shervin_MSS( X, P, inlier_ratio, max_dis, max_iter);
            type_block((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 255;
            Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)=foreground_map;
            Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)=background_rec; 
            end
        elseif( ( std(X(:)) ~= 0 )  &&  ( estimated_inr >= 0.8) )   %%%% mode 4, iterative least square for those that have more than 0.8
            r_previous = 0;
            r_after = r_index;
            X_1d= reshape(X,N*N,1);
            while ( isequal(r_previous, r_after)~=1)
            X_f_inl= X_1d(r_after,:); %%% extracting first inliers from data
            P_f_inl= P(r_after,:);    %%% corresponding polynomial points
            alpha_sr= (P_f_inl'*P_f_inl)\(P_f_inl'*X_f_inl);
            r_previous= r_after;
            Er_final= P*alpha_sr - X_1d;
            Er_final= abs(Er_final);
            [ r_after, c_after, v_after] = find(Er_final < maximum_ls_err);
            estimated_inlier_ratio_f= length(r_after)/(N*N);
            avg_err_f = sum(Er_final(:))/(N*N);
            end
            foreground_map_ils= ones(N*N,1);  %%%ils: iterative least square
            foreground_map_ils(r_after)=0;
            Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)= reshape(foreground_map_ils, N, N);
            back_rec_ils= P*alpha_sr;
            Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)= reshape(back_rec_ils, N, N);
            type_block((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 180;
            least_squa=least_squa+1;
            estimated_inlier_ratio(i,j)= estimated_inlier_ratio_f;
            
        else
        [ poly_weight, estimated_inlier_ratio(i,j), avg_err(i,j), foreground_map, background_rec, no_inner_loop(i,j)] = RANSAC_shervin_MSS( X, P, inlier_ratio, max_dis, max_iter);
        type_block((i-1)*N+1:i*N,(j-1)*N+1:j*N)= 255;
        RANSAC=RANSAC+1;
        Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N)=foreground_map;
        Background_rec((i-1)*N+1:i*N,(j-1)*N+1:j*N)=background_rec;
        end           

    cur_for_ind((i-1)*N+1:i*N,(j-1)*N+1:j*N)= Foreground_map((i-1)*N+1:i*N,(j-1)*N+1:j*N);

    temppp= temppp+1
    end
end

total_time = toc;

SE=[0 1 0; 1 1 1; 0 1 0];
Foreground_map1 = imclose(Foreground_map, SE);
Foreground_map3 = isolated_remove( Foreground_map1, 4 );


cur_for_ind= Foreground_map3;
Background_map= 1- Foreground_map3;
Background_rec(Background_rec<0)=0;
Background_rec(Background_rec>255)=255;
cur_back_ind= 1 - cur_for_ind;


imshow(cur_for_ind)
title('current fore ind');







