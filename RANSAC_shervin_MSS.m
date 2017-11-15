function [ poly_weight, EIR, avg_err, foreground_map, background_rec, no_inner_loop] = RANSAC_shervin_MSS( input_block, polynomial_basis, inlier_ratio, maximum_distortion, maximum_iteration )

%%%%% this is the third version of RANSAC algorithm which is modified on
%%%%% july 9, 2014

X=input_block;
[s1,s2]=size(X);
P=polynomial_basis;
ro=inlier_ratio;
nop=size(P,2);
[s3,s4]=size(P);
alpha=zeros(s4,1);
X_1d = reshape(X,s1*s2,1);
no_iteration=0;
Err=zeros(s1*s2,1);
% kmeans_ini_seg_fore=zeros(size(X)); %%% initial segmentation by k-means
estimated_inlier_ratio=0;
temp=0;
r_max_con=[];
m_max_con=0;
used_iteratoion=0;

while ( (no_iteration < maximum_iteration) && (estimated_inlier_ratio < ro)  )
    
       used_iteratoion=used_iteratoion+1;
       ind=randperm(s1*s2,nop);
       A=P(ind,:);
       while ( abs(det(A)) < 10^(-60) )
       ind=randperm(s1*s2,nop);
       A=P(ind,:);
       end
       b=X_1d(ind,:);
       alpha=A\b;
       Err=P*alpha - X_1d;
       Err=abs(Err);
       [r,c,v] = find(Err < maximum_distortion);
       estimated_inlier_ratio= length(r)/(s1*s2);    %%% here we get initial background pixels using minimum sample set 
       if ( length(r) > m_max_con)
           m_max_con= length(r);
           r_max_con= r;
       end
       no_iteration=no_iteration+1;
       avg_err = sum(Err(:))/(s1*s2);
       temp=temp+1;
end

%%%%% extracting initial inliers
r_previous = 0;
r_after = r_max_con;
while ( isequal(r_previous, r_after)~=1)
X_f_inl= X_1d(r_after,:); %%% extracting first inliers from data
P_f_inl= P(r_after,:);    %%% correspong polynomial points
alpha_sr= (P_f_inl'*P_f_inl)\(P_f_inl'*X_f_inl);
r_previous= r_after;
Er_final= P*alpha_sr - X_1d;
Er_final= abs(Er_final);
[ r_after, c_after, v_after] = find(Er_final < maximum_distortion);
estimated_inlier_ratio_f= length(r_after)/(s1*s2);
avg_err_f = sum(Er_final(:))/(s1*s2);
end


%%%%% extracting inliers
background_index= r_after;
poly_weight= alpha_sr;
foreground_map=zeros(s1,s2);
foreground_map_1d=ones(s1*s2,1);
foreground_map_1d(r_after)=0;
foreground_map = reshape(foreground_map_1d, s1, s2);
background_1D= P*alpha_sr;
background_rec = reshape(background_1D,s1,s2); 
EIR=estimated_inlier_ratio_f;
no_inner_loop= temp;

end

