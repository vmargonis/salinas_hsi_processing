clear
format compact
close all

load Salinas_Data

[p,n,l]=size(Salinas_Image); % Size of the Salinas cube

% Making a two dimensional array whose rows correspond to the pixels and
% the columns to the bands, containing only the pixels with nonzero label.
X_total=reshape(Salinas_Image, p*n,l);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);   %This contains 1 in the positions corresponding to 
                   %pixels with known class label
X=X_total(existed_L,:);
[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)

%labels >0
L0 = L(L>0);

%--------------------------------------------------------------------------
% PCA
n_components = 7;
[~,~,~,Y,~] = pca_fun(X', n_components);

% uniform scaling on Y:
a = mean(mean(abs(Y)));
Y = Y / a;

% FUZZY C-MEANS
m=8;

% min-max init
[~, theta_init] = most_dist_repre(Y, m);
theta_init = theta_init';

% parameters and execution
q=4;
max_iter = 1000;
tic;
[theta, U] = fuzzy(Y, m, theta_init, q, max_iter);
run_time=toc;
[~, bel] = max(U, [], 2);
%--------------------------------------------------------------------------

% cluster and total accuracy
[conf_mat, acc, bel_new] = accuracy(L0, bel, m);
cluster_acc = (diag(conf_mat) ./ sum(conf_mat, 2))';
conf_mat, cluster_acc, acc

% some custom colors...
c = struct('rr', [0.9, 0.19, 0.19], ...  
    'bb', [0.29, 0.54, 0.74], ... 
    'um', [0.08, 0.12, 0.41], ... 
    'br', [0.65, 0.57, 0.34], ... 
    'gl', [0.83, 0.7, 0.78] ); 

% plot in the first two principal components
theta=theta';

figure(2), hold on
figure(2), scatter(Y(1,bel_new==1), Y(2,bel_new==1), 15, 'b')
figure(2), scatter(Y(1,bel_new==2), Y(2,bel_new==2), 15, c.bb)
figure(2), scatter(Y(1,bel_new==3), Y(2,bel_new==3), 15, c.rr)
figure(2), scatter(Y(1,bel_new==4), Y(2,bel_new==4), 15, c.um)
figure(2), scatter(Y(1,bel_new==5), Y(2,bel_new==5), 15, c.br)
figure(2), scatter(Y(1,bel_new==6), Y(2,bel_new==6), 15, c.gl)
figure(2), scatter(Y(1,bel_new==7), Y(2,bel_new==7), 15, 'm')
figure(2), scatter(Y(1,bel_new==8), Y(2,bel_new==8), 15, 'g')
figure(2), scatter(theta(1,:),theta(2,:), 25 ,'MarkerEdgeColor','k',...
              'MarkerFaceColor','r')
figure(2), axis equal
figure(2), hold off

% plot cube
cl_label = bel_new';
cl_label_tot=zeros(p*n,1);
cl_label_tot(existed_L)=cl_label;
im_cl_label=reshape(cl_label_tot,p,n);
figure(3), imagesc(im_cl_label)