clear
format compact
close all

load Salinas_Data

[p,n,l]=size(Salinas_Image); % Size of the Salinas cube

%Depicting the bands of the Salinas cube
for i=1:l
    figure(1); imagesc(Salinas_Image(:,:,i))
    pause(0.1)
end

% Making a two dimensional array whose rows correspond to the pixels and
% the columns to the bands, containing only the pixels with nonzero label.
X_total=reshape(Salinas_Image, p*n,l);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);   %This contains 1 in the positions corresponding to 
                   % pixels with known class label
X=X_total(existed_L,:);
[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)

X_all=reshape(Salinas_Image,p*n,l)';

%Choose the data points that have available ground truth
Q_ground=find(L>0);
X1=X_all(:,Q_ground);
L=L(Q_ground);

%Perform PCA 
[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(X1,3);

% Creating the first three principal components
Y=Y-min(min(Y)); % Making Y>=0
PCmat=zeros(3,p*n);
PCmat(1,Q_ground)=Y(1,:);
PCmat(2,Q_ground)=Y(2,:);
PCmat(3,Q_ground)=Y(3,:);
PCmat=PCmat';
PCmat=reshape(PCmat,p,n,3);
for i=1:3
    PCmat(:,:,i) = (PCmat(:,:,i)-(min(min(PCmat(:,:,i))) )) / ...
    (max(max(PCmat(:,:,i)) - min(min(PCmat(:,:,i)))));
end
% Depicting the first PC's of the salinas image cube
for i=1:3
    figure(i+11), imagesc(PCmat(:,:,i)); colorbar; axis off; axis image
    pause(1)
end

%labels >0
L0 = L(L>0);

% PCA
n_components = 7;
[~,~,E,Y,~] = pca_fun(X', n_components);

% uniform scaling on Y:
a = mean(mean(abs(Y)));
Y = Y / a;

m=8;
% min-max init
[~, theta_init] = most_dist_repre(Y, m);

% some custom colors...
c = struct('rr', [0.9, 0.19, 0.19], ...  
    'bb', [0.29, 0.54, 0.74], ... 
    'um', [0.08, 0.12, 0.41], ... 
    'br', [0.65, 0.57, 0.34], ... 
    'gl', [0.83, 0.7, 0.78] ); 

%PLOT DATA SET IN 2D and CUBE

Y_1 = Y(:, (L0 == 1));
Y_2 = Y(:, (L0 == 2));
Y_3 = Y(:, (L0 == 3));
Y_4 = Y(:, (L0 == 4));
Y_5 = Y(:, (L0 == 5));
Y_6 = Y(:, (L0 == 6));
Y_7 = Y(:, (L0 == 7));
Y_8 = Y(:, (L0 == 8));

figure(1), hold on
figure(1), scatter(Y_1(1,:), Y_1(2,:), 15, 'b', 'filled')
figure(1), scatter(Y_2(1,:), Y_2(2,:), 15, c.bb, 'filled')
figure(1), scatter(Y_3(1,:), Y_3(2,:), 15, c.rr, 'filled')
figure(1), scatter(Y_4(1,:), Y_4(2,:), 15, c.um, 'filled')
figure(1), scatter(Y_5(1,:), Y_5(2,:), 15, c.br, 'filled')
figure(1), scatter(Y_6(1,:), Y_6(2,:), 15, c.gl, 'filled')
figure(1), scatter(Y_7(1,:), Y_7(2,:), 15, 'm', 'filled')
figure(1), scatter(Y_8(1,:), Y_8(2,:), 15, 'g', 'filled')
figure(1), axis equal
figure(1), hold off

cl_label = L0';
cl_label_tot=zeros(p*n,1);
cl_label_tot(existed_L)=cl_label;
im_cl_label=reshape(cl_label_tot,p,n);
figure(2), imagesc(im_cl_label)

% PLOT REPRESENTATIVES

figure(3), hold on
figure(3), scatter(Y(1,:), Y(2,:), 15, c.bb)
figure(3), scatter(theta_init(1,:),theta_init(2,:), 25 , ...
                'MarkerEdgeColor','k','MarkerFaceColor','r')
figure(3), axis equal
figure(3), hold off

