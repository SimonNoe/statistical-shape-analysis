%% Init
% Training data
load('Pairwise HQ.mat')
N = length(reg3);

% Test data
load('testdata.mat');
n = length(testset);

%% SSM
shape = cell(N,n,10);
X = cell(N,n,10);
RMSE = zeros(N,n,10);
vols = zeros(N,n,10);
for k = 1:10
    % Arrange data
    Idx = randperm(N);
    Xdata = reg3{Idx(1)}(:,1);
    Ydata = reg3{Idx(1)}(:,2);
    Zdata = reg3{Idx(1)}(:,3);
    for i = 2:N
        Xdata = [Xdata,reg3{Idx(i)}(:,1)];
        Ydata = [Ydata,reg3{Idx(i)}(:,2)];
        Zdata = [Zdata,reg3{Idx(i)}(:,3)];
    end
    
    options = optimset('Display','off');
    parfor i = 2:N
        disp(['Creating model (',num2str(i),'/',num2str(N),')'])

        [ssmV,~,~,MEAN,~,~] = SSMbuilder(Xdata(:,1:i),Ydata(:,1:i),Zdata(:,1:i));

        % SSM validation
        for j = 1:n
            % Shape vector
            [~,tmp] = procrustes(reshape(MEAN,[],3),testset{j},'reflection',false);
            X{i,j,k} = mldivide(ssmV,reshape(tmp,[],1)-MEAN);

            % Volumes
            tmp = MEAN+(ssmV*X{i,j,k});
            shape{i,j,k} = reshape(tmp,[],3);
            [RMSE(i,j,k),tmp] = procrustes(shape{i,j,k},testset{j},'reflection',false);
            shp1 = alphaShape(shape{i,j,k});
            shp2 = alphaShape(tmp);
            vols(i,j,k) = ((volume(shp2) - volume(shp1))/volume(shp2))*100;
        end
    end
end

%% Figure
avg_RMSE = mean(mean(RMSE,3),2);
avg_volume = mean(mean(vols,3),2);

figure
yyaxis left
plot(avg_RMSE);
set(gcf,'Color',[1 1 1])
ylabel('RMSE (mm)');
xlabel('Number of muscles in model');
xlim([2 N]);

yyaxis right
plot(avg_volume);
ylabel('Volume (%)');
legend('Average RMSE','Average volume difference');

% %%
% figure
% [~,tmp] = procrustes(shape{N,20},testset{20},'reflection',false);
% pcshowpair(pointCloud(tmp),pointCloud(shape{N,20}));