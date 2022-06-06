%% File selection
clearvars
data = readcell('Longitudinal data_shape modelling_subject info.xlsx','Range','A2');
files = data(:,1);
N = size(files,1);
n = N/2;
for i = 1:N
    files{i} = [files{i} '.ply'];
end

%% Read files

% Stradwin settings

% Resolution: medium
% Smoothing: standard
% Strength: medium
% Prioritise: objects
% Do not triangulate ends
% Compounding type: max

tic
poly = cell(N,1);
conn = cell(N,1);
clouds = cell(N,1);
sizes = zeros(N,1);
volumes = zeros(N,1);
edgelength = 1;
parfor i = 1:N
    disp(['Reading files (',num2str(i),'/',num2str(N),')'])
    [B,A] = plyread(files{i},'tri');
    [poly{i},conn{i},~,~] = remesher(A,B,edgelength,1);
    if data{i,10} == 'R'
        % Do nothing
    else
        % Mirror about X axis
        poly{i}(:,1) = -poly{i}(:,1);
    end
    m = min(poly{i},[],1);
    M = max(poly{i},[],1);
    M = M-m;
    volumes(i) = prod(M);
    sizes(i) = length(poly{i});
    clouds{i} = pointCloud(poly{i});
end
[~,reference] = min(sizes);
disp('All files loaded');
toc

%% Mean shape construction
if ~exist('Mean_shape_pairwise_HQ.ply','file')
    % Get template
    template_cloud = pcread('Mean_shape_HQ.ply');
    template_shape = template_cloud.Location;
    load('mean_conn_HQ.mat');
    [template_shape,connectivity,~,~] = remesher(template_shape,connectivity,1.25,1);
    % Nonrigid registration #1
    reg1 = cell(N,1);
    parfor i = 1:N
        disp(['Non-rigid registration 1 (',num2str(i),'/',num2str(N),')']);
        reg1{i} = nonrigidICPv1(clouds{i}.Location,template_shape,conn{i},connectivity,10,0);
    end
    % Rigid registration
    parfor i = 1:N
        [~,reg1{i},trans] = procrustes(template_shape,reg1{i},'scaling',false,'reflection',false);
    end
    % Average shape
    tmp = cat(3,reg1{:});
    mean_shape = mean(tmp,3);
    pcwrite(pointCloud(mean_shape),'Mean_shape_pairwise.ply');
else
    mean_cloud = pcread('Mean_shape.ply');
    mean_shape = mean_cloud.Location;
    load('mean_conn.mat');
end
%% Nonrigid registration #2 (to mean shape)
reg2 = cell(N,1);
parfor i = 1:N
    disp(['Non-rigid registration 2 (',num2str(i),'/',num2str(N),')']);
    reg2{i} = nonrigidICPv1(clouds{i}.Location,mean_shape,conn{i},connectivity,10,0);
end

%% Visualization 
% % Volumetric scaling
% tmp = cell(size(reg2));
% parfor i = 1:N
%     sf(i) = (volumes(reference)/volumes(i))^(1/3);
% %     sf(i) = (lengths(reference)/lengths(i))/sqrt(3);
%     T = [sf(i) 0 0 0; ...
%         0 sf(i) 0 0; ...
%         0 0 sf(i) 0; ...
%         0 0 0 1];
%     tform = affine3d(T);
%     tmp{i} = pctransform(pointCloud(reg2{i}),tform);
%     reg2{i} = tmp{i}.Location
% end
% clear tmp
% Rigid registration
vis = cell(N,1);
parfor i = 1:N
    [~,reg3{i},trans] = procrustes(mean_shape,reg2{i},'reflection',false);
    tmp = trans.b*(clouds{i}.Location)*trans.T + repmat(trans.c(1,:),length(clouds{i}.Location),1);
    vis{i} = pointCloud(tmp);
end

% %% Figures of individual muscles
% close all
% for i = 1:n
%     figure
%     pcshow(vis{i});
%     axis equal
% end

%% Distance calculation
parfor i = 1:N % Closest vertices
    distances(:,i) = distfield(mean_shape,vis{i}.Location);
end
avg_dist = mean(distances,2);

%% vis
dist = pointCloud(mean_shape);
dist.Intensity = avg_dist;

figure('Name','Mean shape')
pcshow(dist)
axis equal

%% Grouping

groups = zeros(N,6);
cats = categorical(groups);
tmp = discretize(cell2mat(data(:,4)),2:2:10,'categorical',{'[2-4]','[4-6]','[6-8]','[8-10]'}); % Age
cats(:,1) = categorical(tmp,'Ordinal',false);
cats(:,2) = categorical(data(:,5)); % Gender
cats(:,3) = categorical(data(:,6)); % CP type
cats(:,4) = categorical(data(:,7)); % GMFCS
cats(:,5) = categorical(data(:,10)); % Leg

groups(:,1) = double(cats(:,1))-1; % Age
groups(:,2) = double(categorical(data(:,5))); % Gender
groups(:,3) = double(categorical(data(:,6))); % CP type
groups(:,4) = double(categorical(data(:,7))); % GMFCS
groups(:,5) = double(categorical(data(:,10))); % Leg
tmp = groups(:,3)-3;
groups(tmp < 0,6) = -1;
groups(:,6) = groups(:,6)+2; % CP & TD
tmp = discretize(groups(:,6),1:0.5:2,'categorical',{'CP','TD'}); % CP & TD
cats(:,6) = categorical(tmp,'Ordinal',false);
varnames = {'Age','Gender','Diagnosis','GMFCS','Leg','Subject group'};

ncat = max(groups);
ngroups = size(groups,2);

% Exclusions
groups(:,[2 3 4 5]) = [];
varnames([2 3 4 5]) = [];
ngroups = size(groups,2);

%% Average shapes

avg_shp = cell(max(ncat),size(groups,2));
count = cell(size(avg_shp));
ages = cell(size(avg_shp));
for i = 1:size(groups,2)
    for j = 1:ncat(i)
        Idx = groups(:,i) == j;
        tmp = reg3(Idx);
        count{j,i} = nnz(Idx);
        tmp = cat(3,tmp{:});
        avg_shp{j,i} = pointCloud(mean(tmp,3)); % Average shape for each group
        ages{j,i} = mean(cell2mat(data(Idx,4))); % Average age for each group
        
    end
end
% Figures
for i = 1:size(groups,2)
    catnames = string(unique(cats(:,i)));
    for j = 1:ncat(i)
        if (i == 6 && j == 2)||(i == 3 && j == 3)||(i == 4 && j == 1) %% TD does not get colour
            color = zeros(1,length(mean_shape));
        else
            color = distfield(avg_shp{3,3}.Location,avg_shp{j,i}.Location); % Colour based on distance
            avg_shp{j,i}.Intensity = color;
        end
        figure('WindowState','maximized')
        trisurf(connectivity,-avg_shp{j,i}.Location(:,1),avg_shp{j,i}.Location(:,2), ...
            -avg_shp{j,i}.Location(:,3),color,'Edgecolor','none');
        title(strcat(catnames(j),", N = ",string(count{j,i}),", Age = ",sprintf('%0.2f',ages{j,i})));
        caxis([-1.5 1.5]);
        set(gca,'XColor','none','YColor','none','ZColor','none')
        set(gca,'View',[80.1,0.19])
        grid off
        axis equal
        set(gcf,'Color',[1 1 1])
        light
        lighting phong;
        colormap('jet')
        saveas(gcf,strcat("avg shp_",catnames(j),".png"));
    end
end

%% ANOVA
etas = zeros(ngroups,length(mean_shape));
tb1 = cell(ngroups+3,7,length(mean_shape));
p = zeros(ngroups,length(mean_shape));
SS = zeros(1,ngroups);
for i = 1:length(mean_shape)
    [p(:,i),tb1(:,:,i)] = anovan(distances(i,groups(:,2)==1),groups(groups(:,2)==1,1),'display','off');
    for j = 1:ngroups
        SS(j) = tb1{1+j,2,i};
    end
    SSE = tb1{2 + ngroups,2,i};
    etas(:,i) = SS./(SS+SSE);
end
p = 1./p;

% Figures
for i = 1:ngroups
    % posterior
    color = p(i,:)';
    figure('WindowState','maximized')
    trisurf(connectivity,-mean_shape(:,1),mean_shape(:,2), ...
        -mean_shape(:,3),color,'Edgecolor','none');
    title(strcat(varnames{i},", posterior"));
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light
    lighting phong;
    colormap('jet')
%     caxis([0.9 1])
    set(gca,'ColorScale','log')
    saveas(gcf,strcat(varnames(i),"_p_posterior.png"));
    
    % anterior
    figure('WindowState','maximized')
    trisurf(connectivity,mean_shape(:,1),-mean_shape(:,2), ...
        -mean_shape(:,3),color,'Edgecolor','none');
    title(strcat(varnames{i},", anterior"));
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    set(gca,'ColorScale','log')
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light('Position',[1 0 1])
    lighting phong;
    colormap('jet')
%     caxis([0.9 1])
    saveas(gcf,strcat(varnames(i),"_p_anterior.png"));
    
    % Etas
    % posterior
    color = etas(i,:)';
    figure('WindowState','maximized')
    trisurf(connectivity,-mean_shape(:,1),mean_shape(:,2), ...
        -mean_shape(:,3),color,'Edgecolor','none');
    title(strcat(varnames{i},", posterior"));
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light
    lighting phong;
    colormap('jet')
    caxis([0 0.25])
    saveas(gcf,strcat(varnames(i),"_eta_posterior.png"));
    
    % anterior
    figure('WindowState','maximized')
    trisurf(connectivity,mean_shape(:,1),-mean_shape(:,2), ...
        -mean_shape(:,3),color,'Edgecolor','none');
    title(strcat(varnames{i},", anterior"));
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light('Position',[1 0 1])
    lighting phong;
    colormap('jet')
    caxis([0 0.25])
    saveas(gcf,strcat(varnames(i),"_eta_anterior.png"));
end


%% SSM
% Training data
Xdata = reg3{1}(:,1);
Ydata = reg3{1}(:,2);
Zdata = reg3{1}(:,3);
for i = 2:N
    Xdata = [Xdata,reg3{i}(:,1)];
    Ydata = [Ydata,reg3{i}(:,2)];
    Zdata = [Zdata,reg3{i}(:,3)];
end

[ssmV,Eval,Evec,MEAN,PCcum,Modes] = SSMbuilder(Xdata,Ydata,Zdata);

%% Figures
nmax = 20;
SD = 3;
for i = 1:nmax
    figure
    trisurf(connectivity,-mean_shape(:,1),mean_shape(:,2), ...
        -mean_shape(:,3),'Edgecolor','none'); % Mean shape
    hold;
    colormap bone
    title(strcat('PC ',string(i)));
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light('Position',[1 0 1])
    lighting phong;
    
    offset = 50;
    shape = MEAN+SD*ssmV(:,i);
    shape = reshape(shape,[],3);
    trisurf(connectivity,-shape(:,1),shape(:,2)+offset, ...
        -shape(:,3),'Edgecolor','none'); % 2 SD
    shape = MEAN-SD*ssmV(:,i);
    shape = reshape(shape,[],3);
    trisurf(connectivity,-shape(:,1),shape(:,2)-offset, ...
        -shape(:,3),'Edgecolor','none'); % -2 SD
%     saveas(gcf,strcat('PC ',string(i),'.png'));
end

figure
plot((Eval(1:nmax)./sum(Eval))*100);
set(gcf,'Color',[1 1 1])
title('Variation explained by PCs');
ylabel('% of total');
xlabel('PC number');
saveas(gcf,'PCcum.png');

figure
plot(PCcum(1:50)*100);
set(gcf,'Color',[1 1 1])
title('Cumulative variation explained by PCs');
ylabel('% of total');
xlabel('PC number');

%% SSM ANOVA
p = zeros(ngroups,nmax);
for i = 1:nmax
    p(:,i) = anovan(Modes(:,i),groups,'display','off');
end
p = 1-p;

%%
figure
hold on;
bar(p')
plot(linspace(0,20.5,nmax),ones(1,nmax)*0.95,'r-');
set(gcf,'Color',[1 1 1])
legend(varnames);
title('1-p for each PC by category');
xlim([0 20.5]);
ylabel('1-p');
xlabel('PC number');