%% File selection
clearvars
close all
data = readcell('subject info.xlsx','Range','A2');
files = data(:,1);
n = size(files,1);
for i = 1:length(files)
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

poly = cell(n,1);
conn = cell(n,1);
clouds = cell(n,1);
sizes = zeros(n,1);
volumes = zeros(n,1);
edgelength = 1;
parfor i = 1:n
    disp(['Reading files (',num2str(i),'/',num2str(n),')'])
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

%% Distortion
tmp = cell(n,1);
norm = logical([1 1 0 0 0 0 1 1 0 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 1]);
sigma = 10;
strength = 100;
for i = 1:n % Apply distortions
    [~,Idx] = min(clouds{i}.Location(:,3),[],1); % Point selection
    point = clouds{i}.Location(Idx,:);
    tmp{i} = pointCloud(Distortion(clouds{i}.Location,conn{i},point,0,norm(i),sigma,strength));
end
clouds = cat(1,clouds,tmp);
data(:,16) = {0}; % Distortion variable
data = cat(1,data,data);
data(n+1:end,16) = {1};
conn = cat(1,conn,conn);
volumes = cat(1,volumes,volumes);
clear tmp
n = n*2;

%% Figures of individual muscles
% close all
% for i = 1:n
%     figure
%     trisurf(conn{i},clouds{i}.Location(:,1),clouds{i}.Location(:,2), ...
%         clouds{i}.Location(:,3),'Edgecolor','none');
%     set(gca,'XColor','none','YColor','none','ZColor','none')
%     set(gca,'View',[80.1,0.19])
%     grid off
%     axis equal
%     set(gcf,'Color',[1 1 1])
%     light
%     lighting phong;
% end

%%
color = zeros(size(clouds{reference}.Location(:,1)));
figure('WindowState','maximized')
trisurf(conn{reference},-clouds{reference}.Location(:,1),clouds{reference}.Location(:,2), ...
    -clouds{reference}.Location(:,3),color,'Edgecolor','none');
colormap bone
set(gca,'XColor','none','YColor','none','ZColor','none')
set(gca,'View',[80.1,0.19])
grid off
axis equal
set(gcf,'Color',[1 1 1])
light
lighting phong;
saveas(gcf,strcat(string(strength),'_',string(sigma),'_1.png'));

figure('WindowState','maximized')
trisurf(conn{reference+n/2},-clouds{reference+n/2}.Location(:,1),clouds{reference+n/2}.Location(:,2), ...
    -clouds{reference+n/2}.Location(:,3),color,'Edgecolor','none');
colormap bone
set(gca,'XColor','none','YColor','none','ZColor','none')
set(gca,'View',[80.1,0.19])
grid off
axis equal
set(gcf,'Color',[1 1 1])
light
lighting phong;
saveas(gcf,strcat(string(strength),'_',string(sigma),'_2.png'));

%% Mean shape construction
% Nonrigid registration #1
reg1 = cell(n,1);
parfor i = 1:n
    disp(['Non-rigid registration 1 (',num2str(i),'/',num2str(n),')']);
    if i == reference
        reg1{i} = clouds{i}.Location;
    else
        reg1{i} = nonrigidICPv1(clouds{i}.Location,clouds{reference}.Location,conn{i},conn{reference},10,0);
    end
end
% Rigid registration
template = clouds{reference}.Location;
parfor i = 1:n
    [~,reg1{i},trans] = procrustes(template,reg1{i},'scaling',false,'reflection',false);
end
% Average shape
tmp = cat(3,reg1{:});
mean_shape = mean(tmp,3);

%% Nonrigid registration #2 (to mean shape)
reg2 = cell(n,1);
parfor i = 1:n
    disp(['Non-rigid registration 2 (',num2str(i),'/',num2str(n),')']);
    reg2{i} = nonrigidICPv1(clouds{i}.Location,mean_shape,conn{i},conn{reference},10,0);
end

%% Rigid registration
vis = cell(n,1);
for i = 1:n
    [~,reg3{i},trans] = procrustes(mean_shape,reg2{i},'reflection',false);
    tmp = trans.b*(clouds{i}.Location)*trans.T + repmat(trans.c(1,:),length(clouds{i}.Location),1);
    vis{i} = pointCloud(tmp);
end
tmp = cat(3,reg3{:});
mean_shape = mean(tmp,3);

%% Distances
disp('Distance calculation')
cdistances = zeros(length(mean_shape),n);
distances = zeros(length(mean_shape),n);

% Check which reference points are inside target
% tri = alphaShape(mean_shape); % Partition into tetrahedrons
% for i = 1:n % Corresponding vertices
%     IsInside = inShape(tri,reg3{i}(:,1),reg3{i}(:,2),reg3{i}(:,3));
%     for j = 1:length(mean_shape)
%         if IsInside(j)
%             cdistances(j,i) = -(norm(mean_shape(j,:)-reg3{i}(j,:)));
%         else
%             cdistances(j,i) = norm(mean_shape(j,:)-reg3{i}(j,:));
%         end
%     end
% end
parfor i = 1:n % Closest vertices
    distances(:,i) = distfield(mean_shape,vis{i}.Location);
end

%% Mean shape
avg_dist = mean(distances,2);
dist = pointCloud(mean_shape);
dist.Intensity = avg_dist;
figure('Name','Mean shape')
pcshow(dist)
axis equal

%% Distances
% tmp = pointCloud(mean_shape);
% for i = 1:n
%     figure
%     tmp.Intensity = distances(:,i);
%     pcshow(tmp);
%     axis equal
% end

%% Grouping
groups = zeros(n,7);
cats = categorical(groups);
tmp = discretize(cell2mat(data(:,4)),2:2:10,'categorical',{'[2-4]','[4-6]','[6-8]','[8-10]'}); % Age
cats(:,1) = categorical(tmp,'Ordinal',false);
cats(:,2) = categorical(data(:,5)); % Gender
cats(:,3) = categorical(data(:,6)); % CP type
cats(:,5) = categorical(data(:,10)); % Leg
% cats(:,7) = categorical(data(:,16)); % Distortion

groups(:,1) = double(cats(:,1))-1; % Age
groups(:,2) = double(categorical(data(:,5))); % Gender
groups(:,3) = double(categorical(data(:,6))); % CP type
groups(:,5) = double(categorical(data(:,10))); % Leg
groups(:,7) = cell2mat(data(:,16))+1; % Distortion
tmp = groups(:,3)-3;
groups(tmp < 0,6) = -1;
groups(:,6) = groups(:,6)+2; % CP & TD
tmp = discretize(groups(:,6),1:0.5:2,'categorical',{'CP','TD'}); % CP & TD
cats(:,6) = categorical(tmp,'Ordinal',false);
varnames = {'Age','Gender','Diagnosis','GMFCS','Leg','Subject group','Distortion'};

ncat = max(groups);
ngroups = size(groups,2);

% Exclusions
groups(:,[2 3 4 5]) = [];
varnames([2 3 4 5]) = [];
ngroups = size(groups,2);

% g1 = data{1,4}; % Age
% g2 = data{1,6}; % CP Type
% g3 = data{1,16}; % Distortion
% for i = 1:n-1
%     g1 = cat(2,g1,data{i+1,4});
%     g2 = char(g2,data{i+1,6});
%     g3 = cat(2,g3,data{i+1,16});
% end
% groups = {g1';g2;g3'};
%% ANOVA
varnames = {'Age','CP Type','Distortion'};
ngroups = size(groups,1);
etas = zeros(ngroups,length(mean_shape));
tb1 = cell(6,7,length(mean_shape));
p = zeros(3,length(mean_shape));
for i = 1:length(mean_shape)
    [p(:,i),tb1(:,:,i)] = anovan(distances(i,:),groups,'display','off','continuous',1);
    for j = 1:ngroups
        SS(j) = tb1{1+j,2,i};
    end
    SSE = tb1{2 + ngroups,2,i};
    etas(:,i) = SS./(SS+SSE);
end
p = 1./p;

%% Figures
% p-values
for i = 1:size(p,1)
    color = p(i,:)';
    figure('WindowState','maximized')
    trisurf(conn{reference},-mean_shape(:,1),mean_shape(:,2), ...
        -mean_shape(:,3),color,'Edgecolor','none');
    title(varnames{i});
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light
    lighting phong;
    colormap('jet')
    set(gca,'ColorScale','log')
%     caxis([0.8 1])
end
saveas(gcf,strcat(string(strength),'_',string(sigma),'_3.png'));
% Distortion from the back
color = p(i,:)';
figure('WindowState','maximized')
trisurf(conn{reference},mean_shape(:,1),-mean_shape(:,2), ...
    -mean_shape(:,3),color,'Edgecolor','none');
title('Back');
set(gca,'XColor','none','YColor','none','ZColor','none')
set(gca,'View',[80.1,0.19])
grid off
axis equal
set(gcf,'Color',[1 1 1])
light('Position',[1 0 0])
lighting phong; 
colormap('jet')
set(gca,'ColorScale','log')
% caxis([0.8 1])
saveas(gcf,strcat(string(strength),'_',string(sigma),'_4.png'));

% etas
for i = 1:size(etas,1)
    color = etas(i,:)';
    figure('WindowState','maximized')
    trisurf(conn{reference},-mean_shape(:,1),mean_shape(:,2), ...
        -mean_shape(:,3),color,'Edgecolor','none');
    title(varnames{i});
    set(gca,'XColor','none','YColor','none','ZColor','none')
    set(gca,'View',[80.1,0.19])
    grid off
    axis equal
    set(gcf,'Color',[1 1 1])
    light
    lighting phong;
    colormap('jet')
    caxis([0 0.2])
end
saveas(gcf,strcat(string(strength),'_',string(sigma),'_5.png'));
% Distortion from the back
color = etas(i,:)';
figure('WindowState','maximized')
trisurf(conn{reference},mean_shape(:,1),-mean_shape(:,2), ...
    -mean_shape(:,3),color,'Edgecolor','none');
title('Back');
set(gca,'XColor','none','YColor','none','ZColor','none')
set(gca,'View',[80.1,0.19])
grid off
axis equal
set(gcf,'Color',[1 1 1])
light('Position',[1 0 0])
lighting phong; 
colormap('jet')
caxis([0 0.2])
saveas(gcf,strcat(string(strength),'_',string(sigma),'_6.png'));

%% SSM
% Training data
Xdata = reg3{1}(:,1);
Ydata = reg3{1}(:,2);
Zdata = reg3{1}(:,3);
for i = 2:n
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
    trisurf(conn{reference},-mean_shape(:,1),mean_shape(:,2), ...
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
    trisurf(conn{reference},-shape(:,1),shape(:,2)+offset, ...
        -shape(:,3),'Edgecolor','none'); % 2 SD
    shape = MEAN-SD*ssmV(:,i);
    shape = reshape(shape,[],3);
    trisurf(conn{reference},-shape(:,1),shape(:,2)-offset, ...
        -shape(:,3),'Edgecolor','none'); % -2 SD
%     saveas(gcf,strcat('PC ',string(i),'.png'));
end

figure
plot((Eval(1:nmax)./sum(Eval))*100);
set(gcf,'Color',[1 1 1])
title('Variation explained by PCs');
ylabel('% of total');
xlabel('PC number');
% saveas(gcf,'PCcum.png');

figure
plot(PCcum(1:50)*100);
set(gcf,'Color',[1 1 1])
title('Cumulative variation explained by PCs');
ylabel('% of total');
xlabel('PC number');

%% SSM ANOVA
p = zeros(ngroups,nmax);
for i = 1:nmax
    [p(:,i),tb] = anovan(Modes(:,i),groups,'display','off');
end
p = 1-p;

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
%%
save(strcat(string(strength),'_',string(sigma)));
%%
beep