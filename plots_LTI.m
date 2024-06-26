%% Plot Figures 5 and 6
clear all; clc

addpath('Functions')

%% Surface plot (Figure 6)
% Load data
load 'Data/20240131 - LTI example'

cap = 3;
Hnorm(Hnorm>cap) = cap;
h = figure;
h.Position=[100 100 600 300];
s = surf(GG1,GG2,(Hnorm'));
xlabel('$G_{(1)}$')
ylabel('$G_{(2)}$')
zlabel('$\gamma$')
s.EdgeColor = 'none';
set(gca,'xscale','log')

c = ['mrcyg'];

hold all

final_idx = [149 165 144 196 210];
for k = 1:5
    load(['Data/20240130 - Data LTI - Init' num2str(k)])
    GG = GG{1};
    GG = cell2mat(GG);
    gamma = cell2mat(gamma);
    GG = GG(:,1:final_idx(k));
    gamma = gamma(1:final_idx(k));
    plot_handle(k) = plot3(GG(1,:),GG(2,:),(gamma)+1,c(k));
    plot3(GG(1,1),GG(2,1),(gamma(1))+1,[c(k) 's']);
    plot3(GG(1,end),GG(2,end),(gamma(end))+1,[c(k) 'o']);
end
xlim([1e-4 100])
ylim([-1 1])
view(0,90)
a = colorbar;
a.Label.String = '\gamma';

set(gca,'fontsize', 14)
set(findall(gcf,'type','line'),'linewidth',1)

load 'Data/20240131 - LTI example'
plot_handle(end+1) = plot3(min_G_psw(1),min_G_psw(2),100,'x','color',[255 155 0]/255,'markersize',1,'linewidth',15);

legend(plot_handle,'Init I','Init II','Init III','Init IV','Init V','Optimal')

%% Iteration history (Figure 5) - Run previous section first
clear plot_handle

c = ['mrcyg'];
h = figure;
h.Position = [100 100 600 300];
for k = 1:5
    load(['Data/20240130 - Data LTI - Init' num2str(k)])
    idxs = 0;
    gamma_ = cell2mat(gamma);
    gamma_tmp{k} = gamma_;
    plot_handle(k) = semilogy(0:length(gamma_)-1,gamma_,c(k));
    hold all
    for i = 1:length(gamma)
        idxe = idxs + length(gamma{i}) - 1;
        if i == 1
            semilogy(0,gamma{1},[c(k) 's'])
        else
            semilogy(idxe,gamma{i}(end),[c(k) 'x'])
        end
        idxs = idxe + 1;
    end
end
xlim([0 70])
xlabel('Iteration index')
ylabel('$\gamma$')
set(gca,'fontsize', 14)

axes('Position',[0.3 0.5 0.35 0.35])
for k = 1:5
    load(['Data/20240130 - Data LTI - Init' num2str(k)])
    idxs = 0;
    gamma_ = cell2mat(gamma);
    semilogy(0:length(gamma_)-1,gamma_,c(k))
    hold all
    for i = 1:length(gamma)
        idxe = idxs + length(gamma{i}) - 1;
        semilogy(idxe,gamma{i}(end),[c(k) 'x'])
        idxs = idxe + 1;
    end
end
xlim([50 60])
set(gca,'fontsize', 14)
set(findall(gcf,'type','line'),'linewidth',1)
set(findall(gcf,'type','line'),'markersize',8)

legend(plot_handle,'Init I','Init II','Init III','Init IV','Init V')