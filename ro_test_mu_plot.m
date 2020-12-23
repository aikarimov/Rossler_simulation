load('R1_data.mat');
dataCell = cell(8,1); %cell array for data

dataCell{1,1} = R1_5k;
dataCell{2,1} = R1_8k;
dataCell{3,1} = R1_9k;
dataCell{4,1} = R1_10k;
dataCell{5,1} = R1_10k5;
dataCell{6,1} = R1_11k;
dataCell{7,1} = R1_12k8;
dataCell{8,1} = R1_20k;

R1Xspan = [5;8;9;10;10.5;11;12.8;20]; %values for R1, kOhms
R1span = [5:0.1:20]; %for smooth plot
muspan = 10./R1span; %exact mu values
muvalspan = 10./R1Xspan; %exact mu values coordiunated with R1Xspan
muXspan = zeros(8,1); %experimental mu values
Nstart = 1000; %transiend time, samples

Vmu = zeros(8,1); %volume of error attractor
Vn = zeros(8,1);
for k = 1:8
    
    data = dataCell{k,1}; %obtain time series
    xd = data.x; yd = data.y; zd = data.z;
    Yd = 4*[xd';yd';zd']; %get Yd - array of master system dynamics
    N = length(xd); %length of data series
    h = 5.18e-2; %discretization time, sec
    
    %slave with mu
    U = ones(7,1); %state variables for system with mu
    Ys = zeros(3,N); % array of slave system dynamics
    Ys(:,1) = U(1:3);
    Ps = zeros(4,N);
    Ps(:,1) = U(4:7);
    
    %lave without mu
    Un = ones(6,1); %state variables for system without mu
    Ysn = zeros(3,N); % array of slave system dynamics
    Ysn(:,1) = U(1:3);
    
    % %Heun's method
    for i = 2:N
        
        %system with mu
        U1 = U + h*ro_mu(U,Yd(:,i-1));
        U = U + 0.5*h*(ro_mu(U,Yd(:,i-1)) + ro_mu(U1,Yd(:,i)));
        Ys(:,i) = U(1:3);
        Ps(:,i) = U(4:7);
        
        %system without mu
        U1n = Un + h*ro(Un,Yd(:,i-1));
        Un = Un + 0.5*h*(ro(Un,Yd(:,i-1)) + ro(U1n,Yd(:,i)));
        Ysn(:,i) = Un(1:3);       
    end
    
    %calculate statistics on mu
    muM = median(Ps(4,Nstart:end));
    muXspan(k) = muM;
    
    %calculate error attractors
    emu = Ys - Yd; emu = emu(:,Nstart:end);
    en = Ysn - Yd; en = en(:,Nstart:end);
    
    emu_min = min(emu,[],2);  emu_max = max(emu,[],2);
    en_min = min(en,[],2);  en_max = max(en,[],2);
    
    Vmu(k) = sum(emu_max - emu_min);
    Vn(k) = sum(en_max - en_min);
    
    %figure with error attractors
    figure(10);
    subplot(2,4,k);    
    plot3(en(1,:),en(2,:),en(3,:));
    xlabel('$e_x$','interpreter','latex','FontSize',10);
    ylabel('$e_y$','interpreter','latex','FontSize',10);
    zlabel('$e_z$','interpreter','latex','FontSize',10);
    title(['$\mu =',num2str(muvalspan(k)),'$'],'interpreter','latex','FontSize',10);
    set(gca,'TickLabelInterpreter','latex','FontSize',10);
    
    hold on 
    plot3(emu(1,:),emu(2,:),emu(3,:));
    xlabel('$e_x$','interpreter','latex','FontSize',10);
    ylabel('$e_y$','interpreter','latex','FontSize',10);
    zlabel('$e_z$','interpreter','latex','FontSize',10);
    title(['$\mu =',num2str(muvalspan(k)),'$'],'interpreter','latex','FontSize',10);
    set(gca,'TickLabelInterpreter','latex','FontSize',10);
    
    axis square
    grid;
    
    hold off
end

%figure 1 is dedicated to system with mu
figure(1);
subplot(1,2,1);
plot(R1span,R1span,R1Xspan,10./muXspan,'xr','MarkerSize',8);
legend('theoretical','experimental','interpreter','latex');
xlabel('$R_1$, kOhms','interpreter','latex');
ylabel('$\hat{R_1}$, kOhms','interpreter','latex');
grid;
set(gca,'TickLabelInterpreter','latex');
title('\textbf{(a)}','interpreter','latex');

figure(1);
subplot(1,2,2);
semilogy(R1Xspan,abs((R1Xspan - 10./muXspan)./R1Xspan),'xr','MarkerSize',8);
xlabel('$R_1$, kOhms','interpreter','latex');
ylabel('$|\frac{\Delta R_1}{R_1}|$','interpreter','latex');
grid;
set(gca,'TickLabelInterpreter','latex');
title('\textbf{(b)}','interpreter','latex');

%figure 2 is dedicated to both slave systems
figure(2);
semilogy(muvalspan,Vn,'-*',muvalspan,Vmu,'-x','MarkerSize',8);
grid;
set(gca,'TickLabelInterpreter','latex');
xlabel('$\mu$','interpreter','latex');
ylabel('$V(\mathcal{A}_e)$','interpreter','latex');
[~, hobj, ~, ~] = legend('without controller on $\hat{\mu}$','with controller on $\hat{\mu}$','interpreter','latex');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',0.5);