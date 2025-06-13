function plotModes(file,Phi,freq,imodes)

% Load data
switch file
    case 'shell'
        load('shell.mat','xn','Tn','Tc');
    case 'wing'
        load('wing.mat','xn','Tn_wb','Tn_rb','Tn_sk','Tc');
        Tn = [Tn_wb;Tn_rb;Tn_sk];
end

% Precompute
scale = 0.25;
x0 = xn(:,1);
y0 = xn(:,2);
z0 = xn(:,3);
Phi_x = Phi(1:6:end,:);
Phi_y = Phi(2:6:end,:);
Phi_z = Phi(3:6:end,:);
Phi_max = max(abs([Phi_x;Phi_y;Phi_z]),[],1);

nm = length(imodes);
nf = ceil(nm/2);
for j = 1:nf
    figure
    for k = 1:2
        i = 2*(j-1)+k;
        if i<=nm
            I = imodes(i);
            x = x0+scale*Phi_x(:,I)/Phi_max(I);
            y = y0+scale*Phi_y(:,I)/Phi_max(I);
            z = z0+scale*Phi_z(:,I)/Phi_max(I);
            subplot(1,2,k)
            hold on
            patch(x0(Tc)',y0(Tc)',z0(Tc)',ones(size(Tc))','facecolor','none','edgecolor',0.5*[1,1,1]);
            patch(x(Tc)',y(Tc)',z(Tc)',ones(size(Tc))','facecolor','none','edgecolor','k');
            patch(x(Tn)',y(Tn)',z(Tn)',ones(size(Tn))','EdgeColor','none','FaceColor',0.5*[1,1,1]);
            view(40,20);
            set(gca,'color','none','xcolor','none','ycolor','none','zcolor','none');
            title(sprintf('f_{%i} = %.3f Hz',I,freq(I)));
            axis equal;
            axis tight;
            axis vis3d;
        end
    end
end