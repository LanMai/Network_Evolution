clear all
close all
clc   

% NeuronProp
data = load('DataSet0\NeuronProp82A2300.txt');
% Isolated Neuron Activity
dataN = load('IsoNeuronActivity\gt0p10\neuronData3104.txt');
% Evolved Network Activity
dataM = load('DataSet0\Final\neuronData2317.txt');
% Initial Connectivity Matrix
InitialConMat0 = load('DataSet0\Initial2300ConMat0_0.txt');
% Evolved Connectivity Matrix
conMatMutated0 = load('DataSet0\Final\Final2300ConMat0.txt');

[m n] = size(data);
NoNe = 80;
Zo = zeros(4,1);
Zdo = zeros(m,1);

L = size(dataN,1);
BurstingNeurons = zeros(m,2); % in iso-popAct
for i = 1:m
    for j = 1:L
        if (BurstingNeurons(i,1)==0 && i-1==dataN(j,2))
            BurstingNeurons(i,1) = 1;
            BurstingNeurons(i,2) = j;
        end
    end
end

L = size(dataM,1);
SilentNeurons = zeros(m,2); % bursting Neurons in Evolved network
for i = 1:m
    for j = 1:L
        if (SilentNeurons(i,1)==0 && i-1==dataM(j,2))
            SilentNeurons(i,1) = 1;
            SilentNeurons(i,2) = j;
        end
    end
end

InitialConMatEE = InitialConMat0(1:NoNe,1:NoNe);
conMatMutatedEE = conMatMutated0(1:NoNe,1:NoNe);

fracEE0 = sum(sum(InitialConMatEE))./(NoNe*NoNe);
fracEE1 = sum(sum(conMatMutatedEE))./(NoNe*NoNe);

[fracEE0 fracEE1]

inSynInitialEE = sum(InitialConMatEE,1)';
outSynInitialEE = sum(InitialConMatEE,2);
inSynFinalEE = sum(conMatMutatedEE,1)';
outSynFinalEE = sum(conMatMutatedEE,2);

MaxOutSynEE = max(outSynFinalEE);
MaxInSynEE = max(inSynFinalEE);

MaxSyn = max([MaxOutSynEE,MaxInSynEE]);

MinOutSynEE = min(outSynFinalEE);
MinInSynEE = min(inSynFinalEE);

MinSyn = min([MinOutSynEE,MinInSynEE]);

% MeanSyn = 0.5*(MaxSyn + MinSyn);


[inSynInitialEEsort,inSynInitialEEsortID] = sort(inSynInitialEE,'descend');
[outSynInitialEEsort,outSynInitialEEsortID] = sort(outSynInitialEE,'descend');
[inSynFinalEEsort,inSynFinalEEsortID] = sort(inSynFinalEE,'descend');
[outSynFinalEEsort,outSynFinalEEsortID] = sort(outSynFinalEE,'descend');

Nmax = 3;

% NPM
x = [0.2 5.0];
y = [0.1 3.4];
z = [0.0 0.0];

% PM
xou = [0.2; 1.3];
you = [0.2; 5.0];
xob = [0.2; 0.2];
yob = [0.2; 0.1];
zo = [0.0 0.0];

%% Out-Syn (Initial)
FigHandle = figure('Position', [100, 100, 1049, 895]);
subplot(2,2,1)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
Cmin = 1*10;
Cmax = 60;
maxSynCount = 0;
for i = 1:NoNe
    SynCount = (outSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    if (SynCount>maxSynCount)
        maxSynCount = SynCount;
    end
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.2,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.1,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
for i = 1:Nmax
    t = outSynInitialEEsortID(i);
    for j = 1:NoNe        
        if ( InitialConMatEE(t,j)>0 )
            SynCount1 = (outSynInitialEE(t,1)-MinSyn)/(MaxSyn-MinSyn);
            SynCount2 = (outSynInitialEE(j,1)-MinSyn)/(MaxSyn-MinSyn);
%             plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[SynCount1 SynCount2],'r-','LineWidth',0.5)
            plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[maxSynCount maxSynCount],'r-','LineWidth',0.5)
        end
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EE, initially)')
axis([0,5,0,5])
axis equal
grid on
Fac = (MaxSyn-MinSyn)/(Cmax-Cmin);
Lo1 = num2str(Fac*(10-Cmin)+MinSyn);
Lo2 = num2str(Fac*(20-Cmin)+MinSyn);
Lo3 = num2str(Fac*(30-Cmin)+MinSyn);
Lo4 = num2str(Fac*(40-Cmin)+MinSyn);
Lo5 = num2str(Fac*(50-Cmin)+MinSyn);
Lo6 = num2str(Fac*(60-Cmin)+MinSyn);
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%------------------------------------------
subplot(2,2,2)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
maxSynCount = 0;
for i = 1:NoNe
    SynCount = (inSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    if (SynCount>maxSynCount)
        maxSynCount = SynCount;
    end
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
for i = 1:Nmax
    t = inSynInitialEEsortID(i);
    for j = 1:NoNe        
        if ( InitialConMatEE(j,t)>0 )
            SynCount1 = (inSynInitialEE(t,1)-MinSyn)/(MaxSyn-MinSyn);
            SynCount2 = (inSynInitialEE(j,1)-MinSyn)/(MaxSyn-MinSyn);
%             plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[SynCount1 SynCount2],'r-','LineWidth',0.5)
            plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[maxSynCount maxSynCount],'r-','LineWidth',0.5)
        end
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (EE, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%---------------------------------------------
subplot(2,2,3)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
maxSynCount = 0;
for i = 1:NoNe
    SynCount = (outSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    if (SynCount>maxSynCount)
        maxSynCount = SynCount;
    end
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
%     if (BurstingNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])       
%     else 
%         plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
%     end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
for i = 1:Nmax
    t = outSynFinalEEsortID(i);
    for j = 1:NoNe        
        if ( conMatMutatedEE(t,j)>0 )
            SynCount1 = (outSynFinalEE(t,1)-MinSyn)/(MaxSyn-MinSyn);
            SynCount2 = (outSynFinalEE(j,1)-MinSyn)/(MaxSyn-MinSyn);
%             plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[SynCount1 SynCount2],'r-','LineWidth',0.5)
            plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[maxSynCount maxSynCount],'r-','LineWidth',0.5)
        end
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EE, finally)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%---------------------------------------------
subplot(2,2,4)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
maxSynCount = 0;
for i = 1:NoNe
   SynCount = (inSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
   if (SynCount>maxSynCount)
        maxSynCount = SynCount;
    end
   [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
%    if (BurstingNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
%     else
%         plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
%     end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
for i = 1:Nmax
    t = inSynFinalEEsortID(i);
    for j = 1:NoNe        
        if ( conMatMutatedEE(j,t)>0 )
            SynCount1 = (inSynFinalEE(t,1)-MinSyn)/(MaxSyn-MinSyn);
            SynCount2 = (inSynFinalEE(j,1)-MinSyn)/(MaxSyn-MinSyn);
%             plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[SynCount1 SynCount2],'r-','LineWidth',0.5)
            plot3([data(t,1) data(j,1)],[data(t,2) data(j,2)],5*[maxSynCount maxSynCount],'r-','LineWidth',0.5)
        end
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (EE, finally)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on


%% Out-Syn (Initial)
FigHandle = figure('Position', [100, 100, 1049, 895]);
subplot(2,2,1)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
Cmin = 1*10;
Cmax = 60;
for i = 1:NoNe
    SynCount = (outSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);    
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.2,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.1,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EE, initially)')
axis([0,5,0,5])
axis equal
grid on
Fac = (MaxSyn-MinSyn)/(Cmax-Cmin);
Lo1 = num2str(Fac*(10-Cmin)+MinSyn);
Lo2 = num2str(Fac*(20-Cmin)+MinSyn);
Lo3 = num2str(Fac*(30-Cmin)+MinSyn);
Lo4 = num2str(Fac*(40-Cmin)+MinSyn);
Lo5 = num2str(Fac*(50-Cmin)+MinSyn);
Lo6 = num2str(Fac*(60-Cmin)+MinSyn);
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%------------------------------------------
subplot(2,2,2)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (inSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (EE, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%---------------------------------------------
subplot(2,2,3)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (outSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
%     if (BurstingNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])       
%     else 
%         plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
%     end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EE, finally)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%---------------------------------------------
subplot(2,2,4)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
   SynCount = (inSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
   [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
%    if (BurstingNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
%     else
%         plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
%     end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',8,'MarkerFaceColor',[R G B])
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (EE, finally)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
