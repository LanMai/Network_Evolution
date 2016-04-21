clear all
close all
clc   

% NeuronProp
data0 = load('DataSet0\NeuronProp82A2300.txt');
% Isolated Neuron Activity
dataN = load('IsoNeuronActivity\gt0p10\neuronData3104.txt');
% Evolved Network Activity
dataM = load('DataSet0\Final\neuronData2317.txt');
% Initial Connectivity Matrix
InitialConMat0 = load('DataSet0\Initial2300ConMat0_0.txt');
% Evolved Connectivity Matrix
conMatMutated0 = load('DataSet0\Final\Final2300ConMat0.txt');


NoNe = 80;
data = data0(1:NoNe,:);
[m n] = size(data);
Zo = zeros(4,1);
Zdo = zeros(m,1);


L = size(dataN,1);
BurstingNeurons = zeros(m,2);
for i = 1:m
    for j = 1:L
        if (BurstingNeurons(i,1)==0 && i-1==dataN(j,2))
            BurstingNeurons(i,1) = 1;
            BurstingNeurons(i,2) = j;
        end
    end
end


L = size(dataM,1);
SilentNeurons = zeros(m,2);
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

MeanSyn = 0.5*(MaxSyn + MinSyn);

% MeanSynOee = 0.5*(MaxSyn2EE+MinSyn2EE);
% MeanSynIee = 0.5*(MaxSyn3EE+MinSyn3EE);
CF = 15; % Contraction factor for colormap (<=32)

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

figure;
plot3(x',y',z','-','LineWidth',2)
hold on
plot3(xou',you',zo','r-','LineWidth',2)
plot3(xob',yob',zo','r-','LineWidth',2)
Cmin = 1*10;
Cmax = 60;
for i = 1:NoNe
    SynCount = (outSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.1,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EE, initially)')
axis([0,5,0,5])
axis equal
grid on
% Lo1 = num2str((MaxSyn-MinSyn)*10/64);
% Lo2 = num2str((MaxSyn-MinSyn)*20/64);
% Lo3 = num2str((MaxSyn-MinSyn)*30/64);
% Lo4 = num2str((MaxSyn-MinSyn)*40/64);
% Lo5 = num2str((MaxSyn-MinSyn)*50/64);
% Lo6 = num2str((MaxSyn-MinSyn)*60/64);
Fac = (MaxSyn-MinSyn)/(Cmax-Cmin);
Lo1 = num2str(Fac*(10-Cmin)+MinSyn);
Lo2 = num2str(Fac*(20-Cmin)+MinSyn);
Lo3 = num2str(Fac*(30-Cmin)+MinSyn);
Lo4 = num2str(Fac*(40-Cmin)+MinSyn);
Lo5 = num2str(Fac*(50-Cmin)+MinSyn);
Lo6 = num2str(Fac*(60-Cmin)+MinSyn);
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on

figure;
plot3(x',y',z','-','LineWidth',2)
hold on
plot3(xou',you',zo','r-','LineWidth',2)
plot3(xob',yob',zo','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (inSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.1,'MarkerSize',12,'MarkerFaceColor',[R G B])
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

figure;
plot3(x',y',z','-','LineWidth',2)
hold on
plot3(xou',you',zo','r-','LineWidth',2)
plot3(xob',yob',zo','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (outSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (BurstingNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.1,'MarkerSize',12,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
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

figure;
plot3(x',y',z','-','LineWidth',2)
hold on
plot3(xou',you',zo','r-','LineWidth',2)
plot3(xob',yob',zo','r-','LineWidth',2)
for i = 1:NoNe
   SynCount = (inSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
   [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
   if (BurstingNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.1,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
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

%% In-Syn (Final)
FigHandle = figure('Position', [100, 100, 1049, 895]);
subplot(2,2,1)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (outSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EE, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%---------------------------------------------
subplot(2,2,2)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (inSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
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
    if (BurstingNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
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
   if (BurstingNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'wo','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
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


break

XoutCentroid = data(:,1)'*outSynFinalEE(:,1)/sum(outSynFinalEE);
YoutCentroid = data(:,2)'*outSynFinalEE(:,1)/sum(outSynFinalEE);

figure;
plot3(x',y',z','-','LineWidth',2)
hold on
plot3(xou',you',zo','r-','LineWidth',2)
plot3(xob',yob',zo','r-','LineWidth',2)
for i = 1:NoNe
    [R, G, B] = jetplot(outSynFinalEE(i,1)/MeanSynOee*CF);
    if (BurstingNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),outSynFinalEE(i,1)/SynMaxEE*5,'ko','LineWidth',0.1,'MarkerSize',15,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i,1),data(i,2),outSynFinalEE(i,1)/SynMaxEE*5,'wo','LineWidth',3,'MarkerSize',15,'MarkerFaceColor',[R G B])
        plot3(data(i,1),data(i,2),outSynFinalEE(i,1)/SynMaxEE*5,'ko','LineWidth',0.1,'MarkerSize',15)
    end
    
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),outSynFinalEE(i,1)/SynMaxEE*5,'ko','LineWidth',2,'MarkerSize',15,'MarkerFaceColor',[R G B])
    end
end
plot3(XoutCentroid,YoutCentroid,5,'ko','LineWidth',3,'MarkerSize',15,'MarkerFaceColor',[0 1 0])
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('No. of out-going synapses (excitatory to excitatory neurons, finally)')
axis([0,6,0,6])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on

break


neuronProp0 = load('I:\Hyogo\July2014\pmProp.txt');
neuronProp = flipud(neuronProp0);
figure;
mesh(neuronProp0)
figure;
mesh(neuronProp)

xa = 0:0.1:5.0;
ya = 0:0.1:5.0;
[Xa,Ya] = meshgrid(xa,ya);
gt_prop = interp2(Xa,Ya,neuronProp,data(:,1),data(:,2));
figure;
plot3(x',y',z','-','LineWidth',2)
hold on
plot3(xou',you',zo','r-','LineWidth',2)
plot3(xob',yob',zo','r-','LineWidth',2)
plot3(data(:,1),data(:,2),gt_prop,'o')
axis([0,5,0,5])
axis equal
grid on
box on

figure;
plot(gt_prop,outSynFinalEE,'bo','LineWidth',2,'MarkerSize',8)
grid on
box on

figure;
plot(gt_prop,inSynFinalEE,'ro','LineWidth',2,'MarkerSize',8)
grid on
box on
