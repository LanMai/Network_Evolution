clear all
close all
clc   

% NeuronProp
data = load('DataSet1\NeuronProp160A7000.txt'); 
% IsoNeuronActivity
dataN = load('IsoNeuronActivity\gt0p34\neuronData3104.txt'); 
% EvolvedNetworkActivity
dataM = load('DataSet1\Final\neuronData7010.txt');
% InitialConMat
InitialConMat0 = load('DataSet1\Initial7000ConMat0_0.txt');
% EvolvedConMat
conMatMutated0 = load('DataSet1\Final\Final7000ConMat0.txt');

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
InitialConMatEI = InitialConMat0(1:NoNe,NoNe+1:m);
conMatMutatedEI = conMatMutated0(1:NoNe,NoNe+1:m);
InitialConMatIE = InitialConMat0(NoNe+1:m,1:NoNe);
conMatMutatedIE = conMatMutated0(NoNe+1:m,1:NoNe);
InitialConMatII = InitialConMat0(NoNe+1:m,NoNe+1:m);
conMatMutatedII = conMatMutated0(NoNe+1:m,NoNe+1:m);

fracEE0 = sum(sum(InitialConMatEE))./(NoNe*NoNe);
fracEE1 = sum(sum(conMatMutatedEE))./(NoNe*NoNe);
fracEI0 = sum(sum(InitialConMatEI))./(NoNe*(m-NoNe));
fracEI1 = sum(sum(conMatMutatedEI))./(NoNe*(m-NoNe));
fracIE0 = sum(sum(InitialConMatIE))./(NoNe*(m-NoNe));
fracIE1 = sum(sum(conMatMutatedIE))./(NoNe*(m-NoNe));
fracII0 = sum(sum(InitialConMatII))./((m-NoNe)*(m-NoNe));
fracII1 = sum(sum(conMatMutatedII))./((m-NoNe)*(m-NoNe));

[fracEE0 fracEE1
    fracEI0 fracEI1
    fracIE0 fracIE1
    fracII0 fracII1]

inSynInitialEE = sum(InitialConMatEE,1)';
inSynInitialEI = sum(InitialConMatEI,1)';
inSynInitialIE = sum(InitialConMatIE,1)';
inSynInitialII = sum(InitialConMatII,1)';

outSynInitialEE = sum(InitialConMatEE,2);
outSynInitialEI = sum(InitialConMatEI,2);
outSynInitialIE = sum(InitialConMatIE,2);
outSynInitialII = sum(InitialConMatII,2);

inSynFinalEE = sum(conMatMutatedEE,1)';
inSynFinalEI = sum(conMatMutatedEI,1)';
inSynFinalIE = sum(conMatMutatedIE,1)';
inSynFinalII = sum(conMatMutatedII,1)';

outSynFinalEE = sum(conMatMutatedEE,2);
outSynFinalEI = sum(conMatMutatedEI,2);
outSynFinalIE = sum(conMatMutatedIE,2);
outSynFinalII = sum(conMatMutatedII,2);

MaxOutSynEE = max(outSynFinalEE);
MaxOutSynEI = max(outSynFinalEI);
MaxOutSynIE = max(outSynFinalIE);
MaxOutSynII = max(outSynFinalII);
MaxInSynEE = max(inSynFinalEE);
MaxInSynEI = max(inSynFinalEI);
MaxInSynIE = max(inSynFinalIE);
MaxInSynII = max(inSynFinalII);

MaxSyn = max([MaxOutSynEE,MaxOutSynEI,MaxOutSynIE,MaxOutSynII,...
    MaxInSynEE,MaxInSynEI,MaxInSynIE,MaxInSynII]);

MinOutSynEE = min(outSynFinalEE);
MinOutSynEI = min(outSynFinalEI);
MinOutSynIE = min(outSynFinalIE);
MinOutSynII = min(outSynFinalII);
MinInSynEE = min(inSynFinalEE);
MinInSynEI = min(inSynFinalEI);
MinInSynIE = min(inSynFinalIE);
MinInSynII = min(inSynFinalII);

MinSyn = min([MinOutSynEE,MinOutSynEI,MinOutSynIE,MinOutSynII,...
    MinInSynEE,MinInSynEI,MinInSynIE,MinInSynII]);

MeanSyn = 0.5*(MaxSyn + MinSyn);


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
for i = 1:NoNe
    SynCount = (outSynInitialEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.2,'MarkerSize',12,'MarkerFaceColor',[R G B])
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
% break
%------------------------------------------
subplot(2,2,2)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (outSynInitialEI(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EI, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%------------------------------------------
subplot(2,2,3)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:m-NoNe
    SynCount = (outSynInitialIE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (IE, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
%------------------------------------------
subplot(2,2,4)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:m-NoNe
    SynCount = (outSynInitialII(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (II, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on

%% In-Syn (Initial)
FigHandle = figure('Position', [100, 100, 1049, 895]);
subplot(2,2,1)
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
subplot(2,2,2)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:m-NoNe
    SynCount = (inSynInitialEI(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (EI, initially)')
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
    SynCount = (inSynInitialIE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (IE, initially)')
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
for i = 1:m-NoNe
    SynCount = (inSynInitialII(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (II, initially)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on


%% Out-Syn (Final)
FigHandle = figure('Position', [100, 100, 1049, 895]);
subplot(2,2,1)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (outSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
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
subplot(2,2,2)
plot(x',y','-','LineWidth',2)
hold on
plot(xou',you','r-','LineWidth',2)
plot(xob',yob','r-','LineWidth',2)
for i = 1:NoNe
    SynCount = (outSynFinalEI(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i,1),data(i,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (EI, finally)')
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
for i = 1:m-NoNe
    SynCount = (outSynFinalIE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
%     plot(data(i+NoNe,1),data(i+NoNe,2),'ko','LineWidth',0.1,'MarkerSize',12,'MarkerFaceColor',[R G B])
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i+NoNe,1)==1)
%         plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('out-going (IE, finally)')
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
for i = 1:m-NoNe
    SynCount = (outSynFinalII(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])       
    else 
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i+NoNe,1)==1)
%         plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('No. of out-going synapses (inhibitory to inhibitory neurons, finally)')
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
   SynCount = (inSynFinalEE(i,1)-MinSyn)/(MaxSyn-MinSyn);
   [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
   if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (EE, finally)')
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
for i = 1:m-NoNe
    SynCount = (inSynFinalEI(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i+NoNe,1)==1)
%         plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-coming (EI, finally)')
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
    SynCount = (inSynFinalIE(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i,1)==1)
        plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i,1),data(i,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i,1)==1)
%         plot3(data(i,1),data(i,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (IE, finally)')
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
for i = 1:m-NoNe
    SynCount = (inSynFinalII(i,1)-MinSyn)/(MaxSyn-MinSyn);
    [R, G, B] = jetplot(Cmin+SynCount*(Cmax-Cmin));
    if (SilentNeurons(i+NoNe,1)==1)
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
    else
        plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ws','LineWidth',2,'MarkerSize',12,'MarkerFaceColor',[R G B])
    end
    
%     if (SilentNeurons(i+NoNe,1)==1)
%         plot3(data(i+NoNe,1),data(i+NoNe,2),5*SynCount,'ko','LineWidth',0.3,'MarkerSize',12,'MarkerFaceColor',[R G B])
%     end
end
xlabel('g_{leak} (nS)')
ylabel('g_{NaP} (nS)')
title('in-comming (II, finally)')
axis([0,5,0,5])
axis equal
grid on
colorbar('YTickLabel',{Lo1,Lo2,Lo3,Lo4,Lo5,Lo6})
box on
