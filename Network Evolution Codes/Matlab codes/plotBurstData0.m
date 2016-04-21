close all
clear all
clc

NoN = 160;
NoNe = 80;
popThreshold = 0.30*NoNe;

%% Initial network activity
% Z  = load('DataSet1\BurstData7010.txt');
% Z2 = load('DataSet1\NeuronBurstData7010.txt');

%% Evolved network activity
Z  = load('DataSet1\Final\BurstData7010.txt');
Z2 = load('DataSet1\Final\NeuronBurstData7010.txt');


[m n] = size(Z);
M = 1;
while Z(M,1)>0
    M = M + 1;
end

if M<m
    m = M - 1;
end

figure;
subplot(5,1,[1:3])
hold on
for i=2:NoNe+1
    plot(Z(1:m,1)./1000,(i).*Z(1:m,i)-1,'.','LineWidth',3,'MarkerSize',8)
end
for i=NoNe+2:n-2
    plot(Z(1:m,1)./1000,(i).*Z(1:m,i)-1,'m.','LineWidth',3,'MarkerSize',8)
end
xlim([40 180])
ylim([0 n-1])
% xlabel('Time (s)')
ylabel('Neuron #')
box on
grid on
subplot(5,1,4)
hold on
plot(Z(1:m,1)./1000,Z(1:m,n),'r-','LineWidth',1)
plot(Z(1:m,1)./1000,popThreshold*ones(m,1),'r--','LineWidth',1)
% plot(Z2(:,1)./1000,Z(1:m,n),'r-')
% plot(popThreshold*NoNe*ones(m,1),'r-')
xlim([40 180])
ylim([0 NoNe])
% xlabel('Time (s)')
ylabel('No. of neurons bursting')
box on
grid on
subplot(5,1,5)
plot(Z2(:,1)./1000,Z2(:,2),'r-','LineWidth',2)
xlim([40 180])
ylim([-68 -30])
xlabel('Time (s)')
ylabel('Avg. population activity (mV)')
box on
grid on

figure;
hold on
for i=2:NoNe+1
    plot(Z(1:m,1)./1000,(i).*Z(1:m,i)-1,'.','LineWidth',3,'MarkerSize',8)
end
for i=NoNe+2:n-2
    plot(Z(1:m,1)./1000,(i).*Z(1:m,i)-1,'m.','LineWidth',3,'MarkerSize',8)
end
xlim([40 180])
ylim([0 n-1])
xlabel('Time (s)')
ylabel('Neuron #')
box on
grid on
