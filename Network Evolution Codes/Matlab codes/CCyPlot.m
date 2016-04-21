close all
clear all
clc

data = load('DataSet0\Cum2300Frac0.txt'); % using MM1
data1 = load('DataSet0\Final\Cum2300Frac0.txt'); % using MM3
Exptdata = load('DataSet0\ExptSyn2300.txt'); % Objective
[n iter] = size(data);

count = 2;
while ( data(n,count)>0 )
    count = count + 1;
end

Nmod = 6;
figure;
hold on
plot(data(:,1),Exptdata,'r.-','LineWidth',2)
for i = 2:count-1
    if ( i==2 )
        plot(data(:,1),data(:,i),'r.-','LineWidth',2) 
    elseif ( i==count-1 )
        plot(data(:,1),data(:,i),'k.-','LineWidth',3)
    end
end
plot(data1(:,1),data1(:,2),'m.--','LineWidth',2)
for i = 2:count-1
    if (mod(i,Nmod)==0)
        plot(data(:,1),data(:,i),'b.--','LineWidth',2) 
        if (i==Nmod) legend('Objective','Initial','Best evolved (MM1)','Best evolved (MM3)','Intermediate',2)
        end
    end
end
xlabel('Fraction of times a neuron bursts with population burst')
ylabel('(Cummulative) fraction of total neurons')
grid on
box on


