function valueode;
% estimate error of one run
% n=10;
% x=0.3e-2;
% y =5e-6;
% p=5e11;
% q= 0;


n=6;
x=18e-4;
y =5e-6;
z=10e3;
zz=5e1;
p=1e5;
q= 5e1; 
prob=90e-6



A_1=0.50;
F_12=0.1;

Y0=zeros(1,n); 
Y0(1)=F_12;
Y0(n+1)=A_1;
theta=[x,y,z,zz,p,q]; 

t_range=linspace(0,337,337); 
[t_val,Y_val]=ode23s(@ode_fibril_LFAO_New7,t_range,Y0,[],n,theta);
Y_val([1:20:337 ],[1 3  n n+1])

signal=Y_val(:,n)*0;
for i=n
signal=signal + Y_val(:,i);
% signal= Y_val(:,i);
end

signal(1:20:337);
smax=max(signal);
signal=signal-smax/5;
signal(signal<0)=0;
signal(1:20:337);

signal2=signal;
signal2 = (signal2 - min(signal2))/(max(signal2) - min(signal2));
signal2(signal2<0)=0;
hold on
plot(t_range,signal2)

signal1=Y_val(:,n)*0;
 for i=2:n-1
     signal1=signal1 +(Y_val(:,i)-smax/5)*prob  ;
 end
 signal1(1:20:337);

signal3=signal1; 
signal3 = (signal3 - min(signal3))/(max(signal3) - min(signal3));
signal3(signal3<0)=0;
hold on
plot(t_range,signal3)
 
signal=signal+signal1;
signal = (signal - min(signal))/(max(signal) - min(signal));
signal(signal<0)=0;
hold on
plot(t_range,signal)

A = [t_range', signal];
fileID = fopen('N_6_2_data.txt','w');
fprintf(fileID,'%6s %12s\n','Time','Signal');
fprintf(fileID,'%6.2f %12.8f\n',A');
fclose(fileID);
hold on

X=[0
0.07652
0.11945
0.19411
0.24835
0.55645
1
0.87959
0.94942];
Y=signal([1,24,48,96,144,168,192,216,336]);
axis=[1,24,48,96,144,168,192,216,336];
hold on
plot(axis,X)
mdl = fitlm(Y,X)
                

