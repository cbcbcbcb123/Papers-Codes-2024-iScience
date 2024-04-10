clear
clc

% R1: Da->integrin “正常态基因转录整合素蛋白”
% R2: YAP+Da->Da* “YAP转录因子结合整合素基因使之处于激活态”
% R3: Da*->Da+YAP “激活态基因变为正常态基因”
% R4: Da*->integrin “激活态基因转录整合素蛋白”
% R5: integrin->null “整合素蛋白降解”
% R6: integrin->YAP “整合素促进YAP入核” NOTE THAT only this change!!!
% R7: YAP->null “YAP蛋白降解”

t = 0;
k = 10; % Fix 10 转录活性提高倍数
alpha = 0.1; % 硬度权重 Fix 0.1
beta = 1; % TGF-beta 权重 Fix 1
K_H = 0.003; % 最大速率 Fix 0.003
K_M = 10; % 希尔常数 Fix 10

% soft (S=-1;T=-1); soft+TGF-b (S=-1;T=0); stiff (S=0;T=-1); stiff+TGF-b (S=0;T=0);
S = 0; % stiffness
T = 0; % TGF-b

r1 = 0.3; % Fix 0.3 低转录活性
r2 = 0.0001; % Fix 0.0001
r3 = 0.001;  % Fix 0.001
r4 = k*r1; % Fix 高转录活性
r5 = 0.1; % Fix 0.1
r7 = 0.1; % Fix 0.1

C_A = 0; % 整合素总浓度
C_A_c = 0;  % 整合素成簇浓度
C_yap_nuc = 0; % YAP入核浓度
C_yap_cyto = 1000; % YAP出核浓度 Fix 1000
C_Da = 1;  % 正常态整合素基因浓度
C_Da_a = 0; % 激活态整合素基因浓度

Ncyc = 5000000;
TimeArry = zeros(Ncyc,1);
C_A_Arry = zeros(Ncyc,1);
C_yap_Arry = zeros(Ncyc,1);
C_Da_Arry = zeros(Ncyc,1);
C_Da_a_Arry = zeros(Ncyc,1);

for i=1:Ncyc
    %======
    % Calculate clustered integrin
    % C_A_c = (10^S)/(10^S+1)*C_A;
    C_A_c = (10^S+C_yap_nuc)/(10^S+C_yap_nuc+1)*C_A;

    % Calculate R6 rate
    r6 = K_H*(alpha*10^S+beta*10^T)/(K_M+alpha*10^S+beta*10^T);

    % Calculate time
    Tlink_R1 = -log(rand)/(r1*C_Da);
    Tlink_R2 = -log(rand)/(r2*C_yap_nuc*C_Da);
    Tlink_R3 = -log(rand)/(r3*C_Da_a);
    Tlink_R4 = -log(rand)/(r4*C_Da_a);
    Tlink_R5 = -log(rand)/(r5*C_A);
    Tlink_R6 = -log(rand)/(r6*C_yap_cyto*C_A_c);
    Tlink_R7 = -log(rand)/(r7*C_yap_nuc);
    %======
    % Min time
    t1 = min(Tlink_R1,Tlink_R2);
    t2 = min(Tlink_R3,Tlink_R4);
    t3 = min(Tlink_R5,Tlink_R6);
    t4 = min(t1,t2);
    t5 = min(t3,t4);
    tmin = min(Tlink_R7,t5);
    %======
    % Go on
    t = t+tmin;
    %======
    % Things happen
    if tmin==Tlink_R1
        C_A = C_A+1;
    elseif tmin==Tlink_R2
        C_Da = C_Da-1;C_Da_a = C_Da_a+1;
    elseif tmin==Tlink_R3
        C_Da = C_Da+1;C_Da_a = C_Da_a-1;
    elseif tmin==Tlink_R4
        C_A = C_A+1;
    elseif tmin==Tlink_R5
        C_A = C_A-1;
    elseif tmin==Tlink_R6
        C_yap_nuc = C_yap_nuc+1;
        C_yap_cyto = C_yap_cyto-1;
    elseif tmin==Tlink_R7
        C_yap_nuc = C_yap_nuc-1;
        C_yap_cyto = C_yap_cyto+1;
    end
    %======
    TimeArry(i,1)=t;
    C_A_Arry(i,1)=C_A_c;
    C_yap_Arry(i,1)=C_yap_nuc;
    C_Da_Arry(i,1)=C_Da;
    C_Da_a_Arry(i,1)=C_Da_a;
end

% figure (1)
% plot(TimeArry,C_A_Arry)

figure (2);
histogram(C_A_Arry,50,'EdgeColor','none','Normalization','probability');
ax = gca;
% set(ax,'Position',[0.3 0.3 0.05 0.05],'FontSize',5)
set(ax,'Position',[0.3 0.3 0.3 0.3],'FontSize',10)
xlabel(ax,'Clustered integrin size');ylabel(ax,'Probability');
pbaspect(ax,[1 1 1])
hold on

exportgraphics(ax,'Peppers300.jpg','Resolution',300);

% save('data.txt','C_A_Arry','-ascii')