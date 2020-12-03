%% Introduction
% Samuel Chernov(sc1823)
% Daniella Chung(djc375)
% Andrew Manfredi(ajm418) 

% Orbital Mechanics: Project 3 
% Professor Mookerjee
% Fall 2020 

% Problem #1

%% Part 0: Initializing constants 
theta=74; % Degress, Ang. btwn E Velo & Pos. Vctr
Vinf=10290; % m/s 
rE=6378e3; % Radius of Earth, m 
mu=3.986e14; % Gravitational Constant 

%% Part A: Epsilon Values

% Initial Calculations 
E=Vinf^2/2; % Energy m^2/s^2
a=-mu/(2*E); % Semimajor axis, m 
vL=0; % initial position angle 
rl=rE; % From Project3Help

% Initial Arrays 
posEps=[]; % Epsilon values
vLNew=[]; % New vl values
vLDif=[]; % Difference of New-Old vl terms 

% While loop 
i=1;
t=1; % Truth value, initialization 
vLOld=vL; % Initialization for old vL
vLTemp=vL; % Initialization for vL

while t==1
    %% Initial Values
    if i~=1
        vLTemp=vLNew(i-1);
    end
    %% Epsilon Solution 
    
    % Coefficients of the Eqn. 
    term1=1;
    term2=(rl/a)*cosd(vLTemp);
    term3=(rl/a)-1;
    
    % Use the 'roots' fn. & save to array  
    tempEqn=[term1 term2 term3];
    tempEps=roots(tempEqn); 
    loc=find(tempEps>0); % Use the 'find' for '+' results
    posEps(i)=tempEps(loc); 
    
    %% vL solution
    tempvL=acosd(-1/posEps(i))-theta;
    vLNew(i)=tempvL;
    
    %% vL Difference
    if i~=1
        vLOld=vLNew(i-1);
    end
    vLDif(i)=vLNew(i)-vLOld; % Difference in vL
    
    %% Breaking Condition
    bCond=abs(vLDif(i)); % Absolute value of the difference
    
    if bCond<0.4 % break when difference <0.4
        break
    end
    
    i=i+1;
end

%% Part B

ce=['Converged Epsilon is: ',num2str(posEps(length(posEps)))];
cv=['Converged vL is: ',num2str(vLNew(length(vLNew))),' °'];
ci=['Converged Iteration Amount: ',num2str(length(posEps))];
disp(ce);
disp(cv);
disp(ci); 

%% Part C
itr=linspace(1,length(posEps),3); % # of iterations
figure 
hold on 
xticks(1:1:length(posEps));
yyaxis left
plot(itr,posEps,'LineWidth',2);
ylabel('Epsilon', 'FontWeight','Bold');
yyaxis right
plot(itr,vLNew,'LineWidth',2);
ylabel('Position Angle(vL) (°)', 'FontWeight','Bold');
title('Epsilion/vL vs. Iteration', 'FontWeight','Bold');
xlabel('Iteration #', 'FontWeight','Bold');
legend('Epsilon', 'vL (°)');

%% Part D 
% Orbit Perigee 
rp=a*(1-posEps(length(posEps)));
rpd=['Perigee Radius is: ',num2str(rp),' m'];
disp(rpd);

%% Part E 
vLnch=sqrt(2*(E+mu/rE)); % Radius of the Earth, m/s, 2-8-12
H=sqrt(mu*a*(1-posEps(length(posEps))^2)); % Formula 2-8-2
phi=acosd(H/(rE*vLnch)); % Formula 2-8-6 

vLnchd=['Launch Velocity is: ',num2str(vLnch),' m/s'];
phid=['Corresponding Elevation Angle is: ',num2str(phi),' °'];
disp(vLnchd);
disp(phid);


