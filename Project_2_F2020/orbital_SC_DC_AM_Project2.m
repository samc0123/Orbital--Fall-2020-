%% Introduction
% Samuel Chernov(sc1823)
% Daniella Chung(djc375)
% Andrew Manfredi(ajm418) 

% Orbital Mechanics: Project 2 
% Professor Mookerjee
% Fall 2020 

% Problem #1

%% Part 0: Initializing constants 
rE=6.378e6; % Earth's Radius, m
mu=3.986e14; % m^3/s^2-> value for Earth 
h1=4663e3; % m 
h2=6352e3; % m 
theta=78.5; % deg 
t_break=44.5-0.033; % min 

%% Part A %%
r1=rE+h1; %m 
r2=rE+h2; %m 
a=(r1+r2)/2; %m 

% Initialize space for values 
TOF=zeros(1,25); 
aPltVals=zeros(1,25);

% Main Calculation 
for i=1:25
    aTemp=0; % a value for this iteration 
    
    if i==1
       aTemp=a; % calculate the maximum value 
    elseif(i~=1)
        factor=(i-1)*0.043;
        aTemp=factor*a+a;
    end
    
    % Save 'aTemp' value for plotting 
    aPltVals(i)=aTemp;
    
    % Constant 'P'
    P=sqrt(aTemp^3/mu);

    % Constant 'd' 
    d=sqrt(r1^2+r2^2-2*r1*r2*cosd(theta));

    % Finding alpha 
    tempSqrtA=sqrt((r1+r2+d)/(aTemp));
    tempSqrtA=tempSqrtA*0.5; % Finalizes right side of equal sign 

    alpha=2*asin(tempSqrtA); % alpha value, degrees

    % Finding beta
    tempSqrtB=sqrt((r1+r2-d)/(aTemp));
    tempSqrtB=tempSqrtB*0.5; % Finalizes right side of equal sign 

    beta=2*asin(tempSqrtB); % degrees value, degrees

    % Obtain TOF 
    tempTOF=P*((alpha-sin(alpha))-(beta-sin(beta))); % s 
    TOF(i)=tempTOF/60; % min
end

% Plot the Figure 
figure
hold on 
plot(aPltVals,TOF,'LineWidth',2);
plot(aPltVals,TOF,'g*','LineWidth',2);
xlabel('Semi-major Axis Value (m)');
ylabel('TOF (min)');
title('TOF vs. Semi-major axis length');


%% Part B

% Initialize space for the variables 
TOF_break=[];
aPltVals_b=[];

% Initilaize the first TOF & 'a' value from Part (a) 
TOF_break(1)=TOF(1); % initial TOF is the same 
aPltVals_b(1)=a; % intial 'a' is the same 

%{
* While Loop
* Decrease magnitude of 'a' by 0.0043
* Break when TOF falls bellow 't_break' 
* First iteration already calcualated, start loop
from iteration 2
%}
i=1; % Counter
while TOF_break(i)>t_break
   
   
   % Get 'a' for the current calculation
   aTemp_b=0.0043*a+aPltVals_b(i); 
   aPltVals_b(i+1)=aTemp_b;
   
    % Constant 'P'
    P_b=sqrt(aTemp_b^3/mu);

    % Constant 'd' 
    d_b=sqrt(r1^2+r2^2-2*r1*r2*cosd(theta));

    % Finding alpha 
    tempSqrtA_b=sqrt((r1+r2+d_b)/(aTemp_b));
    tempSqrtA_b=tempSqrtA_b*0.5; % Finalizes right side of equal sign 

    alpha_b=2*asin(tempSqrtA_b); % alpha value, degrees

    % Finding beta
    tempSqrtB_b=sqrt((r1+r2-d_b)/(aTemp_b));
    tempSqrtB_b=tempSqrtB_b*0.5; % Finalizes right side of equal sign 

    beta_b=2*asin(tempSqrtB_b); % degrees value, degrees

    % Obtain TOF 
    tempTOF_b=P_b*((alpha_b-sin(alpha_b))-(beta_b-sin(beta_b))); % s 
    TOF_break(i+1)=tempTOF_b/60; % hr
    
    % Check time condition 
    if TOF_break(i+1)<t_break
        break
    end
    
    % Increment the counter
    i=i+1;
   
end

% Plot on the same curve
plot(aPltVals_b(length(aPltVals_b)-1),TOF_break(length(TOF_break)-1),'rx','LineWidth',3)
legend('TOF Curve','Individual TOF Pts.','TOF Break Pt.')
hold off

%% Part C

% Obtain the value of the breakpoint for 'A' 
a_c=aPltVals_b(length(aPltVals_b)-1); 

% Use this value to get the new alpha_c and beta_c (Form 4-2-2 a&b)
alpha_c=1-((r1+r2+d)/(2*a_c)); 
beta_c=1-((r1+r2-d)/(2*a_c)); 

alpha_c=acos(alpha_c);
beta_c=acos(beta_c);

% Utilize Formula 4-5-7 to find psi 
psi=alpha_c-beta_c;

% Calculate ua from formula 4-5-4
fracTemp_ua=(a_c-r2)/(a_c-r1);
multTemp_ua=cos(psi)-fracTemp_ua;

ua=atan((1/sin(psi))*multTemp_ua); 

% Calculate ub using 4-5-7
ub=psi+ua;

% Get the eccentricity from 4-5-15
% Will be same for any pt., so use pt. A to calculate it
% ua-> r1 &&&& ub->r2, therefore: 
eps_c=(a_c-r1)/(a_c*cos(ua)); 

% Get Energy from Ch.2 formulas 
E_c=-mu/(2*a_c); % m^2/s^2

% Get the specific angular momentum 'H' from Ch. 2 formulas 
H_c=sqrt(a_c*mu*(1-eps_c^2)); %m^2/s

% For Va, must calculate 'ra'
% Use this & 'H' to get Va 
ra_c=a_c*(1+eps_c); % m 
Va_c=H_c/ra_c; % m/s 

%% Part D
%{
* Look @ the bottom of pg 75 in text 
* Calculate the 'p' value, and then 'v' from there
* That is used to calculate the position angle 
* Velocity is simply a chap. 2 formula
%}

% Calculate the velocities based on Chapter 2 formulas 
V1_d=sqrt(2*(E_c+mu/r1)); % m/s 
V2_d=sqrt(2*(E_c+mu/r2)); % m/s 

% Calculate the semilatus rectum 'P'
P_d=a_c*(1-eps_c^2); %m 

% Calculate the position 'v' 
v1_d=acos((1/eps_c)*(P_d/r1-1)); %rad
v2_d=acos((1/eps_c)*(P_d/r2-1)); %rad

% Calculate the elevation angles, phi from Ch.2 
phi1_d=atan((eps_c*sin(v1_d))/(1+eps_c*cos(v1_d))); % rad
phi2_d=atan((eps_c*sin(v2_d))/(1+eps_c*cos(v2_d))); % rad 

%% Part E 

% Find Ma_e & Mb_e from middle of pg 75
Ma_e=ua-eps_c*sin(ua);
Mb_e=ub-eps_c*sin(ub);

% TOF eqn. as sqrt(a^3/mu)*(Mb-Ma)
TOF_e=sqrt(a_c^3/mu)*(Mb_e-Ma_e);
TOF_e=TOF_e/60; % put into minutes

%% Part F

% Angle between the two observations
% Look at Fig. 4-5-1 on pg 72 in the text
theta_f=(v2_d-v1_d);
theta_f=theta_f*180/pi; % Convert to degrees, like given in assignment

%% Part G 
% Orbital Drawing

% Step 0- initialize values for orbit plotting
gammaPlt=linspace(0,2*pi,50);
figure 
hold on

% Step 1- draw the earth
center=[0,0];
ePlt=viscircles(center,rE,'LineWidth',3,'Color','b');

% Step 2- get values necessary to create orbit plot 

% 'r' value from Eqn. 2-5-5
eps_g=eps_c; % eccentricity value
a_g=a_c; % value of 'a' used for subsequent calcualtions 

rxPlt=zeros(1,length(gammaPlt));
ryPlt=zeros(1,length(gammaPlt));

for i=1:length(rxPlt)
    rTemp=(a_g*abs(1-eps_g^2))/(1+eps_g*cos(gammaPlt(i)));
    
    rxPlt(i)=cos(gammaPlt(i))*rTemp;
    ryPlt(i)=sin(gammaPlt(i))*rTemp;
end
orbPlt=plot(rxPlt,ryPlt,'LineStyle','--','Color','g','LineWidth',3);

% Step 3- Mark the observations 

%{
* Observations
    * Occur @ v1 and v2 (angles)
    * Occur @ r1 and r2 (positions) 
* Plotting the points of those observations 
%}

% Observation 1 
obs1x_g=cos(v1_d)*r1;
obs1y_g=sin(v1_d)*r1;

% Observation 2 
obs2x_g=cos(v2_d)*r2;
obs2y_g=sin(v2_d)*r2;

% Plot the Observations 
obs1x_plt_g=[0 obs1x_g];
obs1y_plt_g=[0 obs1y_g]; 
obs2x_plt_g=[0 obs2x_g]; 
obs2y_plt_g=[0 obs2y_g];

obs1Plt=plot(obs1x_plt_g,obs1y_plt_g,':','Color','r','LineWidth',3);
obs2Plt=plot(obs2x_plt_g,obs2y_plt_g,':','Color','m','LineWidth',3);
plot(obs1x_g,obs1y_g,'x','Color','r','LineWidth',4);
plot(obs2x_g,obs2y_g,'o','Color','m','LineWidth',4);

% Step 4- label the figure 

% Axes
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% Labels
title('Earth and Orbit Trajectory Plot');
xlabel('Location (m)');
ylabel('Location (m)');
legend([ePlt orbPlt obs1Plt obs2Plt],{'Earth','Orbit Trajectory',...
    'Observation 1','Observation 2'},...
    'Location','SouthEast');
hold off







