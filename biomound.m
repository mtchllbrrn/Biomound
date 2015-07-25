pause on;
clear all;
 
graph = 2;   % 1 = show movie, 0 = graph at end, 2 = don't graph at all
avi = 0; % Create .avi file
shrubs = 1;
fig = figure('position',[250 200 800 700]);
 
%Basic Model Parameters
dta = 1; % shrub population time step (yr)
dtb = 0.02; % geomorphology time step (yr)
tmax = 250; % run time (yr)
kmax = tmax/dtb;
 
%Physical Dimensions of Modeled Area
X = 10; % downslope length (m)
Y = 10; % cross-slope length (m)
Area = X*Y; % Total area (m^2)
 
% Nx = 101; % number of length elements
% Ny = 101; % number of width elements
Nx = X*20;
Ny = Y*20;
 
dx = X/(Nx - 1); % element spacing (m)
dy = Y/(Ny - 1); % element spacing
dZdx = .0035*X; % background slope
 
%Site specific climatic parameters
WPmu = 0.125; %Long-term averge annual winter precipitation (m/yr)
WPsigma = 0.025; %Standard deviation of long-term average annual winter precipitation (m/yr)
phi = 0.15; %Constant 0<phi<1 describing the soil's propensity to retain moisture from year to year
SPmu = 0.125; %Long-term average summer precipitation (m/yr)
SPsigma = 0.025; %Standard deviation of long-term average summer precipitation (m/yr)
ds = 180; %Length of growing season (180 days)
MAP = WPmu + SPmu; %Site's Mean Annual Precipitation (Parameter used to estimate rooting depth)
 
%Site-specific, a priori, estimate of rooting depth using climatic parameters and allometric relationships (From Schneck 2002)
RDMax = 10^(-0.3857 + 0.2412*log10(MAP*1000)); %Max rooting depth (m) where MAP is in mm
 
%Soil profile differentiation and parameters
UL = 0.20; %Depth of upper soil layer (m), pertinent to recruitment and juvenile mortality
LL = RDMax - UL; %Depth of lower layer soil compartment (m), pertinent to adult mortality
 
%Description of initial population and population specific parameters
shrubdensitymu = 0.5; % Average number of individuals per square meter (DEFAULT .50)
shrubdensitysigma = 0.00; %Standard deviation in the number of individuals per square meter (SHRUBDENSITY VARIABLES ARE ONLY USED IN THE INITIAL RANDOMIZATION OF SHRUB LOCATIONS)
minshrubspacing = 0.05; %This is the minimum distance that extant shrubs can be from one another (DEFAULT .20) (Used after initial randomization)
 
%Parameters for characterizing individual shrubs
 
%Assuming a gamma distribution, parameters used to assign radius' to individuals comprising the initial population
A = 6; %Parameters estimated from a GamFit on rabbitbursh radius field data (DEFAULT 6)
B = .05;
%Canopy dimensions and parameters for logistic growth model Furbish 2009
R0 = 0.05; %Initial Radius (m) assigned to shrubs with age = 1 yr
Rf = 0.5; %"Final" Radius (m) that can be achieved by a shrub
t90 = 10; %Age when a shrub reaches 90% of maximum radius (yr)
T = (-t90)/log(1-((0.90*Rf - R0)/(Rf - R0))); % Characteristic growth rate constant
RDCVRatio = RDMax/((4/3)*pi*(Rf)^3); %Rooting depth to canopy volume ratio (m^-2), where it is assumed that RD/CV is constant throughout shrub ontogeny
 
%Mortality Parameters
FC = 0.14; %
PWP = 0.06; %Permanent wilting point (Volumetric Water Content, m/m)
m = 2; %Parameter describing the species' sensitivity to water-stress
betaJ = .40; %Juvenile mortality rate constant, where 0 < beta < 1
betaA = .30; %Adult mortality rate constant, where 0 < beta < 1
 
%Recruitment Parameters
smu = 10000; %Average number of seeds produced per shrub
ssigma = 0; %Annual variability in seed production
alpha = 0.0055; %Recruitment rate constant, what is the range on this?
 
% Transport variables
K = 0.004; % transport coefficient (m/yr)
D = 0.0004; % diffusivity (m^2/yr)
 
% Finite difference constants
Kx = K/(2*dx^2);
Dx = D/(dx^2);
Ky = K/(2*dy^2);
Dy = D/(dy^2);
 
T90 = 3;   % time to 90% of canopy maturation (yr)
TT = -T90/(log(1-(0.9*Rf-R0)/(Rf-R0)));
Maxcover = 0.8; % maximum cover at center of shrub
 
%% Description of Initial Population
shrubdensity = shrubdensitymu;
N0 = round(Area * shrubdensity); % initial number of shrubs
 
%% Pre-allocating memory
% Mov = zeros(1,kmax,'struct');
maxrad = 0; %initializing max-radius querying script to zero (script is at bottom)
Mov = moviein(kmax,gcf);
% x = zeros(X*30);
% y = zeros(Y*30);
% zold = zeros(X*30,Y*30);
status = zeros(1,N0);
radius = zeros(1,N0);
canopyvolume = zeros(1,N0);
RD = zeros(1,N0);
age = zeros(1,N0);
xc = zeros(1,N0);
yc = zeros(1,N0);
time = zeros(1,kmax);
g = zeros(Nx,Ny);
% znew = zeros(X*30,Y*30);
Nplot = zeros(1,kmax);
AverageAge = zeros(1,kmax);
% finalagecount = zeros(1,kmax);
% CrestHeight = zeros(kmax,1);
SP = zeros(1,kmax);
WS = zeros(1,kmax);
WP = zeros(1,kmax);
excelcount = 0; % Counter for writing flux data to excel
 
% Sets initial coordinates and surface.  This assigns each node a real-world length value based on its distance from the origin.
for i=1:Nx
    for j=1:Ny
        x(i,j) = (i-1)*dx;
        y(i,j) = (j-1)*dy;
    end
end
 
heightadjust = -dZdx*(x(1,1)-X/2)^2;
 
for i=1:Nx
    for j=1:Ny
        %                 zold(i,j) = dZdx*x(i,j);
        zold(i,j) = -dZdx*(x(i,j)-X/2)^2-heightadjust;
    end
end
 
if shrubs == 1
    N = N0;
    % INITIAL SHRUB POSITIONING
    for n=1:N0;
        
        %Shrub attributes
        status(n) = 1; % 1 = live, 0 = dead
        
        %Generates random radius based on gamma distribution.  If the radius is
        %too large (That is, radius(n) > Rf), a new value is generated.
        radiusflag = 1;
        while radiusflag == 1
            radius(n) = gamrnd(A,B); % Assigns each shrub a radius from known(fitted) gamma distribution of the population's radius (m)
            if radius(n) >= Rf
                radiusflag = 1;
            else
                radiusflag = 0;
            end
        end
        
        canopyvolume(n) = (4/3)*pi*(radius(n))^3; %Assigns canopy volume based on radius, assuming the shrub is a sphere (m^3)
        RD(n) = canopyvolume(n)*RDCVRatio;  %Assigns the shrub a rooting depth from allometric scaling relationship (Schenck 2002)
        age(n) = round(log(1-((radius(n) - R0)/(Rf - R0)))*(-T));  %Assigns each shrub an age from the logistic growth function
        
        %Positioning of shrub
        xc(n) = random('unif',0,X); % x-coordinate of shrub center
        yc(n) = random('unif',0,Y); % y-coordinate of shrub center
        temp = n;
        if n > 1; %This tests positioning with respect to previously planted shrubs
            b = 1;
            attempt = 1;
            while b <= n - 1;
                positionflag = 1;
                d = sqrt((xc(n)-xc(b))^2 + (yc(n)-yc(b))^2)- radius(n) - radius(b);
                if d <= minshrubspacing; %Repositions shrub that are too close to extant shrubs
                    xc(n) = random('unif',0,X);
                    yc(n) = random('unif',0,Y);
                    positionflag = 0;
                end
                b = b + 1;
                if positionflag == 0;
                    b = 1;
                    attempt = attempt +1;
                    if attempt > 500;
                        error('cannot position shrub')
                    end
                end
            end
        end
        
        
    end
end
 
%%
 
 
N = N0; %This is the initial population count
v = 0; %Movie frame counter
moviecounter = 0;
fileprefix = 'Biomound';
g(:,:) = 1;
 
%TIME LOOP
for k=1:kmax
    tic
    time(k) = k*dtb;
 
    %% Geomorph
    
    for i=1:Nx   % nodes' x-values
        for j=1:Ny  % nodes' y-values
            g(i,j) = 1;
            for n = 1:N
                if abs(x(i,j)-xc(n)) <= 0.75 && abs(y(i,j)-yc(n)) <= 0.75  % if node lies somewhere near this shrub (note that we're not finding distance.  This is a faster way to find which shrubs should be considered)
                    dist = sqrt((x(i,j)-xc(1,n))^2 + (y(i,j)-yc(1,n))^2);
                    if dist <= radius(n)                                   % if current node is within canopy
                        C = 1 - dist^2/radius(n)^2;                        % parabolic function representing canopy cover (always a value between 0 and 1 because dist < radius.  Larger values indicate denser cover)
                        g(i,j) = g(i,j) - Maxcover*C;                           % sets activity based on corresponding canopy cover
                    end
                    if g(i,j) < 0.2
                        g(i,j) = 0.2;
                    end
                end
            end
        end
    end
    
    
    for i=2:Nx-1  % finite differences
        for j=2:Ny-1
            Tx1 = Kx*((g(i,j)+g(i+1,j))*(zold(i+1,j)-zold(i,j))...
                -(g(i,j)+g(i-1,j))*(zold(i,j)-zold(i-1,j)));
            Tx2 = Dx*(g(i+1,j)-2*g(i,j)+g(i-1,j));
            Ty1 = Ky*((g(i,j)+g(i,j+1))*(zold(i,j+1)-zold(i,j))...
                -(g(i,j)+g(i,j-1))*(zold(i,j)-zold(i,j-1)));
            Ty2 = Dy*(g(i,j+1)-2*g(i,j)+g(i,j-1));
            znew(i,j) = zold(i,j) + dtb*(Tx1+Tx2+Ty1+Ty2);
        end
    end
    
    % Boundary Conditions
    for i = 2:Nx-1
        Txb1 = Kx*((g(i,j)+g(i+1,j))*(zold(i+1,j)-zold(i,j))...
            -(g(i,j)+g(i-1,j))*(zold(i,j)-zold(i-1,j)));
        Txb2 = Dx*(g(i+1,j)-2*g(i,j)+g(i-1,j));
        Tyb1 = Ky*((g(i,j)+g(i,j+1))*(zold(i,j+1)-zold(i,j))...
            -(g(i,j)+g(i,j+1))*(zold(i,j)-zold(i,j+1)));
        Tyb2 = Dy*(2*g(i,j+1)-2*g(i,j));
        znew(i,j) = zold(i,j) + dtb*(Txb1+Txb2+Tyb1+Tyb2);
    end
    
    for i = 2:Nx-1
        j = Ny;
        Txb3 = Txb1;
        Txb4 = Txb2;
        Tyb3 = Ky*((g(i,j)+g(i,j-1))*(zold(i,j-1)-zold(i,j))...
            -(g(i,j)+g(i,j-1))*(zold(i,j)-zold(i,j-1)));
        Tyb4 = Dy*(2*g(i,j-1)-2*g(i,j));
        znew(i,j) = zold(i,j) + dtb*(Txb3+Txb4+Tyb3+Tyb4);
    end
    
    for i=2:Nx-1   % replacement - Interior nodes
        for j=2:Ny-1
            zold(i,j) = znew(i,j);
        end
    end
    
    for i = 2:Nx-1  % boundary nodes along y = 1
        for j = 1
            zold(i,j) = znew(i,j);
        end
    end
    
    for i = 2:Nx-1  % Boundary nodes along y = Ny (ymax)
        for j = Ny
            zold(i,j) = znew(i,j);
        end
    end
    
    % Canopy growth
    if shrubs == 1
        for i = 1:N
            age(i) = age(i) + dtb; %The age of this aged shrub is simply it's old age plus the time step...
            dRdtmu = (exp(1)^(-age(i)/T))/T*(Rf-R0); %Average age specific annual growth rate
            dRdtsigma = .15 * dRdtmu; %Average age specific standard deviation in growth rate
            dRdt = normrnd(dRdtmu,dRdtsigma); %Actualized growth drawn from normal distribution
            radius(i) = radius(i) + dRdt*dtb; %Grow canopy
            canopyvolume(i) = (4/3)*pi*(radius(i))^3;  %Calculate new canopy volume
            RD(i) = canopyvolume(i)*RDCVRatio;
        end
    end
    
    %Stat tracking
    if shrubs == 1
        Nplot(k) = N;
        xx = zeros(N,1);
        yy = zeros(N,1);
        rr = zeros(N,1);
        gg = zeros(N,1);
        AgeCount = 0;
        for n = 1:N;
            xx(n) = xc(n);
            yy(n) = yc(n);
            rr(n) = radius(n);
            gg(n) = age(n);
            AgeCount = AgeCount + gg(n);
        end
        AverageAge(k) = AgeCount/N;
        
        agemax = 75;
        for w = 1:agemax;
            AGE = w;
            agecount = 0;
            for n = 1:N;
                if gg(n) >= AGE && gg(n) < AGE+1;
                    agecount = agecount + 1;
                end
            end
            finalagecount(w) = agecount;
        end
 
    end
    
    
    %BEGIN STACEY'S SECTION:  This portion of the code only runs every 50th
    %step (with dtb = 0.02).
    if shrubs == 1
        if rem(k,1/dtb) == 0
            %Recruitment/Juvenile Mortality Variables
            
            
            %Different Possible Summer Precipitation Scenarios:  Turn on/of as needed
            SP(k)= SPmu; % Summer precipitation
            %SP(k) = normrnd(SPmu,SPsigma); %Summer precipitation received which will be deliver to the upper soil layer (m)
            %SP(k) = .05*sin((pi()/100)*k)+.125;
            
            ULthetaave = SP(k)/(UL*ds) + PWP; %Average daily water content of upper soil layer assuming water is distributed uniformly with depth of the compartment (m/m, unitless)
            
            if SP(k) < 0;
                SP(k) = 0;
            end
            
            %Adult Mortality Variables
            
            %Different Possible Winter Precipitation Scenarios: Turn on/off as needed
            WP(k) = WPmu;
            %WP(k) = normrnd(WPmu,WPsigma);
            %WP(k) = .05*sin((pi()/100)*k)+.125;
            %WP(k) = phi*WP(k-1) + normrnd(WPmu,WPsigma); % Winter precipitation recieved in year k that will be contributed to the lower soil layer = some fraction "phi" of k-1's precipitation + average WP + some random shock (First-order autoregressive process) (m)
            if WP(k) < 0;
                WP(k) = 0;
            end
            LLthetaave = (WP(k)/LL); %Average water content of lower soil layer, if water is distributed uniformly with depth of the compartment (m/m, unitless)
            
            %Boundary Conditions:
            if LLthetaave < 0; %Precipitation cannot be a negative value
                LLthetaave = 0;
            end
            if LLthetaave > FC-PWP; %Precipitation input cannot contribute more than the soil's water holding capacity
                LLthetaave = FC - PWP;
            end
            
            %Gamma describes how the lower layer's average moisture-content is distributed with depth (m^-1)
            gamma = (((ULthetaave - PWP)/LLthetaave)-1)*(-1/LL);
            
            %% Juvenile and Adult Mortality
            
            PJM = betaJ*(PWP/(ULthetaave))^m; %Probability of Juvenile Mortality in the year k is proportional to the adequacy of the upper soil layer's water content 0 < PJM < 1
            
            for n=1:N;
                %For "juvenile" shrubs, dependent on upper level moisture/summer precipitation
                if RD(n) <= UL; %If n's rooting depth is within the confines of the upper soil layer
                    Death = rand();
                    if Death <= PJM;
                        status(n) = 0;
                    end
                    if Death > PJM;
                        status(n) = 1;
                    end
                end
                
                %For "adult" shrubs, dependent on lower level moisture/winter precipitation
                if RD(n) > UL; %If n has a root system long enough to access the lower soil layer
                    theta(n) = LLthetaave*(1 + gamma*((RD(n)-UL) - LL)) + PWP; %#ok<SAGROW> %This is the water content at n's rooting depth, adjusting for the fact that water is not distributed uniformly with depth
                    PMA(n) = betaA*(PWP/theta(n))^m; %#ok<SAGROW> %Probability of mortality for adult
                    Death = rand();
                    if Death <= PMA(n);
                        status(n) = 0;
                    end
                    if Death > PMA(n);
                        status(n) = 1;
                    end
                end
            end
            
            
            %% Population count: Renumbering and aging all survivors, canopy growth
            % This loop eliminates all dead shrubs from the position, age, radius,
            % volume, and rooting depth arrays.  It then 'tidies' the arrays such that
            % any extraneous values left over from dead shrubs are pruned.
            
            count = 0;
            for n=1:N; %For all shrubs...
                temp = status(n);
                if temp == 1 %if n was a survivor
                    count = count + 1; %Count it
                    xc(count) = xc(n); %The x-position of the surviving shrub is simply it's original position
                    yc(count) = yc(n); %The y-position of the surviving shrub is simple it's original position
                    age(count) = age(n);
                    radius(count) = radius(n);
                    canopyvolume(count) = canopyvolume(n);
                    RD(count) = RD(n);
                    
                    %                 age(count) = age(n) + dta; %The age of this aged shrub is simply it's old age plus the time step...
                    %                 dRdtmu = (exp(1)^(-age(count)/T))/T*(Rf-R0); %Average age specific annual growth rate
                    %                 dRdtsigma = .15 * dRdtmu; %Average age specific standard deviation in growth rate
                    %                 dRdt = normrnd(dRdtmu,dRdtsigma); %Actualized growth drawn from normal distribution
                    %                 radius(count) = radius(n) + dRdt*dta; %Grow canopy
                    %                 canopyvolume(count) = (4/3)*pi*(radius(count))^3;  %Calculate new canopy volume
                    %                 RD(count) = canopyvolume(count)*RDCVRatio;
                    
                end
            end
            Nsurvivors = count; %Number of survivors
            
            %Tidying up the shrub description vectors
            status = [];
            for n = 1:count;  %Tidies status vector
                status = [status 1]; %#ok<AGROW>
            end
            xc(Nsurvivors+1:end) = [];  %Tidies shrub position vectors
            yc(Nsurvivors+1:end) = [];
            age(Nsurvivors+1:end) = [];
            radius(Nsurvivors+1:end) = [];
            canopyvolume(Nsurvivors+1:end) = [];
            RD(Nsurvivors+1:end) = [];
            
            
            %% Recruitment Variables
            
            %Calculate the number of seed producing shrubs ("Adults")
            Nseed = 0;
            for n = 1:Nsurvivors;
                if RD(n) > UL; %Does the shrub have access to LL moisture?
                    Nseed = Nseed + 1;
                end
            end
            
            %Stochastic variable:  Number of seeds produced per shrub
            s = round(normrnd(smu,ssigma));
            
            %Calculates the proportion of the modeled area occuppied by the extant population; "effective area" (EA) precludes the addition of new individuals (m^2)
            EAPop = 0;
            for n = 1:Nsurvivors;
                EA(n) = 1.5 * pi()*(radius(n))^2; %#ok<SAGROW> %Effective area of shrub n (m^2)
                EAPop = EAPop + EA(n);
            end
            
            SPACE = (Area - EAPop)/Area; %Area available for new recruits
            WATER = (ULthetaave-PWP)/ULthetaave; %Water available to seedlings; This should have it's own formulation?
            
            %% Final Recruitment Statement
 
                Nrecruited = round(alpha * WATER * SPACE * Nseed* s); %Number of shrubs that we are adding is simply some proportion of those already there.
            overflow = 0;
            
            for n = Nsurvivors+1:Nsurvivors+Nrecruited;
                
                %Assigns all new recruits dimensions
                status(n) = 1; %#ok<SAGROW>
                age(n) = 0;
                radius(n) = R0;
                canopyvolume(n) = (4/3)*pi*(radius(n))^3;
                RD(n) = canopyvolume(n)*RDCVRatio;
                
                %Positioning of new recruit
                attempts = 0;
                flag = 0;
                while flag == 0
                    success = 0;
                    xc(n) = random('unif',0,X); % x-coordinate of shrub center
                    yc(n) = random('unif',0,Y); % y-coordinate of shrub center
                    for b = 1:n-1
                        d = sqrt((xc(n)-xc(b))^2 + (yc(n)-yc(b))^2)- radius(n) - radius(b);
                        if d < minshrubspacing; %Repositions shrub that are too close to extant shrubs
                            attempts = attempts + 1;
                            break
                        else
                            success = success + 1;
                        end
                    end
                    if success == n-1
                        flag = 1;
                    end
                    if attempts > 1000
                        overflow = 1;
                        overflowloc = n;
                        flag = 1;
                    end
                end
                
                if overflow == 1
                    xc = xc(1:overflowloc-1);
                    yc = yc(1:overflowloc-1);
                    status = status(1:overflowloc-1);
                    age = age(1:overflowloc-1);
                    radius = radius(1:overflowloc-1);
                    canopyvolume = canopyvolume(1:overflowloc-1);
                    RD = RD(1:overflowloc-1);
                    Nrecruited = N-overflowloc-1;
                    break
                end
            end
            
            N = Nsurvivors + Nrecruited;
            
        end
    end
    
    % QUERYING FLUX
    % Bottom of slope, near side
    for j = 2:Ny-1
        i = 1;
        q(i,j) = K*g(i,j)*(zold(i,j)-zold(i+1,j))/dx - D*(g(i,j)-g(i+1,j))/dx;
    end
    
    % Middle of hillslope
    for j = 2:Ny-1
        i = ceil(Nx/4);
        q(i,j) = K*g(i,j)*(zold(i,j)-zold(i+1,j))/dx - D*(g(i,j)-g(i+1,j))/dx;
    end
    
    % Crest of hillslope
    for j = 2:Ny-1
        i = ceil(Nx/2);
        q(i,j) = K*g(i,j)*(zold(i,j)-zold(i+1,j))/dx - D*(g(i,j)-g(i+1,j))/dx;
    end
    
    
   % Assigning queried fluxes to arrays
    if rem(k,25) == 0
        excelcount = excelcount + 1;
        avgfluxbottom(excelcount,1) = mean(q(1,2:Ny-1));
        avgfluxmiddle(excelcount,1) = mean(q(ceil(Nx/4),2:Ny-1));
        avgfluxtop(excelcount,1) = mean(q(ceil(Nx/2),2:Ny-1));
        CrestHeight(excelcount,1) = mean(znew(Nx/2-1,1:Ny-1));
    end
    
    %Querying largest shrub diameter (checking for realistic sizes)
    for i = 1:size(radius,2)
        if radius(1,i) > maxrad
            maxrad = radius(1,i);
        end
    end
    
    
    %%
    if graph == 1
        if shrubs == 1
            subplot(4,2,1:4)
            surf(x,y,zold,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
            axis equal
            set(gca,'DataAspectRatio',[1 1 1])
            view(-40,20)
            camlight right
            
            subplot(4,2,5);
            scatter(xx,yy,round(rr*100),'markerfacecolor','g');
            axis equal;
            axis([0,X,0,Y]);
            
            subplot(4,2,6);
            plot(time, Nplot);
            axis([1,tmax,0,6*N0]);
            xlabel('Year')
            ylabel('shrub count')
            
            subplot(4,2,7);
            bar(1:agemax,finalagecount);
            xlabel('age')
            ylabel('shrub count')
            axis([0 50 0 20]);
            
            subplot(4,2,8);
            plot(time,CrestHeight);
            axis([1,tmax,-.2,.2]);
            xlabel('Year')
            ylabel('Avg Height')
            
            %Activity display
            %             subplot(4,2,[5:8])
            %             surf(x,y,-g,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
            %             axis equal
            %             set(gca,'DataAspectRatio',[1 1 1])
            %             view(-40,20)
            %             camlight right
            
        else
            subplot(4,2,1:4)
            surf(x,y,zold,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
            axis equal
            set(gca,'DataAspectRatio',[1 1 1])
            view(-40,20)
            camlight right
            
            subplot(4,2,8);
            plot(time,CrestHeight);
            axis([1,tmax,-.2,.2]);
            xlabel('Year')
            ylabel('Avg Height')
            
        end
        if avi == 1
            if rem(k,10) == 0
                v = v + 1;
                Mov(v) = getframe(fig);
            end
        end
    end
    
    if rem(k,25) == 0
        timeleft = dispTimeLeft(1,1,kmax,k);
    end
    
end  %end time loop
 
if graph == 0
    if shrubs == 1
        clf('reset')
        
        subplot(4,2,1:4)
        surf(x,y,zold,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
        axis equal
        set(gca,'DataAspectRatio',[1 1 1])
        view(-40,20)
        camlight right
        
        subplot(4,2,5);
        scatter(xx,yy,round(rr*100),'markerfacecolor','g');
        axis equal;
        axis([0,X,0,Y]);
        
        subplot(4,2,6);
        plot(time, Nplot);
        axis([1,tmax,0,3*N0]);
        xlabel('Year')
        ylabel('shrub count')
        
        subplot(4,2,7);
        bar(1:agemax,finalagecount);
        xlabel('age')
        ylabel('shrub count')
        axis([0 50 0 20]);
        
        subplot(4,2,8);
        plot(time, CrestHeight);
        axis([1,200,-.2,.2]);
        xlabel('Year')
        ylabel('Avg Height')
        
    else
        subplot(4,2,1:4)
        surf(x,y,zold,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
        axis equal
        set(gca,'DataAspectRatio',[1 1 1])
        view(-40,20)
        camlight right
        
        subplot(4,2,8);
        plot(time, CrestHeight);
        axis([1,200,-.2,.2]);
        xlabel('Year')
        ylabel('Avg Height')
    end
end
 
if avi == 1
    movie2avi(Mov,'BioMound','compression','CinePak','quality', 10, 'fps', 30)
end
 
% Writing data to Excel sheets
xlswrite('HillslopeStatistics',CrestHeight,1)
xlswrite('HillslopeStatistics',avgfluxbottom,2)
xlswrite('HillslopeStatistics',avgfluxmiddle,3)
xlswrite('HillslopeStatistics',avgfluxtop,4)
 
