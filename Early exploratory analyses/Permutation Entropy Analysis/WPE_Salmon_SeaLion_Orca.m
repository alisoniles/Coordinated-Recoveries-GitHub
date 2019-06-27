%Objective: Calculate the weighted permutation entropy of salmon, orca, sea
%lion time series to gauge information content. 

%By: Alison Iles

%Date initiated: November 15, 2018


%% Load Data

clearvars
close all
cd('/Users/alisoniles/Google Drive/Coordinated Recoveries//Data/analysis ready/FPC annual dam counts/');
addpath(genpath('/Users/alisoniles/Google Drive/Coordinated Recoveries/Analysis/Coordinated-Recoveries-GitHub/Permutation Entropy Analysis'));
addpath(genpath('/Users/alisoniles/Google Drive/Coordinated Recoveries/Data/analysis ready/'));

headerlines=1;

% Fish Passage Center data in files by dam, then by species in columns

d=dir('*.csv');  % all 14 .csv files (one for each Dam) in working directory
N=length(d);     % how many did we find?
C=cell([N,12]);  % cell array to save raw daily count data in
Dam=strings(N,1);

for k=1:N
  [fid]=fopen(d(k).name);  % open for read permission
    
  C(k,:)=textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f', ...
            'delimiter',',',...
            'headerlines',1);
      
  fid=fclose(fid);
  
  Dam_temp=d(k).name;
  Dam(k)=Dam_temp(1:3); %Create a list of dam code names for later plotting
end


% remove WAN and WFA dams, not enough data to run analyses
C(14,:)=[]; %WFA : only 16 years
C(12,:)=[]; %WAN : only 12 years
Dam(14)=[];
Dam(12)=[];

% sort dams in order from mouth to headwaters. Columbia River then Snake River
Dam_No=[1;9;3;12;11;10;4;5;6;7;2;8];
[Dam_sorted, Dam_order] = sort(Dam_No);
C = C(Dam_order,:);
Dam = Dam(Dam_order,:);

        % Orca Data
        fid = fopen('analysis ready/SRKW.csv');
            [data] = textscan(fid,'%f%f','delimiter',',','headerlines',headerlines); fclose(fid);
SRKWy=cell2mat(data(1)); % year
SRKW=cell2mat(data(2)); % number of individuals in the southern resident orca population

        % Sea Lion Data
        fid = fopen('analysis ready/sealion.csv');
            [data] = textscan(fid,'%f%f%f','delimiter',',','headerlines',headerlines); fclose(fid);
CSLy=cell2mat(data(1)); % year
CSLp=cell2mat(data(2)); % number of pups counted or inferred in the Santa Catalina Islands
CSLm=cell2mat(data(3)); % number of males from populaiton model based on pup count

        % PDO data
        fid = fopen('analysis ready/PDO.csv');
            [data] = textscan(fid,'%f%f','delimiter',',','headerlines',headerlines); %fclose(fid);
PDOy=cell2mat(data(1)); % year
PDO=cell2mat(data(2)); % average PDO for year



%% WPE Analysis

WPE_results=zeros(size(C)); %Matrix to store WPE results
No_perm_results=zeros(size(C)); %Matrix to store the number of permutations behind each WPE value

wl=3; %word length of permutation, i.e. the length of the permutation order of successive time points
tau=1; %time delay between the 'letters' of each 'word'

for i=1:size(C,1)
    for j=2:size(C,2)
        if sum(cell2mat(C(i,j))==0)/length(cell2mat(C(i,j))==0)<0.1 %only proceed if there are less than 10% of zeros in the data
                [ WPE_results(i,j), No_perm_results(i,j) ] = WPE_gap_func(cell2mat(C(i,j)), cell2mat(C(i,1)), wl, tau); %WPE of Each salmon species for each Ddam
        end
    end
end

WPE_results(:,1)=[];
No_perm_results(:,1)=[];

%WPE_gap_func(CSL_p, CSL_y, wl, tau); %WPE of sea lion pup counts



%% Visualize WPE Results
f=figure(1);
hold all
set(f, 'Position', [1000, 1000, 600, 600]);
set(f, 'DefaultAxesFontSize', 18);
set(f, 'DefaultLineLineWidth',1.5);
Fish=["Chinook", "Coho","Steelhead","Unclipped Steelhead","Sockeye","Pink","Chum","Lamprey","Shad","Chinook SS","Chinook F"];  

plot(Dam_sorted,WPE_results(:,[1;10;11;2;3;5]),'.','MarkerSize',30);
xlabel('Dam No.')
ylabel('WPE')         
set(gca, ...
            'YLim', [0.1,1], ...
            'ytick',0.1:0.1:1, ...
            'xtick',1:1:length(Dam), ...
            'xticklabel',Dam, ...
            'Linewidth', 1, ...
            'Box','on'     ) 
       
%Legend
legend(Fish([1;10;11;2;3;5]), 'Location','southwest','box','on' );  



%% WPE for Killer whales, sea lions and PDO

[WPE_CSLp, No_perm_CSLp]=WPE_gap_func(CSLp, CSLy, wl, tau); %WPE of sea lion pup counts
[WPE_CSLm, No_perm_CSLm]=WPE_gap_func(CSLm, CSLy, wl, tau); %WPE of sea lion male estimates
[WPE_SRKW, No_perm_SRKW]=WPE_gap_func(SRKW, SRKWy, wl, tau); %WPE of killer whales
[WPE_PDO, No_perm_PDO]=WPE_gap_func(PDO, PDOy, wl, tau); %WPE of PDO

f=figure(2);
hold all
set(f, 'Position', [1000, 1000, 600, 600]);
set(f, 'DefaultAxesFontSize', 18);
set(f, 'DefaultLineLineWidth',1.5);

plot([1;2;3;4],[WPE_CSLp; WPE_CSLm; WPE_SRKW; WPE_PDO],'.','MarkerSize',30);
ylabel('WPE')         
set(gca, ...
            'YLim', [0.1,1], ...
            'ytick',0.1:0.1:1, ...
            'XLim', [0.8,4.2], ...
            'xtick',1:1:4, ...
            'xticklabel',["CSL pups", "CSL males","SR orcas","PDO"], ...
            'Linewidth', 1, ...
            'Box','on'     ) 
 
