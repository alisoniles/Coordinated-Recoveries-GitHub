%Objective: Calculate the weighted permutation entropy of salmon, orca, sea
%lion time series to gauge information content. 

%By: Alison Iles

%Date initiated: November 15, 2018


%% Load Data

clearvars
close all
cd('/Users/alisoniles/Google Drive/Coordinated Recoveries/Analysis/Coordinated-Recoveries-GitHub-master/Permutation Entropy Analysis');
addpath(genpath('/Users/alisoniles/Google Drive/Coordinated Recoveries/Data'));
addpath(genpath('/Users/alisoniles/Google Drive/Coordinated Recoveries/Analysis/EDM Analysis'));

headerlines=1;

%% Orca Data
fid = fopen('analysis ready/SRKW.csv');
    [data] = textscan(fid,'%f%f','delimiter',',','headerlines',headerlines); fclose(fid);
SRO_y=cell2mat(data(1)); % year
SRO_n=cell2mat(data(2)); % number of individuals in the southern resident orca population

%% Sea Lion Data
fid = fopen('raw/CA_SeaLionPupCounts_Laake2018.csv');
    [data] = textscan(fid,'%f%f%f%f%f%f','delimiter',',','headerlines',headerlines); fclose(fid);
CSL_y=cell2mat(data(1)); % year
CSL_n=cell2mat(data(2)); % number of pups counted or inferred in the Santa Catalina Islands
CSL_f=cell2mat(data(3)); % number of females from population model based on pup count
CSL_m=cell2mat(data(4)); % number of males from populaiton model based on pup count
CSL_t=cell2mat(data(5)); % total number of adults from population model based on pup count

%% Chinook Salmon data from COE passage through both ladders at Bonneville Dam
fid = fopen('raw/Fish Passage Center COE/Fish_Passage_Counts_BON_dam_1938_2018.csv');
    [data] = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',','headerlines',headerlines); %fclose(fid);
m=cell2mat(data(2)); % month
d=cell2mat(data(3)); % day
y=cell2mat(data(4)); % year
mtime=datenum(y,m,d); % Convert sampling dates to a serial date number - the number of days from January 0, 0000.

CHA=cell2mat(data(5)); % chinook adult
CHJ=cell2mat(data(6)); % chinook jack

i=(m>=3 & m<=11); %index the Spring-Summer-Fall months to keep; March-November, because of zeros during the winter
y=y(i);
mtime=mtime(i);
CHA=CHA(i);
CHJ=CHJ(i);

year=unique(y);
CHAy=grpstats(CHA,y,'sum');
CHJy=grpstats(CHJ,y,'sum');

% CohoA=cell2mat(data(7)); CohoA=CohoA(i); % coho adult
% CohoJ=cell2mat(data(8)); CohoJ=CohoJ(i); % coho jack
% Stealhead=cell2mat(data(9)); Stealhead=Stealhead(i); % stealhead
% Stealhead_unclip=cell2mat(data(10)); Stealhead_unclip=Stealhead_unclip(i); % steal head unclipped
% Sockeye=cell2mat(data(11)); Sockeye=Sockeye(i); % sockeye
% Pink=cell2mat(data(12));  Pink=Pink(i); % pink
% Chum=cell2mat(data(13));  Chum=Chum(i); % chum
% Lamprey=cell2mat(data(14));  Lamprey=Lamprey(i); % lamprey
% Shad=cell2mat(data(15));  Shad=Shad(i); %Shad

figure(1)
doy = day(datetime(y,m,d),'dayofyear');
plot(doy,CHA+CHJ, '.k'); hold on
plot(


%% Chinook Salmon data from COE passage at Ice Harbor Dam
fid = fopen('raw/Fish Passage Center COE/Fish_Passage_Counts_IHR_dam_1962_2018.csv');
    [data] = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',','headerlines',headerlines); %fclose(fid);
m=cell2mat(data(2)); % month
d=cell2mat(data(3)); % day
y=cell2mat(data(4)); % year
IHR_mtime=datenum(y,m,d); % Convert sampling dates to a serial date number - the number of days from January 0, 0000.
IHR_CHA=cell2mat(data(5)); % chinook adult
IHR_CHJ=cell2mat(data(6)); % chinook jack

figure(1)
doy = day(datetime(y,m,d),'dayofyear');
plot(doy,IHR_CHA, '.k')

IHR=horzcat(IHR_mtime,m,d,y,IHR_CHA,IHR_CHJ);

%Split into spring/summer chinook and fall chinook and sum by year
%For Ice Harbor Dam, Spring/summer chinook are before 8/11, fall chinook
%are 8/12 and later
Si=logical((m==8 & d<=11)+(m<=7)); %index the Spring-Summer run;
IHR_S=IHR(Si,:);
year=unique(IHR_S(:,4));
IHR_Stotal=grpstats(IHR_S(:,5:6),IHR_S(:,4),'sum');
IHR_Stotal=horzcat(year,IHR_Stotal);
csvwrite('IHR_SprSum_Chinook.csv',IHR_Stotal) %writes to file as comma-separated values.

Fi=logical((m==8 & d>=12)+(m>=9)); %index the Fall run
IHR_F=IHR(Fi,:);
IHR_Ftotal=grpstats(IHR_F(:,5:6),IHR_F(:,4),'sum');
IHR_Ftotal=horzcat(year,IHR_Ftotal);
csvwrite('IHR_Fall_Chinook.csv',IHR_Ftotal) %writes to file as comma-separated values.

figure(2)
plot(year, IHR_Ftotal(:,2), '.k')



%% PDO data
fid = fopen('Users/alisoniles/Google Drive/Coordinated Recoveries/Data/raw/PDO_monthly.csv');
    [data] = textscan(fid,'%f%f%f','delimiter',',','headerlines',headerlines); %fclose(fid);
PDO=cell2mat(data);
PDOy=unique(PDO(:,1));
PDOa=grpstats(PDO(:,3),PDO(:,1),'mean');
PDOya=horzcat(PDOy,PDOa); %Yearly average PDO for export
csvwrite('PDOya.csv',PDOya) 


%% WPE Analysis

WPE_results=zeros(6,7); %Matrix to store WPE results
No_perm_results=zeros(6,7); %Matrix to store the number of permutations behind each WPE value

wl=3; %word length of permutation, i.e. the length of the permutation order of successive time points
for tau=1:6 %time delay between the 'letters' of each 'word'

[ WPE_results(tau,1), No_perm_results(tau,1) ] = WPE_gap_func(SRO_n, SRO_y, wl, tau); %WPE of Southern Resident Orcas
[ WPE_results(tau,2), No_perm_results(tau,2) ] = WPE_gap_func(CSL_n, CSL_y, wl, tau); %WPE of sea lion pup counts
[ WPE_results(tau,3), No_perm_results(tau,3) ] = WPE_gap_func((CHA), mtime, wl, tau); %WPE of Chinook salmon (adults) with daily counts
[ WPE_results(tau,4), No_perm_results(tau,4) ] = WPE_gap_func((CHAy), year, wl, tau); %WPE of Chinook salmon (adults) for annual sums
[ WPE_results(tau,5), No_perm_results(tau,5) ] = WPE_gap_func(IHR_Stotal(:,2)+IHR_Stotal(:,3), IHR_Stotal(:,1), wl, tau); %WPE of Chinook salmon (adults) with daily counts
[ WPE_results(tau,6), No_perm_results(tau,6) ] = WPE_gap_func(IHR_Ftotal(:,2)+IHR_Ftotal(:,3), IHR_Ftotal(:,1), wl, tau); %WPE of Chinook salmon (adults) for annual sums
[ WPE_results(tau,7), No_perm_results(tau,7) ] = WPE_gap_func(PDO(:,3), (1:length(PDO))', wl, tau); %WPE of Chinook salmon (adults) for annual sums

end

%% Visualize WPE Results
f=figure(3);
hold all
set(f, 'Position', [1000, 1000, 600, 600]);
set(f, 'DefaultAxesFontSize', 18);
set(f, 'DefaultLineLineWidth',1.5);
   

fp=plot((1:6),WPE_results);
legend(fp, 'orcas', 'CA sea lions', 'BON chinook - daily','BON chinook - annual', 'IHR chinook Spring/Summer','IHR chinook Fall','PDO', 'Location', 'SouthWest');
   legend('Boxoff');
ylabel('WPE')         
xlabel('tau')
set(gca, ...
            'XLim', [1, 6]    , ... 
            'XTick', 1:1:6      , ...
            'Linewidth', 1, ...
            'Box','on'     ) 

saveas(f,'WPE_Salmon_SeaLion_Orca','jpg')           


