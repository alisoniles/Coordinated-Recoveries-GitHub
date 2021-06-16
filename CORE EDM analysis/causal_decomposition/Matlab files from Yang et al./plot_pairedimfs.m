% Function for plotting and comparisons of intrinsic mode functions (IMFs) 
% of a pair of time series
%
% Syntax: plot_pairedimfs(t,s1,s2,rnoise,nensemble,imfindex)
% Parameters:
% t:  time vector
% s1: time series 1
% s2: time series 2
% rnoise: the level of noise used in ensemble empirical mode decomposition
%         the rnoise value represents the fraction of standard deviation of
%         time series (e.g., 0.1)
% nensemble: the number of ensemble to average the results of noise-assisted
%            empirical mode decomposition. 
%            The error of IMFs equals to rnoise/sqrt(nensemble)        
% imfindex: selection of IMFs to be shown in the figures (default: []; all IMFs)
% 
% Example: load ecosystem_data
%          plotimf(time_DIDIPARA,DIDI,PARA,0.35,1000,[]);
%
% Ver 1.0: Albert C. Yang, MD, PhD 7/11/2018
%
% Referece

function  output = plot_pairedimfs(t,s1,s2,rnoise,nensemble,imfindex)

% normalization
s1r=s1-mean(s1);
s1r=s1r/std(s1);
s2r=s2-mean(s2);
s2r=s2r/std(s2);

% do EMD for each signal
if rnoise>0
    c1=eemd(s1r,rnoise,nensemble,0)';
    c2=eemd(s2r,rnoise,nensemble,0)';
else
    nensemble=1;
    c1=eemd(s1r,rnoise,nensemble,0)';
    c2=eemd(s2r,rnoise,nensemble,0)';
    c1=c1(:,2:end);
    c2=c2(:,2:end);
end

% determine IMF size
imfsize=size(c1);
imfsize=imfsize(2);

% determine time data
tn=length(t);
t1=min(t(1),t(tn));
t2=max(t(1),t(tn));

if isempty(imfindex)
    imfindex=1:1:imfsize;
    rowsize=imfsize+1;
else
    rowsize=length(imfindex)+1;
end

% plot raw data
subplot(rowsize,5,(1:1:5));
plot(t,s1r,'Color','b','LineWidth',3);
hold on
plot(t,s2r,'Color','r','LineWidth',3);
hold off
set(gca,'FontSize',14)
minc1=min(s1r);minc2=min(s2r);
maxc1=max(s1r);maxc2=max(s2r);
axis([t1 t2 min(minc1,minc2)*0.8 max(maxc1,maxc2)*1.2]);	
ylabel('Normalized Unit');

% plot IMFs and trend
for j=1:length(imfindex);
   subplot(rowsize,5,(1:1:5)+j*5);
   plot(t,c1(:,imfindex(j)),'Color','b','LineWidth',3);
   hold on
   plot(t,c2(:,imfindex(j)),'Color','r','LineWidth',3);
   hold off
   set(gca,'FontSize',14)
   if imfindex(j)==imfsize
       minc1=min(c1(:,imfindex(j)));
       minc2=min(c2(:,imfindex(j)));
       maxc1=max(c1(:,imfindex(j)));
       maxc2=max(c2(:,imfindex(j)));
       axis([t1 t2 min(minc1,minc2)*0.8 max(maxc1,maxc2)*1.2]);
       imflabel='Trend';
   else
       minc1=min(c1(:,imfindex(j)));
       minc2=min(c2(:,imfindex(j)));
       maxc1=max(c1(:,imfindex(j)));
       maxc2=max(c2(:,imfindex(j)));
       yminmax=(max([abs(minc1),abs(minc2),abs(maxc1),abs(maxc2)]))*1.2;
       axis([t1 t2 -yminmax yminmax]);  
       imflabel='IMF ';
       imflabel=[imflabel num2str(imfindex(j))];
   end
   ylabel(imflabel); 
end

 
