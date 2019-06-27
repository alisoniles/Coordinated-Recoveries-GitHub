close all
clearvars

addpath(genpath('/Users/alisoniles/Google Drive/Coordinated Recoveries/CORE_EDM/Data/Raw/'))

PDOmonthly = csvread('PDO_monthly.csv', 1, 0);
PDO = zeros(72,12);
for i=1117:size(PDOmonthly,1) %from 1947 to 2018
    PDO(PDOmonthly(i,1)-1946,PDOmonthly(i,2))=PDOmonthly(i,3);
end

NPGOmonthly = csvread('NPGO_monthly.csv', 1, 0);
NPGO = zeros(69,12);
for i=1:size(NPGOmonthly,1) %from 1950 to 2018
    NPGO(NPGOmonthly(i,1)-1949,NPGOmonthly(i,2))=NPGOmonthly(i,3);
end