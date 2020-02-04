% Load Data for Starr et al 2019 'A Southern Ocean Lead'
data = xlsread('starr_2019_data.xlsx');
ag = data(:,2);
[age,n]=unique(ag);
depth = data(n,1);
bd18o = data(n,3);
bd13c = data(n,4);
ird = data(n,5);
irdmar = data(n,6);
clear('ag');


