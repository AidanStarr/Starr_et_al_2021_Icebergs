%%% script: analysis 4_b (to create output files from analysis_4
%% Split into intervals
enwx = 1:1:100';
for i = 1:23
    newird_cum(:,i) = interp1(time_cum_norm{i},ird_cum{i},enwx,'linear');
end

for i = 1:23
    newtime_cum(:,i) = interp1(time_cum_norm{i},time_cum_norm{i},enwx,'linear');
end

for i = 1:23
    newd18o_cum(:,i) = interp1(time_cum_norm{i},d18o_cum{i},enwx,'linear');
end

% mean accumulation for IRD 
A = newird_cum(:,1:12);
B = newird_cum(:,13:18);
C = newird_cum(:,19:23);
pre_mpt = mean(A,2);
mpt = mean(B,2);
post_mpt = mean(C,2);


% calculate confidence intervals for IRD 
SEM_A = std(A,1,2)/sqrt(40); % get standard error of the mean
SEM_B = std(B,1,2)/sqrt(40);
SEM_C = std(C,1,2)/sqrt(40);
CI95_A = bsxfun(@plus, mean(A,2), bsxfun(@times, [-1  1]*1.96, SEM_A));   % 95% Bootstrap Confidence Intervals
CI95_B = bsxfun(@plus, mean(B,2), bsxfun(@times, [-1  1]*1.96, SEM_B));   
CI95_C = bsxfun(@plus, mean(C,2), bsxfun(@times, [-1  1]*1.96, SEM_C));   % 95% Confidence Intervals

% mean accumulation for d18O
A2 = newd18o_cum(:,1:11);
B2 = newd18o_cum(:,2:18);
C2 = newd18o_cum(:,19:23);
pre_mptd = mean(A2,2);
mptd = mean(B2,2);
post_mptd = mean(C2,2);

% calculate confidence intervals for d18O 
SEM_A2 = std(newd18o_cum(:,1:11),1,2)/sqrt(40);
SEM_B2 = std(newd18o_cum(:,12:18),1,2)/sqrt(40);
SEM_C2 = std(newd18o_cum(:,19:23),1,2)/sqrt(40);
CI95_A2 = bsxfun(@plus, mean(A2,2), bsxfun(@times, [-1  1]*1.96, SEM_A2));   % 95% Confidence Intervals
CI95_B2 = bsxfun(@plus, mean(B2,2), bsxfun(@times, [-1  1]*1.96, SEM_B2));   % 95% Confidence Intervals
CI95_C2 = bsxfun(@plus, mean(C2,2), bsxfun(@times, [-1  1]*1.96, SEM_C2));   % 95% Confidence Intervals


figure
subplot(131);
plot(enwx',pre_mpt,'r');
hold on
plot(enwx',CI95_A(:,1),'r','linestyle','--');
plot(enwx',CI95_A(:,2),'r','linestyle','--');
plot(enwx',pre_mptd,'b');
plot(enwx',CI95_A2(:,1),'b','linestyle','--');
plot(enwx',CI95_A2(:,2),'b','linestyle','--');
ylim([0 100])
set(gca,'xdir','reverse')

subplot(132);
plot(enwx',mpt,'r');
hold on
plot(enwx',CI95_B(:,1),'r','linestyle','--');
plot(enwx',CI95_B(:,2),'r','linestyle','--');
plot(enwx',mptd,'b');
plot(enwx',CI95_B2(:,1),'b','linestyle','--');
plot(enwx',CI95_B2(:,2),'b','linestyle','--');
ylim([0 100])
set(gca,'xdir','reverse')

subplot(133);
plot(enwx',post_mpt,'r');
hold on
plot(enwx',CI95_C(:,1),'r','linestyle','--');
plot(enwx',CI95_C(:,2),'r','linestyle','--');
plot(enwx',post_mptd,'b');
plot(enwx',CI95_C2(:,1),'b','linestyle','--');
plot(enwx',CI95_C2(:,2),'b','linestyle','--');
ylim([0 100])
set(gca,'xdir','reverse')

%%%% output %%%%
output_cumtrap(:,1) = pre_mpt;
output_cumtrap(:,2) = mpt;
output_cumtrap(:,3) = post_mpt;

output_cumtrap(:,4) = CI95_A(:,1);
output_cumtrap(:,5) =  CI95_A(:,2);
output_cumtrap(:,6) =  CI95_B(:,1);
output_cumtrap(:,7) =  CI95_B(:,2);
output_cumtrap(:,8) =  CI95_C(:,1);
output_cumtrap(:,9) =  CI95_C(:,2);

output_cumtrap(:,10) = pre_mptd;
output_cumtrap(:,11) = mptd;
output_cumtrap(:,12) = post_mptd;

output_cumtrap(:,13) = CI95_A2(:,1);
output_cumtrap(:,14) =  CI95_A2(:,2);
output_cumtrap(:,15) =  CI95_B2(:,1);
output_cumtrap(:,16) =  CI95_B2(:,2);
output_cumtrap(:,17) =  CI95_C2(:,1);
output_cumtrap(:,18) =  CI95_C2(:,2);

output_cumtrap(:,19) = enwx';
output_4 = array2table(output_cumtrap,'variablenames',{'IRD_mean_pre','IRD_mean_mpt','IRD_mean_post','IRDu_pre','IRDl_pre','IRDu_mpt','IRDl_mpt','IRDu_post','IRDl_post','d18o_mean_pre','d18o_mean_mpt','d18o_mean_post','d18ou_pre','d18ol_pre','d18ou_mpt','d18ol_mpt','d18ou_post','d18ol_post','time'});
writetable(output_4,'outputs/cumulative_glacials.csv')
