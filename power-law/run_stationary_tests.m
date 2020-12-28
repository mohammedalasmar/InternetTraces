clc, clear all
close all

%   eval(sprintf('trace_%d = GoodPutVsByteSorted',changeable_var(ii) ))

%   A = importdata('BC-pAug89.txt');
% all_time_points = A(:,1);
% all_data_points = A(:,2);
% load('all_twente_traces_time_vs_bytes_RANKED.mat');

% matObj = matfile('all_twente_traces_time_vs_bytes_RANKED.mat');
% details = whos(matObj)
% length(details)

for ii = 1:1
    
    trace_name = ['trace_', num2str(ii)];
    load('all_twente_traces_time_vs_bytes_RANKED.mat' , trace_name);
    trace = eval(trace_name);
    
    all_time_points = trace(:,1);
    all_data_points = trace(:,2);
%       histogram(dataRate)
    %%
    
    time_values = all_time_points;
    packet_length =   all_data_points;  % bits
    
    %  time_scales_T = [0.005, 0.01,0.05,0.1, 0.5, 1, 2, 4,8,12,20,25,30,50,60];
%     
%     time_scales_T = [ 0.1 , 0.5, 1, 2, 4, 8, 10];
    time_scales_T = [0.01];
    
    for t_i=1:length(time_scales_T)
        
        t_i
        aggregationTime = time_scales_T(t_i);
        % getting data rates from the input trace
        picked_times_index=0;
        dataRate = 0;
        aggData = 0;
        j=1;
        T = aggregationTime;
        % aggregation based on selected aggregationTime
        for i=1:length(time_values)
            if time_values(i)>=T
                picked_times_index(j)=i;
                T = T + aggregationTime;
                j=j+1;
            end
        end
        %
        % picked_times_index=[1 picked_times_index length(time_values)+1];
        % aggreate packets
        for u=1:length(picked_times_index)-1
            aggData(u)=sum(packet_length(picked_times_index(u):picked_times_index(u+1)-1));
        end
        
        T = aggregationTime;
        t=0:T:time_values(end);
        aaggDataLen= length(aggData);
        aggData(aaggDataLen:length(t))=  aggData(aaggDataLen);
        % data rate
        dataRate=8*aggData/T; % bits/sec
        
        len_data_rate(t_i) = length(dataRate);
        
        %     if t_i == 1
        %         figure, histogram(dataRate);
        %     end
         
             figure , subplot(121), plot(t,dataRate, '-r', 'LineWidth', 2), set(gca, 'FontSize',28), set(gca,'FontName','Times'),  grid on, box on;
        %     xlabel('Time (seconds)'), ylabel(' Data rate (bps)'), imtitl = [' time scale T = ', num2str(T) , ' sec'];  title(imtitl)
        
        
        %    lags = [0:sqrt(length(dataRate))]';
        
         lenn = length(dataRate);
%     lags = floor(sqrt(lenn));
    lags_adf = floor(12*(lenn/100)^(1/4));
    lags_kspp = floor(lenn^(1/2));
%     trunc(12*(n/100)^0.25)
    
    lags_pp = floor(lenn^(1/4));
%         lags = floor(sqrt(length(dataRate)));
        
%            subplot(122), 
%            
 figure, autocorr(dataRate,'NumLags',90000), set(gca, 'FontSize',28) ,  set(gca, 'FontSize',32)
     diff_data = diff(dataRate,1);  set(gca, 'FontSize',32)
     figure, histogram(dataRate,200);  set(gca, 'FontSize',32)
          figure, plot(aggData); ylabel('dataRate'), xlabel('time (sec)'), set(gca, 'FontSize',32)
 figure, autocorr(diff_data,'NumLags',50), set(gca, 'FontSize',28) ,  set(gca, 'FontSize',32)

     figure, plot(diff_data); ylabel('diff. data'), xlabel('time (sec)'), set(gca, 'FontSize',32)
        [h_kpss,pValue_kpss , stat_kpss,cValue_kpss] =kpsstest(dataRate, 'trend' , false , 'alpha',0.05, 'lags' , lags_kspp) % 'alpha' significance levels for the hypothesis tests. 95% confidence interval
%         [h_lmc,pValue_lmc, stat_lmc,cValue_lmc] = lmctest(dataRate, 'trend' , false , 'alpha',0.05, 'lags' , lags_kspp , 'test', 'var2') 
        [h_lmc,pValue_lmc, stat_lmc,cValue_lmc] = kpsstest(dataRate, 'trend' , false , 'alpha',0.05, 'lags' ,lags_kspp) 
        
%         sh = 10;
%         x = 1:10:3000;
%        figure, plot(x(sh:end) , stat_kpss(sh:end),'-r'), hold on , plot(x(sh:end), cValue_kpss(sh:end),'.-b'), hold on , plot(x(sh:end), pValue_kpss(sh:end),'.--m'), 
%        hold on , plot(0.05*ones(1,3000),'--g')
%         legendInfo{1} = 'stat kpss', legendInfo{2} = 'cValue kpss', legendInfo{3} = 'pValue kpss'
%   lgd = legend;
%   legend(legendInfo)
%   grid on, xlabel('lag')
        
        %         [h_adf,pValue_adf , stat_adf,cValue_adf] = adftest(dataRate, 'alpha',0.05 ,'model',{'ARD'}, 'lags' , 15 , 'test','t2')
%         k = 1000
%         figure, plot(aggData(1:end-k), aggData(k+1:end), 'o')
%         tittt = [' lag k = ', num2str(k)]
%         title(tittt)
%         set(gca, 'FontSize',32)
%         
%         figure, plot((diff_data(1:end-k)), (diff_data(k+1:end)), 'o')
%         tittt = [' lag k = ', num2str(k)]
%         title(tittt)
%         set(gca, 'FontSize',32)
        
        
        [h_adf,pValue_adf , stat_adf,cValue_adf] = adftest(dataRate, 'alpha',0.05 ,'model',{'ARD'}, 'lags' , lags_adf  ) 

        [h_pp,pValue_pp , stat_pp ,cValue_pp ] = pptest(dataRate, 'alpha',0.05, 'model',{'ARD'},'lags' , lags_pp) 
        [h_vratio,pValue_vratio] = vratiotest(dataRate, 'alpha',0.05);
        
        h_kpss_all(t_i) = h_kpss ;pValue_kpss_all(t_i)= pValue_kpss;
        h_lmc_all(t_i) = h_lmc; pValue_lmc_all(t_i)= pValue_lmc;
        h_adf_all(t_i) = h_adf; pValue_adf_all(t_i)= pValue_adf;
        h_pp_all(t_i) = h_pp; pValue_pp_all(t_i)= pValue_pp;
        h_vratio_all(t_i) = h_vratio; pValue_vratio_all(t_i)= pValue_vratio;
        
        
    end
    
    
    for i = 1:length(time_scales_T)
        if h_kpss_all(i) == false && pValue_kpss_all(i)>=0.05
            kpss_results{i} ='S';
        else
            kpss_results{i} ='Non';
        end
    end
    
    
    for i = 1:length(time_scales_T)
        if h_lmc_all(i) == false && pValue_lmc_all(i)>=0.05
            lmc_results{i} ='S';
        else
            lmc_results{i} ='Non';
        end
    end
    
    for i = 1:length(time_scales_T)
        if h_adf_all(i) == true && pValue_adf_all(i)<=0.05
            adf_results{i} ='S';
        else
            adf_results{i} ='Non';
        end
    end
    
    for i = 1:length(time_scales_T)
        if h_pp_all(i) == true && pValue_pp_all(i)<=0.05
            pp_results{i} ='S';
        else
            pp_results{i} ='Non';
        end
    end
    
    for i = 1:length(time_scales_T)
        if h_vratio_all(i) == true && pValue_vratio_all(i)<=0.05
            vratio_results{i} ='S';
        else
            vratio_results{i} ='Non';
        end
    end
    
%     a = [kpss_results;lmc_results;adf_results;pp_results]
    
   a = [kpss_results;adf_results;pp_results]

%     len_data_rate
    
%       resu_mat_name = [trace_name,'_res.mat']
%       save(resu_mat_name, 'a', 'time_scales_T', 'len_data_rate')
    
end


% for i = 1:1
%     
%     
%     resu_mat_name = ['trace_',num2str(i), '_res.mat']
% 
%     load(resu_mat_name)
%     b{i} = a 
%     len{i} = len_data_rate;
% 
% end


  

% legendInfo{1} = 'flow 1'
% lgd = legend;
% lgd.NumColumns = 2;
% legend(legendInfo)
% legend( legendInfo ,'FontSize',27)
% ax = gca;
% ax.YAxis.Exponent = 0;
% ax.XAxis.Exponent = 0;
% set(gca, 'FontSize',32)
% savefig('src.fig');
% h1 =openfig('src.fig');
% set(h1,'Units','Inches');
% pos = get(h1,'Position');
% set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% saveas(h1,'src.pdf');