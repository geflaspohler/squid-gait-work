function ReadRaw(recdir, prefix, df,save_flag)

% ReadRaw(recdir, prefix, df,save_flag)
% recdir = 'L:\DQ_All_Day_Data';
% prefix = 'tt13_271c';
% df = 10;
if nargin<3,
    save_flag = 0;
end

deploy_name = prefix;
filename = sprintf('%s_Data', prefix);

% Notify Matlab of where processed data should be stored
startup(recdir);
% Read the sensor data:
X = d3readswv(recdir,prefix,df);
X_ORI = X;
% [ch_names,descr,ch_nums,type] = d3channames(X.cn);
% Register the deployment:
[CAL,~] = d3deployment(recdir,prefix,prefix) ;
CAL_ORI = CAL;

%% 
p = d3calpressure(X, CAL);
p_ORI = p;
A = d3calacc(X,CAL);
A_ORI = A;
M = d3calmag(X,CAL);
M_ORI = M;

%%
% cd('C:\Users\Yunli Shao\Desktop\Dropbox\SpeedSensor\DTAG\Data');
% cd('.\Data')
if save_flag
    save(['.\Data\', filename], 'X', 'X_ORI', 'CAL', 'CAL_ORI', 'A_ORI', 'M_ORI', 'p_ORI', 'df')
end

fprintf('\nDone\n')