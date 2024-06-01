%this is a modified version of the script provided by professor Lanconelli
%for the Physics of Medical imaging course labs A.Y: 23/24, Applied
%Physics, UNIBO
clear all
directory_path = uigetdir();
numpart = 10000000;
if ~isequal(directory_path, 0)  
    % Get the list of files in the selected directory
    file_list = dir(fullfile(directory_path, '*.hsb'));  

    % Initialize arrays to store p-values and PSFs
    all_p_values_Light = zeros(10,1);
    all_p_values_NoLight = zeros(10,1);
    num_detected = zeros(10,1);
    all_sens =  zeros(10,1);
    all_PSF_Light = zeros(10,1);
    all_PSF_NoLight = zeros(10,1);
    all_PSF_Light_ANI = NaN(10,2); % (sigmax,sigmay)
    all_PSF_NoLight_ANI = NaN(10,2); % (sigmax,sigmay)
for i = 1:numel(file_list)
    % Open the file
    filename = file_list(i).name;
    fid = fopen(filename);

    % Read and reshape data 
    [A,count] = fread(fid,'*single','n');
    fclose(fid);
    B = reshape(A,15,[]);
    ntuple = B';
% set the labels for the output variables 
     ntuple_label(1)="source_num";
     ntuple_label(2)="history";
     ntuple_label(3)="xorig";
     ntuple_label(4)="yorig";
     ntuple_label(5)="zorig";
     ntuple_label(6)="a1orig";
     ntuple_label(7)="a2orig";
     ntuple_label(8)="a3orig";
     ntuple_label(9)="eslab";
     ntuple_label(10)="xslab";
     ntuple_label(11)="yslab";
     ntuple_label(12)="zslab";
     ntuple_label(13)="eslabesr";
     ntuple_label(14)="xslabesr";
     ntuple_label(15)="yslabesr";

% do the histogram of one variable from the ntuple
%    for i=1:15
%        histogram(ntuple(:,i));
%        title(ntuple_label(i));
%        pause;
%    end
% fit the histogram with a gaussian
     %histfit(ntuple(:,14));
     %histfit(ntuple(:,15));
     %histfit(ntuple(:,10));
     %histfit(ntuple(:,11));
% get the parameters of the fitting curve
    % Compute and store p-values and PSFs for LIGHT
     num_detected(i,1) = size(ntuple,1);
     all_sens(i,1) = size(ntuple,1)/numpart;
     cols = [14 15];
     pd = fitdist(ntuple(:,cols(1)),'Normal');
     pd1 = fitdist(ntuple(:,cols(2)),'Normal');
     C = [ntuple(:,cols(1)), ntuple(:,cols(2))];
     p = vartestn(C,'Display','off');
     all_p_values_Light(i,1) = p;  

     if p > 0.05
        PSF = 2.3548 * pd.std;
        all_PSF_Light(i,1) = PSF;
     else
        all_PSF_Light(i,1) = NaN;
        PSFx = 2.3548 * pd.std;
        PSFy = 2.3548 * pd1.std;
        all_PSF_Light_ANI(i,1) = PSFx;
        all_PSF_Light_ANI(i,2) = PSFy;
     end

% Compute and store p-values and PSFs for NO LIGHT
     cols = [10 11];
     pd = fitdist(ntuple(:,cols(1)),'Normal');
     pd1 = fitdist(ntuple(:,cols(2)),'Normal');
     C = [ntuple(:,cols(1)), ntuple(:,cols(2))];
     p = vartestn(C,'Display','off');
     all_p_values_NoLight(i,1) = p; 

     if p > 0.05
        PSF = 2.3548 * pd.std;
        all_PSF_NoLight(i,1) = PSF;
     else
        all_PSF_Light(i,1) = NaN;
        PSFx = 2.3548 * pd.std;
        PSFy = 2.3548 * pd1.std;
        all_PSF_NoLight_ANI(i,1) = PSFx;
        all_PSF_NoLight_ANI(i,2) = PSFy;
     end

end
% Save p-values and PSFs to a file in the same directory
save('computed_results.mat', 'all_p_values_Light','all_p_values_NoLight', 'all_PSF_Light', 'all_PSF_NoLight','num_detected','all_sens','all_PSF_Light_ANI','all_PSF_NoLight_ANI');
all_PSF_NoLight = rmmissing(all_PSF_NoLight);
all_PSF_Light = rmmissing(all_PSF_Light);
all_PSF_NoLight_ANI = rmmissing(all_PSF_NoLight_ANI);
all_PSF_Light_ANI = rmmissing(all_PSF_Light_ANI);
mean_PSF_Light_ANI = zeros(1,2);
std_PSF_Light_ANI = zeros(1,2);
mean_PSF_NoLight_ANI = zeros(1,2);
std_PSF_NoLight_ANI = zeros(1,2);
mean_PSF_NoLight = mean(all_PSF_NoLight);
std_PSF_NoLight = std(all_PSF_NoLight);
mean_PSF_Light = mean(all_PSF_Light);
std_PSF_Light = std(all_PSF_Light);
mean_PSF_NoLight_ANI(1,1) = mean(all_PSF_NoLight_ANI(:,1));
std_PSF_NoLight_ANI(1,1) = std(all_PSF_NoLight_ANI(:,1));
mean_PSF_NoLight_ANI(1,2) = mean(all_PSF_NoLight_ANI(:,2));
std_PSF_NoLight_ANI(1,2) = std(all_PSF_NoLight_ANI(:,2));
mean_PSF_Light_ANI(1,1) = mean(all_PSF_Light_ANI(:,1));
std_PSF_Light_ANI(1,1) = std(all_PSF_Light_ANI(:,1));
mean_PSF_Light_ANI(1,2) = mean(all_PSF_Light_ANI(:,2));
std_PSF_Light_ANI(1,2) = std(all_PSF_Light_ANI(:,2));
mean_numdetected = mean(num_detected);
std_numdetected = std(num_detected);
mean_pvalue_Light = mean(all_p_values_Light);
std_pvalue_Light = std(all_p_values_Light);
mean_pvalue_NoLight = mean(all_p_values_NoLight);
std_pvalue_NoLight = std(all_p_values_NoLight);

save('mean_and_stds.mat', 'mean_PSF_NoLight', 'std_PSF_NoLight', 'mean_PSF_Light', ...
    'std_PSF_Light', 'mean_PSF_NoLight_ANI', 'std_PSF_NoLight_ANI', ...
    'mean_PSF_Light_ANI', 'std_PSF_Light_ANI', 'mean_numdetected', ...
    'std_numdetected', 'mean_pvalue_Light', 'std_pvalue_Light', ...
    'mean_pvalue_NoLight', 'std_pvalue_NoLight');
% 2D scatter plot
end
%scatter(ntuple(:,10),ntuple(:,11),'.')
% 3D scatter plot
%scatter3(ntuple(:,10),ntuple(:,11),ntuple(:,12),'.')
% 