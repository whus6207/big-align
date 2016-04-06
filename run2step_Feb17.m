function [ ] = run2step_Feb17( nodesNo, evalsNo )
%% 


filenameEnd = {'_.mat', ... 
    '_10percent_noise.mat', '_20percent_noise.mat', ...
    '_30percent_noise.mat'} 
noise = [ 0, 5, 10, 20, 30]; 

folder = sprintf('movie_genre_%dn_uniqueSum', nodesNo)
filename = sprintf('movie_genre_%dn_unique', nodesNo)

epsilon = 0.00001; 
lamdaVector =  [ 1, 0.5, 0.1, 0.01, 0.001, 0 ];
runs = 10;

if nargin < 2
    evalsNo = 6;
end

%% run the algorithm 10 times per noise level + lamda, and report the avg
for i = 1 : length(filenameEnd)
    for eps = 1 : length(lamdaVector)
        eps1 = 2*eps-1;
        eps2 = 2*eps;
        for r = 1 : runs
            name = sprintf('data_mat/%s/%s_bundle_%d%s', folder, filename, r, filenameEnd{i})
            [ accPhung(i, eps1), accQhung(i, eps1), accTotalHung(i, eps1), ...
              accP(i, eps1), accQ(i, eps1), accTotal(i, eps1), ...
              matchingPercentQhungInit(i, eps1), matchingPercentQInit(i, eps1), elapsedTime(i, eps1) ] = ...
                twoStepApproach_Feb17('bundle', name, 10000, lamdaVector(eps), noise(i), evalsNo, 1);
            accP_perRun(r) = accP(i,eps1);
            accQ_perRun(r) = accQ(i,eps1);
            time_perRun(r) = elapsedTime(i,eps1);
            accTotal_perRun(r) = accTotal(i,eps1);
            accPhung_perRun(r) = accPhung(i,eps1);
            accQhung_perRun(r) = accQhung(i,eps1);
            accTotalhung_perRun(r) = accTotalHung(i,eps1);
            accQinit_hung_perRun(r) = matchingPercentQhungInit(i, eps1);
            accQinit_perRun(r) = matchingPercentQInit(i, eps1);
        end
        accP(i, eps1) = mean(accP_perRun);
        accP(i,eps2) = std(accP_perRun);
        accQ(i, eps1) = mean(accQ_perRun);
        accQ(i,eps2) = std(accQ_perRun);
        accTotal(i, eps1) = mean(accTotal_perRun);
        accTotal(i,eps2) = std(accTotal_perRun);
        elapsedTime(i, eps1) = mean(time_perRun);
        elapsedTime(i, eps2) = std(time_perRun);
        accPhung(i, eps1) = mean(accPhung_perRun);
        accPhung(i,eps2) = std(accPhung_perRun);
        accQhung(i, eps1) = mean(accQhung_perRun);
        accQhung(i,eps2) = std(accQhung_perRun);
        accTotalhung(i, eps1) = mean(accTotalhung_perRun);
        accTotalhung(i,eps2) = std(accTotalhung_perRun);
        matchingPercentQhungInit(i, eps1) = mean(accQinit_hung_perRun);
        matchingPercentQhungInit(i,eps2) = std(accQinit_hung_perRun);
        matchingPercentQInit(i, eps1) = mean(accQinit_perRun);
        matchingPercentQInit(i,eps2) = std(accQinit_perRun);
    end
    eps2 = eps2 + 1;
    [~, idx] = max(accP(i,1:(eps2-1)), [], 2);
    accP(i, eps2) = idx;
    [~, idx ] = max(accQ(i,1:(eps2-1)), [], 2);
    accQ(i, eps2) = idx;
    [~, accTotal(i, eps2)] = max(accTotal(i,1:(eps2-1)), [], 2);
end


output = sprintf('01twoStep_nodes_%d_lamda_VAR_Feb17.mat', nodesNo);
save(output , ...
    'accP', 'accQ', 'accTotal',  ...
    'elapsedTime')

output2 = sprintf('01twoStep_ACCp_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %d\n', [accP]');
fclose(fid);

output2 = sprintf('01twoStep_ACCq_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %d\n', [ accQ]');
fclose(fid);

output2 = sprintf('01twoStep_ACCtOtal_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %d\n', [accTotal]');
fclose(fid);

output2 = sprintf('01twoStep_ACCpHung_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f\n', [accPhung]');
fclose(fid);

output2 = sprintf('01twoStep_ACCqHung_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f\n', [accQhung]');
fclose(fid);

output2 = sprintf('01twoStep_ACCtOtalHung_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f\n', [accTotalhung]');
fclose(fid);

output2 = sprintf('01twoStep_ACCqINITHung_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f\n', [matchingPercentQhungInit]');
fclose(fid);

output2 = sprintf('01twoStep_ACCqINIT_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f\n', [matchingPercentQInit]');
fclose(fid);

output2 = sprintf('01twoStep_Runtime_nodes_%d_lamda_VAR_Feb17.txt', nodesNo);
fid = fopen(output2,'w');
fprintf(fid, '%.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f %.3f %.6f \n', elapsedTime');
fclose(fid);


end
