% find the actual max voxels for the D_KL model
%
EXPT = contextExpt();
modeldir = EXPT.modeldir;
%V = spm_vol(fullfile(modeldir, ['model59'], ['subj1'], sprintf('beta_%04d.nii',1)));
%V = spm_vol(fullfile(modeldir, ['model53'], ['con1'], sprintf('con_%04d.nii',1)));
V = spm_vol(fullfile(modeldir, ['model53'], ['con1'], 'spmT_0001.nii'));
Y = spm_read_vols(V);

cor = mni2cor([34 -64 48],V.mat)
Y(cor(1), cor(2), cor(3)) % sanity check -- should be 8.8376


% obtained from Show Results Table from ccnl_view(contextExpt(), 53,
% 'surprise');
%
rois = {'Angular_R', 'Frontal_Inf_Oper_R', 'Frontal_Mid_2_R_dorsal', ...
        'Frontal_Mid_2_R_ventral', 'OFCant_R', 'Frontal_Mid_2_L', ...
        'Frontal_Sup_2_L', 'Frontal_Inf_Tri_L'};
max_voxels = {[34 -64 48], [48 20 34], [34 12 54], ...
              [36 56 -2],  [20 48 -16], [-42 56 2], ...
              [-24 60 -10], [-44 20 22]};
rois_max_voxels = containers.Map(rois, max_voxels);

for j = 1:numel(max_voxels)
    roi = rois{j};
    voxel = max_voxels{j};
    cor = mni2cor(voxel,V.mat);
    fprintf('%s --> %.4f %.4f %.4f = %.4f\n', roi, voxel(1), voxel(2), voxel(3), Y(cor(1), cor(2), cor(3)));
end
            
sss = getGoodSubjects();

%{
n_trials_per_run = 20;
EXPT = contextExpt();
prev_trials = zeros(numel(sss) * 9 * 19, numel(max_voxels)); % col = ROI, row = activation of max voxel for given trial
idx = 0;
for subj = sss
    modeldir = fullfile(EXPT.modeldir,['model59'],['subj',num2str(subj)]);
    fprintf('subj = %d\n', subj);
    for run = 1:9
        fprintf('     run = %d\n', run);
        for i = (1:19) + (run - 1) * (n_trials_per_run + 6)
            % for each feedback onset on trials 1..19
            %
            V = spm_vol(fullfile(modeldir, sprintf('beta_%04d.nii', i)));
            Y = spm_read_vols(V);
            idx = idx + 1;
            % for each ROI, get the activation at the max voxel
            %
            for j = 1:numel(max_voxels)
                voxel = max_voxels{j};
                cor = mni2cor(voxel, V.mat);
                value = Y(cor(1), cor(2), cor(3));
                prev_trials(idx, j) = value;
            end
        end
    end
end
prev_trials_act = prev_trials;

analyze;

%}


% load('DKL.mat');
%bla = prev_trials_surprise > 0.2;
%scatter(prev_trials_surprise(bla), prev_trials_act(bla, 6));
