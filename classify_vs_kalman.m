%
%
% One thing that would be useful to know is whether the hippocampal classifier is related to 
% how consistent test trial behavior is with a modulatory structure. What you can do is generate
% test trial predictions for each single-structure model, and compute the likelihood of the test
% trials for each of these structures, normalized across the different structures. 
% Then you can take the normalized likelihood for the modulatory model (which expresses how 
% "modulatory" the behavior looks on each block) and correlate it with the probability of the 
% modulatory class from the hippocampal classifier on the (held-out) last training trial of 
% each block. So for each subject you'll compute a correlation (across blocks) and then report
% the average correlation. This tells us how much the hippocampal beliefs about structure 
% influence predict subsequent behavior. You can do the same thing with other ROIs and for 
% other classes (additive, irrelevant) but start with the hippocampus.

% obtained from classify_scp_from_ncf.m
%
close all;
clear all;
rois = {'hippocampus', 'ofc', 'striatum', 'vmpfc', 'rlpfc', 'bg', 'pallidum'};

held_out_trials = [18 19 20];

folders = {
    'classify_heldout_trial_1', ...
    'classify_heldout_trial_2', ...
    'classify_heldout_trial_3', ...
    'classify_heldout_trial_4', ...
    'classify_heldout_trial_5', ...
    'classify_heldout_trial_6', ...
    'classify_heldout_trial_7', ...
    'classify_heldout_trial_8', ...
    'classify_heldout_trial_9', ...
    'classify_heldout_trial_10', ...
    'classify_heldout_trial_11', ...
    'classify_heldout_trial_12', ...
    'classify_heldout_trial_13', ...
    'classify_heldout_trial_14', ...
    'classify_heldout_trial_15', ...
    'classify_heldout_trial_16', ...
    'classify_heldout_trial_17', ...
    'classify_heldout_trial_18', ...
    'classify_heldout_trial_19', ...
    'classify_heldout_trial_20'};


% held out trial 1
train_files = [
    {'classify_train_cvglmnet_hippocampus_condition_ZAGAABEXLM.mat', ...
    'classify_train_cvglmnet_ofc_condition_NTGYFJTHFC.mat', ...
    'classify_train_cvglmnet_striatum_condition_RRVUXGTXGW.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_KNXPEFSXTH.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_DNOSVOODWD.mat', ...
    'classify_train_cvglmnet_bg_condition_UHIVBSAPBW.mat', ...
    'classify_train_cvglmnet_pallidum_condition_SXMSPRVPEE.mat'}];
% held out trial 2
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_HMYTSSARRH.mat', ...
    'classify_train_cvglmnet_ofc_condition_JZZYRXUIGV.mat', ...
    'classify_train_cvglmnet_striatum_condition_ZDBEJKHSHE.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_RMQIWKZAMB.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_FCDOEBYQXK.mat', ...
    'classify_train_cvglmnet_bg_condition_HGSVSLRONT.mat', ...
    'classify_train_cvglmnet_pallidum_condition_QRYFRZFWSM.mat'}];
% held out trial 3
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_WIURLWSOPP.mat', ...
    'classify_train_cvglmnet_ofc_condition_ZQUCOJIWHZ.mat', ...
    'classify_train_cvglmnet_striatum_condition_XAQJJOSNSU.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_TWINLPYCPS.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_LJGZRZDWAQ.mat', ...
    'classify_train_cvglmnet_bg_condition_FSVIHJFTQK.mat', ...
    'classify_train_cvglmnet_pallidum_condition_GHROKGPEJW.mat'}];
% held out trial 4
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_QYKIIEGRPK.mat', ...
    'classify_train_cvglmnet_ofc_condition_RKZDJSKAII.mat', ...
    'classify_train_cvglmnet_striatum_condition_YESBSFUZKZ.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_TLNZHTBUMN.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_LHQRDOXQQN.mat', ...
    'classify_train_cvglmnet_bg_condition_IPWSKKQRIC.mat', ...
    'classify_train_cvglmnet_pallidum_condition_NJCWWAPJUF.mat'}];
% held out trial 5
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_VGJXGHHFSL.mat', ...
    'classify_train_cvglmnet_ofc_condition_KDQAHXOQWI.mat', ...
    'classify_train_cvglmnet_striatum_condition_FCFJDEFAVG.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_PCYEUIVYZW.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_DEUEEXQTVE.mat', ...
    'classify_train_cvglmnet_bg_condition_ULEQJPJFLY.mat', ...
    'classify_train_cvglmnet_pallidum_condition_TLJWNGQCQZ.mat'}];
% held out trial 6
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_UMMMITCTYA.mat', ...
    'classify_train_cvglmnet_ofc_condition_NFAMGFMSQV.mat', ...
    'classify_train_cvglmnet_striatum_condition_ABVVLUXJUJ.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_TMUBJJSXHB.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_EGIIZGLLJN.mat', ...
    'classify_train_cvglmnet_bg_condition_UWKCZONEXW.mat', ...
    'classify_train_cvglmnet_pallidum_condition_GHWVTCDJNB.mat'}];
% held out trial 7
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_WLOAIRFOUJ.mat', ...
    'classify_train_cvglmnet_ofc_condition_PZBSZUISPC.mat', ...
    'classify_train_cvglmnet_striatum_condition_PIGWUYNZAE.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_FQSXSSULQU.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_GLXKTFUWEO.mat', ...
    'classify_train_cvglmnet_bg_condition_UXWIAHETEC.mat', ...
    'classify_train_cvglmnet_pallidum_condition_TWGTORCSRB.mat'}];
% held out trial 8
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_PTQGFFFNWB.mat', ...
    'classify_train_cvglmnet_ofc_condition_NXWDAUGTHH.mat', ...
    'classify_train_cvglmnet_striatum_condition_TPMUIYODJF.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_CISOQYHCNZ.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_DUJZKROUQV.mat', ...
    'classify_train_cvglmnet_bg_condition_WFRRJKATJF.mat', ...
    'classify_train_cvglmnet_pallidum_condition_NMLJFZBAKF.mat'}];
% held out trial 9
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_HWBJDXGXBV.mat', ...
    'classify_train_cvglmnet_ofc_condition_MOEBJDKRQT.mat', ...
    'classify_train_cvglmnet_striatum_condition_DJERIELAKN.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_GPGHWBDXXG.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_AVLZFFUZNB.mat', ...
    'classify_train_cvglmnet_bg_condition_SZYBWGIZMJ.mat', ...
    'classify_train_cvglmnet_pallidum_condition_IXICYVXZUX.mat'}];
% held out trial 10
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_JKXIIKNFWA.mat', ...
    'classify_train_cvglmnet_ofc_condition_VVXDAHQDIG.mat', ...
    'classify_train_cvglmnet_striatum_condition_GDJXJLRFLO.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_DVGGKTITTX.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_QOYLFSMNKQ.mat', ...
    'classify_train_cvglmnet_bg_condition_CWLWODIABQ.mat', ...
    'classify_train_cvglmnet_pallidum_condition_HCKFWTPQXB.mat'}];
% held out trial 11
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_JCPOOYUQRK.mat', ...
    'classify_train_cvglmnet_ofc_condition_CGWYLACFRU.mat', ...
    'classify_train_cvglmnet_striatum_condition_CGWYLACFRU.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_CGWYLACFRU.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_QPCOVZBOEM.mat', ...
    'classify_train_cvglmnet_bg_condition_JKXIIKNFWA.mat', ...
    'classify_train_cvglmnet_pallidum_condition_JKXIIKNFWA.mat'}];
% held out trial 12
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_IBMXXJXMKT.mat', ...
    'classify_train_cvglmnet_ofc_condition_IBMXXJXMKT.mat', ...
    'classify_train_cvglmnet_striatum_condition_IBMXXJXMKT.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_SKZANMDSVV.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_XQPHITPKBT.mat', ...
    'classify_train_cvglmnet_bg_condition_QFGJLJBXEH.mat', ...
    'classify_train_cvglmnet_pallidum_condition_BOLYLLFWAN.mat'}];
% held out trial 13
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_RFBTSDJRBJ.mat', ...
    'classify_train_cvglmnet_ofc_condition_BHLJIEBPGG.mat', ...
    'classify_train_cvglmnet_striatum_condition_RFBTSDJRBJ.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_RFBTSDJRBJ.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_RFBTSDJRBJ.mat', ...
    'classify_train_cvglmnet_bg_condition_KQSHRVCGZR.mat', ...
    'classify_train_cvglmnet_pallidum_condition_RFBTSDJRBJ.mat'}];
% held out trial 14
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_ZZQFTEPBWV.mat', ...
    'classify_train_cvglmnet_ofc_condition_ZZQFTEPBWV.mat', ...
    'classify_train_cvglmnet_striatum_condition_ONGRIVOUEL.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_ONGRIVOUEL.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_WGRBGCYDBP.mat', ...
    'classify_train_cvglmnet_bg_condition_WGRBGCYDBP.mat', ...
    'classify_train_cvglmnet_pallidum_condition_WGRBGCYDBP.mat'}];
% held out trial 15
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_YBLXIJCTJK.mat', ...
    'classify_train_cvglmnet_ofc_condition_YBLXIJCTJK.mat', ...
    'classify_train_cvglmnet_striatum_condition_UVYSNKPHGD.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_UVYSNKPHGD.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_ZZQFTEPBWV.mat', ...
    'classify_train_cvglmnet_bg_condition_ONGRIVOUEL.mat', ...
    'classify_train_cvglmnet_pallidum_condition_UZZSVENUHW.mat'}];
% held out trial 16
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_THUZJZNXQX.mat', ...
    'classify_train_cvglmnet_ofc_condition_ANWCHKGLOA.mat', ...
    'classify_train_cvglmnet_striatum_condition_ANWCHKGLOA.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_ULWBRUDAWX.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_ULWBRUDAWX.mat', ...
    'classify_train_cvglmnet_bg_condition_ULWBRUDAWX.mat', ...
    'classify_train_cvglmnet_pallidum_condition_AIXAOWATHY.mat'}];
% held out trial 17
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_TYVKCKYPKQ.mat', ...
    'classify_train_cvglmnet_ofc_condition_TYVKCKYPKQ.mat', ...
    'classify_train_cvglmnet_striatum_condition_TYVKCKYPKQ.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_TYVKCKYPKQ.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_KLXQRUCSYA.mat', ...
    'classify_train_cvglmnet_bg_condition_XBTTGATLYH.mat', ...
    'classify_train_cvglmnet_pallidum_condition_SSEEAGMFUW.mat'}];
% held out trial 18
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_DTJYVAPPHH.mat', ...
    'classify_train_cvglmnet_ofc_condition_ZRGZEKCRKX.mat', ...
    'classify_train_cvglmnet_striatum_condition_ZRGZEKCRKX.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_FSFRDGKWFR.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_FSFRDGKWFR.mat', ...
    'classify_train_cvglmnet_bg_condition_QAMTDFVUBY.mat', ...
    'classify_train_cvglmnet_pallidum_condition_QAMTDFVUBY.mat'}];
% held out trial 19
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_LUBGSMICBJ.mat', ...
    'classify_train_cvglmnet_ofc_condition_MSGHNNRZNE.mat', ...
    'classify_train_cvglmnet_striatum_condition_MSGHNNRZNE.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_MSGHNNRZNE.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_BJXGPGEKME.mat', ...
    'classify_train_cvglmnet_bg_condition_XSLBHVRBBC.mat', ...
    'classify_train_cvglmnet_pallidum_condition_BNGQOOQVXF.mat'}];
% held out trial 20
train_files = [train_files;
    {'classify_train_cvglmnet_hippocampus_condition_KQISMCFOKU.mat', ...
    'classify_train_cvglmnet_ofc_condition_SHCJEIGUIE.mat', ...
    'classify_train_cvglmnet_striatum_condition_JDAADPFEQL.mat', ...
    'classify_train_cvglmnet_vmpfc_condition_OSIJBKMSGN.mat', ...
    'classify_train_cvglmnet_rlpfc_condition_FOUWXRXSEY.mat', ...
    'classify_train_cvglmnet_bg_condition_FOUWXRXSEY.mat', ...
    'classify_train_cvglmnet_pallidum_condition_OBGQFGVQFB.mat'}];


held_out_runs = 1:9;
n_training_trials_per_run = 20;
n_test_trials_per_run = 4;

load_data;
[subjects, ~, ~] = contextGetSubjectsDirsAndRuns();
assert(isequal(subjects',unique(participant)));
isTest = ~isTrain;
newTrialId = trialId + isTest * n_training_trials_per_run; % number trials 1..24
sss = getGoodSubjects();

% summary stats of corr_coefs for each ROI, collapsed across subjects
sem = @(x) std(x) / sqrt(length(x));
means = [];
sems = [];

roi_corr_coefs = [];

for i = 1:length(rois)
    roi = rois{i};

    % for each subject, how well does the classifier's predition for each
    % condition (based on ROI activation) correlate with the Kalman filter's prediction
    % for the corresponding causal structure (based on test trial behavior)
    corr_coefs = []; % rows = subject, cols = causal structure / condition
    
    disp(roi);
        
    %
    % Get Y = P(condition), computed by the classifier on the held out trail
    % for all blocks / subjects
    % note that we average it across a bunch of held-out trials
    % TODO dedupe with classify_test
    %    
    outputs_from_all_heldout_trials = [];
    for t=1:length(held_out_trials)
        train_file = fullfile(folders{t}, train_files{t, i});
        load(train_file, 'CVfit');

        % IMPORTANT -- the parameters here must correspond to the ones that the classifier
        % was trained with
        mask = [roi, '.nii'];
        held_out_trial = held_out_trials(t);
        [inputs, targets, subjIds, runIds, trialIds] = classify_get_inputs_and_targets(held_out_trial, held_out_runs, sss, mask, 'condition', true, true);
        
        outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response');    
        accuracy = classify_get_accuracy(outputs, targets);
        % TODO sanity check with output frmo classify_scp_from_ncf
        fprintf('Success rate (held out trial = %d, lambda = %.4f) is %.2f%%\n', held_out_trial, CVfit.lambda_1se, accuracy);
        
        if isempty(outputs_from_all_heldout_trials)
            outputs_from_all_heldout_trials = outputs;
        else
            outputs_from_all_heldout_trials = outputs_from_all_heldout_trials + outputs;
        end
    end
    % average P(condition) across all held out runs 
    % TODO sanity check sems
    outputs = outputs_from_all_heldout_trials / length(held_out_trials);

    for subj = unique(subjIds)'
        fprintf('Kalman for subj %s\n', subjects{subj});
        %
        % Get X = P(structure | test trials) = P(test trials | structure) / sum P(test trials | structure)
        %       = P(test trial 1 | structure) * P(test trial 2 | structure) * P(test trial 3 | structure) * P(test trial 4 | structure) / sum ... 
        % as given by the Kalman filter, assuming a uniform P(structure)
        % for each subject, for each block
        %
        P_structure_all_runs = []; % = X
        P_condition_all_runs = outputs(subjIds == subj,:); % = Y
        actual_conditions_all_runs = []; % sanity check
        
        for run = runIds(subjIds == subj)'
            which_train = strcmp(participant, subjects{subj}) & roundId == run & isTrain;
            which_test =  strcmp(participant, subjects{subj}) & roundId == run & ~isTrain;

            which_models = [1 1 1 0];
            
            actual_condition = contextRole(which_train);
            actual_conditions_all_runs = [actual_conditions_all_runs; actual_condition(1)];
            
            % TODO dedupe with analyze.m
            
            % For a given run of a given subject, run the model on the same
            % sequence of stimuli and see what it does.
            %
            cues = cueId(which_train);
            N = length(cues); % # of trials
            assert(N == n_training_trials_per_run);
            D = 3; % # of stimuli
            prev_trials_surprise = zeros(N, D);
            prev_trials_surprise(sub2ind(size(prev_trials_surprise), 1:N, cues' + 1)) = 1;
            c = contextId(which_train) + 1;
            r = strcmp(sick(which_train), 'Yes');
            [choices, P_n, ww_n, P, ww, values] = train(prev_trials_surprise, c, r, prior_variance, inv_softmax_temp, which_models, false);
            
            % See what the model predicts for the test trials of that run
            %
            test_cues = cueId(which_test);
            test_N = length(test_cues); % # of trials
            assert(test_N == n_test_trials_per_run);
            D = 3; % # of stimuli
            test_x = zeros(test_N, D);
            test_x(sub2ind(size(test_x), 1:test_N, test_cues' + 1)) = 1;
            test_c = contextId(which_test) + 1;
            
            [test_choices, test_values, test_valuess, predict] = test(test_x, test_c, P_n, ww_n, inv_softmax_temp);

            % Compute X for the given run
            %
            P_response_given_structure = nan(n_test_trials_per_run, 3);
            for m = 1:3
                response_probs = predict(test_valuess(:, m));
                assert(length(response_probs) == n_test_trials_per_run);
                % P_blabla(i, j) = P(test trial i | structure j)
                P_response_given_structure(:,m) = binopdf(strcmp(response.keys(which_test), 'left')', 1, response_probs');
            end
            % P_blabla(i) = P(test trials | structure i)
            P_test_trials_given_structure = prod(P_response_given_structure);
            P_structure = P_test_trials_given_structure / sum(P_test_trials_given_structure);
            
            P_structure_all_runs = [P_structure_all_runs; P_structure];            
        end
        
        % Correlate X and Y for all blocks for that subject
        % for each structure / condition
        % i.e. see how well the classifier's prediction for 'modulatory' on
        % that run correlates with the probability that the causal structure is M1 
        % according to the Kalman filter and the subject's responding on
        % the test trials.
        %
        cs = [];
        for m1=1:3
            for m2=1:3
                c = corrcoef(P_condition_all_runs(:,m1), P_structure_all_runs(:,m2));
                cs = [cs c(2,1)];
            end
        end
        corr_coefs = [corr_coefs; cs];
        disp(cs);
    end
    
    means = [means; mean(corr_coefs)];
    sems = [sems; sem(corr_coefs)];
    
    roi_corr_coefs{i} = corr_coefs;
    
   % break; % TODO other ROIs too
end

figure;
barweb(means, sems);
xticklabels(rois);
%legend('irr-M1', 'mod-M2', 'add-M3');
legend('irr-M1', 'irr-M2', 'irr-M3', 'mod-M1', 'mod-M2', 'mod-M3', 'add-M1', 'add-M2', 'add-M3');
ylabel('correlation coefficient');
title('Corr(P_{classifier}(condition), P_{Kalman}(structure | test choices)), averaged across subjects');

save('classify_vs_kalman.mat');