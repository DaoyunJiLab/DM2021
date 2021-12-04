function hmain = DataManager_DataManager_Callback
%%This is a big standalone function for analyzing spike and EEG data

[~, ~, ~, ~, ~, MAname] = CurrentVersion;

hf = figure('Name', MAname, 'NumberTitle', 'off', 'NextPlot', 'add',...
    'MenuBar', 'figure', 'Units', 'normalized', 'Position', [0.05 0.2 0.9 0.7]);
uimenu('Parent', hf, 'Label', '          ');

hfile = uimenu('Parent', hf, 'Label', '  DataBase  ');
    uimenu(hfile, 'Label', 'Generate New SpikeDatabase', 'Callback', 'DataManager_GenerateSpikeDatabase_Callback');
    uimenu(hfile, 'Label', 'Generate New BehavDatabase', 'Callback', 'DataManager_GenerateBehavDatabase_Callback');
    uimenu(hfile, 'Label', 'Generate New EEGDatabase', 'Callback', 'DataManager_GenerateEEGDatabase_Callback');

    uimenu(hfile, 'Label', 'Combine Database', 'Callback', 'DataManager_CombineDatabase_Callback', 'Separator', 'on');

    uimenu(hfile, 'Label', 'Open Database', 'Callback', 'DataManager_OpenDatabase_Callback', 'Separator', 'on');
    uimenu(hfile, 'Label', 'SaveAs', 'Callback', 'DataManager_SaveAs_Callback', 'Separator', 'on');
    uimenu(hfile, 'Label', 'SaveSelect', 'Callback', 'DataManager_SaveSelect_Callback');
    uimenu(hfile, 'Label', 'SaveStripped', 'Callback', 'DataManager_SaveStripped_Callback');
    uimenu(hfile, 'Label', 'Save', 'Callback', 'DataManager_Save_Callback');
    uimenu(hfile, 'Label', 'Close', 'Callback', 'DataManager_Close_Callback');
    
hcross = uimenu('Parent', hf, 'Label', '  Cross Analysis  ');  
   %%%%%%%%%%%sequence analysis not worked out yet!!!
   hseq = uimenu('Parent', hcross, 'Label', 'Sequence/Replay');
        hgtemp = uimenu(hseq, 'Label', 'Generate Template');
               uimenu(hgtemp, 'Label', '2D Active cells', 'Callback', 'DM_Seq_GenerateTemplate_2DActiveCells_Callback');
               %uimenu(hgtemp, 'Label', '2D Place Fields', 'Callback', 'DM_Seq_GenerateTemplate_2DPlaceFields_Callback'); %%this option does not make sense
               uimenu(hgtemp, 'Label', '1D Active cells', 'Callback', 'DM_Seq_GenerateTemplate_1DActiveCells_Callback', 'Separator', 'on');
               uimenu(hgtemp, 'Label', '1D Place Fields', 'Callback', 'DM_Seq_GenerateTemplate_PlaceFields_Callback');
               uimenu(hgtemp, 'Label', 'Cross Correlation', 'Callback', 'DM_Seq_GenerateTemplate_CrossCrr_Callback', 'Separator', 'on');
               
        uimenu(hseq, 'Label', 'Display Templates', 'Callback', 'DE_Seq_DisplayTemplate_Callback', 'Separator', 'on');
        uimenu(hseq, 'Label', 'Modify Templates', 'Callback', 'DE_Seq_ModifyTemplate_Callback');
        uimenu(hseq, 'Label', 'Break Templates', 'Callback', 'DE_Seq_BreakTemplate_Callback');
        uimenu(hseq, 'Label', 'Combine Templates', 'Callback', 'DE_Seq_CombineTemplate_Callback');
        uimenu(hseq, 'Label', 'Equalize Templates', 'Callback', 'DE_Seq_EqualizeTemplate_Callback');
        uimenu(hseq, 'Label', 'Sequence Matching', 'Callback', 'DM_Seq_SequenceMatching_events_Callback', 'Separator', 'on');
        %uimenu(hseq, 'Label', 'Event Sequence Matching', 'Callback', 'DM_Seq_SequenceMatching_events_Callback');
        uimenu(hseq, 'Label', 'Break SequenceDB', 'Callback', 'DM_Seq_BreakSequenceDB_Callback');
        %uimenu(hseq, 'Label', 'Combine SequenceDB', 'Callback', 'DM_Seq_CombineSequenceDB_Callback');
        %uimenu(hseq, 'Label', 'Template Matching', 'Callback', 'DM_Seq_TemplateMatching_Callback'); %, 'Separator', 'on');
        uimenu(hseq, 'Label', '1D Bayesian Decoding - template itemized', 'Callback', 'DM_Seq_BayesianDecoding_events_Callback', 'Separator', 'on');
        uimenu(hseq, 'Label', '1D Bayesian Decoding - event itemized', 'Callback', 'DM_Seq_BayesianDecoding_eventsitemized_Callback');
        uimenu(hseq, 'Label', '2D Bayesian Decoding', 'Callback', 'DM_Seq_2DBayesianDecoding_events_Callback');
    hpv = uimenu(hcross, 'Label', 'Population vector correlation');
        uimenu(hpv, 'Label', 'Compute PV corr', 'Callback', 'DataManager_computePVcrr_Callback');
        uimenu(hpv, 'Label', 'Display example PVcrr', 'Callback', 'DataManager_displayPVcrr_Callback');
        uimenu(hpv, 'Label', 'Plot PVcrr distributions', 'Callback', 'DataManager_plotPVcrrDistribution_Callback');
    hcrr = uimenu(hcross, 'Label', 'Spike train cross-correlation', 'Callback', 'DataManager_GenerateCrrDatabase_Callback');
    
    %hsync = uimenu(hcross, 'Label', 'Synchronization', 'Callback', 'DataManager_GenerateSyncDatabase_Callback');

hrecompute = uimenu('Parent', hf, 'Label', 'Recompute/Modify');    
    uimenu(hrecompute, 'Label', 'ReComputeDataBase_NoStops', 'Callback', 'DataManager_ReComputeDatabase_NoStops');
    uimenu(hrecompute, 'Label', 'ReComputeDataBase_GoodLapsOnly', 'Callback', 'DataManager_ReComputeDatabase_GoodLaps');
    uimenu(hrecompute, 'Label', 'Behavior-modified place field dynamics', 'Callback', 'DataManager_BehaviorModifyPFdynamics_Callback', 'Separator', 'on');
    uimenu(hrecompute, 'Label', 'LapConsistency_PlaceField', 'Callback', 'DataManager_FindLapConsistency_Rate', 'Separator', 'on');
    uimenu(hrecompute, 'Label', 'LapConsistency_CrossCrr', 'Callback', 'DataManager_FindLapConsistency_CrossCrr');
    uimenu(hrecompute, 'Label', 'LapConsistency_SpatialCrossCrr', 'Callback', 'DataManager_FindLapConsistency_SpatialCrossCrr');
    uimenu(hrecompute, 'Label', 'SplitSession_FieldDynamics', 'Callback', 'DataManager_SplitSessLapConsistency_PlaceFields', 'Separator', 'on');
    uimenu(hrecompute, 'Label', 'SplitSession_PVcorrelation', 'Callback', 'DataManager_SplitSess_PVcrr');
hppro = uimenu('Parent', hf, 'Label', '  Projects  ');
    hV1CA1 = uimenu('Parent', hppro, 'Label', 'V1-CA1 pre-SWS interaction');
       uimenu(hV1CA1, 'Label', 'pre-SWS events', 'Callback', 'DataManager_V1CA1_FindpreSWSevents_Callback');
       uimenu(hV1CA1, 'Label', 'pre-HVS events', 'Callback', 'DataManager_V1CA1_FindpreHVSevents_Callback');
       uimenu(hV1CA1, 'Label', 'Ripple Incidence', 'Callback', 'DataManager_V1CA1_FindRippleIncidence_Callback');
       uimenu(hV1CA1, 'Label', 'DownState Incidence', 'Callback', 'DataManager_V1CA1_FindDownStateIncidence_Callback');
       uimenu(hV1CA1, 'Label', 'MUAlevel', 'Callback', 'DataManager_V1CA1_FindMUAlevel_Callback');
       uimenu(hV1CA1, 'Label', 'Cl0Crr', 'Callback', 'DataManager_V1CA1_FindCluster0Crr_Callback');
       uimenu(hV1CA1, 'Label', 'Triggered Cl0 Averages', 'Callback', 'DataManager_V1CA1_FindCluster0TriggeredAverages_Callback');
       uimenu(hV1CA1, 'Label', 'SWS SpindleRipple Correlation', 'Callback', 'DataManager_V1CA1_FindRippleSpindleCrr_SWS_Callback');
       uimenu(hV1CA1, 'Label', 'SWS UPstate PV Correlation', 'Callback', 'DataManager_V1CA1_SWS_UPstate_PVCrr_Callback');
    hTauM = uimenu('Parent', hppro, 'Label', 'Tau Mice');
       uimenu(hTauM, 'Label', 'ReComputeDataBase_GoodLapsOnly', 'Callback', 'DataManager_ReComputeDatabase_GoodLaps');
       uimenu(hTauM, 'Label', 'ReComputeDataBase_NoStops', 'Callback', 'DataManager_ReComputeDatabase_NoStops');
       uimenu(hTauM, 'Label', 'FreeSequencing', 'Callback', 'DM_Seq_FreeSequencing', 'Separator', 'on');
       uimenu(hTauM, 'Label', 'FreeSequencingWithDownSampledSpikes', 'Callback', 'DM_FreeSequencingWithDownSampledSpikes');
       uimenu(hTauM, 'Label', 'PlotSequenceResults', 'Callback', 'DM_Seq_PlotSequenceResults', 'Separator', 'on');
       uimenu(hTauM, 'Label', 'PlotSequencePositionStability', 'Callback', 'DM_Seq_PlotSequencePositionStability');
       uimenu(hTauM, 'Label', 'PlotSequenceVelocity', 'Callback', 'DM_Seq_PlotSequenceVelocity');
       uimenu(hTauM, 'Label', 'PlotSequenceInterval', 'Callback', 'DM_Seq_PlotSequenceInterval');
       uimenu(hTauM, 'Label', 'PlotSequenceLength', 'Callback', 'DM_Seq_PlotSequenceLength');
       uimenu(hTauM, 'Label', 'PlotSequencePreferredDirection', 'Callback', 'DM_Seq_PlotSequenceLocation');
       uimenu(hTauM, 'Label', 'CheckHighCrrCell_Rate', 'Callback', 'DataManager_CheckHighCrrCell_Rate', 'Separator', 'on');
       uimenu(hTauM, 'Label', 'CheckHighCrrConnectivity', 'Callback', 'DataManager_CheckHighCrrConnectivity');
       uimenu(hTauM, 'Label', 'CompareSessionCrr', 'Callback', 'DataManager_CompareSessionCrr');
       uimenu(hTauM, 'Label', 'CompareEventCrr', 'Callback', 'DataManager_CompareEventCrr');
       uimenu(hTauM, 'Label', 'CompareSessionCrr_tmp', 'Callback', 'DataManager_CompareSessionCrr_tmp');
       uimenu(hTauM, 'Label', 'CompareEventCrr_tmp', 'Callback', 'DataManager_CompareEventCrr_tmp');
       uimenu(hTauM, 'Label', 'CompareSessionRate', 'Callback', 'DataManager_CompareSessionRate');
       uimenu(hTauM, 'Label', 'CompareEventRate', 'Callback', 'DataManager_CompareEventRate');
       uimenu(hTauM, 'Label', 'SpikeSpectrum', 'Callback', 'DataManager_FindSpikeSpectrum', 'Separator', 'on');
       uimenu(hTauM, 'Label', 'EEGSpectrumWithinSequences', 'Callback', 'DataManager_FindEEGSpectrumWithinSeqs');
       
       %uimenu(hTauM, 'Label', 'Show1DPF_GoodLapsOnly', 'Callback', 'DataManager_Show1DPlaceFields_GoodLaps');
       %uimenu(hTauM, 'Label', 'LapConsistency_GoodLapsOnly', 'Callback', 'DataManager_FindLapConsistency_GoodLaps');
       %uimenu(hTauM, 'Label', 'PF_NoStops', 'Callback', 'DataManager_FindSpatialInfoFieldProp_NoStops');
       %uimenu(hTauM, 'Label', 'Show1DPF_NoStops', 'Callback', 'DataManager_Show1DPlaceFields_NoStops');
       %uimenu(hTauM, 'Label', 'LapConsistency_NoStops', 'Callback', 'DataManager_FindLapConsistency_NoStops');
    hTauR = uimenu('Parent', hppro, 'Label', 'TauRippleDeficits');
       uimenu(hTauR, 'Label', 'Ripple_SharpWave Calibration', 'Callback', 'DataManager_EEG_RippleSharpWaveCalibration');
       uimenu(hTauR, 'Label', 'Ripple_TrigAvg_quantification', 'Callback', 'DataManager_EEG_RippleTrigAvgQuantification');
       uimenu(hTauR, 'Label', 'Crr within another interval', 'Callback', 'DataManager_EEG_AddAnotherCrrInterval');
%     hVFF = uimenu('Parent', hppro, 'Label', 'Visual Firing Field');
%     hADM = uimenu('Parent', hppro, 'Label', 'AD Mouse Models');
%     hOL = uimenu('Parent', hppro, 'Label', 'Observational Learning');
%     hSI = uimenu('Parent', hppro, 'Label', 'Social Interaction');
%           uimenu(hSI, 'Label', 'Detect HighTheta Events', 'Callback', 'DataManager_SI_DetectThetaEvents_Callback');
    
    hRTT = uimenu('Parent', hppro, 'Label', 'Rett Mice');
       uimenu(hRTT, 'Label', 'Identify Pairs with Overlapping Place Fields', 'Callback', 'DataManager_RettMice_OverlappingPlaceCellPairs');
       uimenu(hRTT, 'Label', 'Quantify Event Participation', 'Callback', 'DataManager_RettMice_EventParticipation');
       uimenu(hRTT, 'Label', 'Revise Place Field Properties', 'Callback', 'DataManager_RettMice_RevisePlaceFieldProperties');   
       uimenu(hRTT, 'Label', 'DSandLapStatistics', 'Callback', 'DataManager_RettMice_DSandLapStatistics');  
    
    hV1CA1Track = uimenu('Parent', hppro, 'Label', 'V1-CA1 track interaction');
       uimenu(hV1CA1Track, 'Label', 'Add shuffle-normalized spatial information', 'Callback', 'DataManager_V1CA1_AddNormalizedSpatioInfo');
       uimenu(hV1CA1Track, 'Label', 'Day by day field/firing dynamics', 'Callback', 'DataManager_V1CA1_DayDynamics', 'Separator', 'on');
       uimenu(hV1CA1Track, 'Label', 'Lap by lap field dynamics', 'Callback', 'DataManager_V1CA1_LapDynamics');
       uimenu(hV1CA1Track, 'Label', 'Lap field dynamics with speed correction', 'Callback', 'DataManager_V1CA1_LapDynamics_speedcorrection');
       uimenu(hV1CA1Track, 'Label', 'Day by day crr dynamics', 'Callback', 'DataManager_V1CA1_CrrDayDynamics');
       uimenu(hV1CA1Track, 'Label', 'Lap by lap crr dynamics', 'Callback', 'DataManager_V1CA1_CrrLapDynamics');
       uimenu(hV1CA1Track, 'Label', 'Lap crr dynamics with speed correction', 'Callback', 'DataManager_V1CA1_CrrLapDynamics_speedcorrection');
       uimenu(hV1CA1Track, 'Label', 'Pairwise lap noise correlation', 'Callback', 'DataManager_V1CA1_CrrLapNoiseDynamics');
       uimenu(hV1CA1Track, 'Label', 'Lap noise correlation with speed correction', 'Callback', 'DataManager_V1CA1_CrrLapNoiseDynamics_speedcorrection');
       uimenu(hV1CA1Track, 'Label', 'Day by day behavior dynamics', 'Callback', 'DataManager_V1CA1_DaySpeedDynamics', 'Separator', 'on');
       uimenu(hV1CA1Track, 'Label', 'Lap by lap speed dynamics', 'Callback', 'DataManager_V1CA1_LapSpeedDynamics');
       uimenu(hV1CA1Track, 'Label', 'Lap by lap field_speed correlation', 'Callback', 'DataManager_V1CA1_LapSpeedFieldNoiseCorrelation');
       uimenu(hV1CA1Track, 'Label', 'Lap by lap crr_speed correlation', 'Callback', 'DataManager_V1CA1_LapSpeedCrrNoiseCorrelation');
       uimenu(hV1CA1Track, 'Label', 'Crr comparison between halves of place fields', 'Callback', 'DataManager_V1CA1_CrrHalfPlaceFields', 'Separator', 'on');
       uimenu(hV1CA1Track, 'Label', 'Compare ripple correlation', 'Callback', 'DataManager_V1CA1_CompareRippleCorrelation', 'Separator', 'on');
       uimenu(hV1CA1Track, 'Label', 'Light_dark correlation', 'Callback', 'DataManager_V1CA1_LightDarkCorrelation', 'Separator', 'on');
       uimenu(hV1CA1Track, 'Label', 'Field Distribution', 'Callback', 'DataManager_V1CA1_FieldDistribution');
       uimenu(hV1CA1Track, 'Label', 'RetroProspectCode_CrrTime_OverlapPairs', 'Callback', 'DataManager_V1CA1_RetroProspectCode_CrrTime_OverlapPairs', 'Separator', 'on');

   hmr = uimenu('Parent', hppro, 'Label', 'MouseRatComparison');
       uimenu(hmr, 'Label', 'SpatialInfo dynamics on all active trajectories', 'Callback', 'DataManager_SIDynamOnAllActiveTraj_Callback');
       uimenu(hmr, 'Label', 'Decoding accuracy dynamics on template trajectories', 'Callback', 'DataManager_DecodeDynamOnTmpTraj_Callback');
       uimenu(hmr, 'Label', 'Within periods replay statistics and significance', 'Callback', 'DataManager_WithinPeriodReplayStatsSig_Callback');
       uimenu(hmr, 'Label', 'Chop events', 'Callback', 'DataManager_ChopEvents_Callback');
   
   hobs = uimenu('Parent', hppro, 'Label', 'Observation');
       uimenu(hobs, 'Label', 'Add short run events', 'Callback', 'DataManager_Obs_AddShortRunEvents_Callback');
          
   hlsd = uimenu('Parent', hppro, 'Label', 'LSD');
       uimenu(hlsd, 'Label', 'Place field responses to (HT) event times', 'Callback', 'DataManager_PFchangesInEvents_Callback');
       uimenu(hlsd, 'Label', '-- Plot results', 'Callback', 'DataManager_plotPFPVchangesInEvents_Callback');
       uimenu(hlsd, 'Label', '-- Evalute shuffle significance', 'Callback', 'DataManager_evaluateP_PFPVchangesInEvents_Callback');
       uimenu(hlsd, 'Label', 'PV responses to (HT) event times', 'Callback', 'DataManager_PVchangesInEvents_Callback');
       uimenu(hlsd, 'Label', '-- Plot results', 'Callback', 'DataManager_plotPFPVchangesInEvents_Callback');
       uimenu(hlsd, 'Label', 'Crr_DownSample pairs with equal rates', 'Callback', 'DataManager_Crr_DSPairsEqualRates_Callback');
       uimenu(hlsd, 'Label', 'Speed_ThetaPeakFreqPower_running', 'Callback', 'DataManager_EEG_speedThetaPeakPowerRun_Callback');
       
   hexp = uimenu('Parent', hppro, 'Label', 'Experience Tracking');
    
   hfear = uimenu('Parent', hppro, 'Label', 'Contextual fear memory');
       uimenu(hfear, 'Label', 'Identify position patches', 'Callback', 'DataManager_IdentifyPosPatches_Callback');
   
if nargout == 1
   hmain = hf;
end