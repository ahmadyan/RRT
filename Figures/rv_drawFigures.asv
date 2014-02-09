%Script to generate figures for the rv-journal

%TDO Figures
%Figure 1: Monte Carlo without any deviation
tdovcmin=-0.2;
tdovcmax=1.2;
tdoilmin=-0.02;
tdoilmax=0.1;
tdotmax=2e-6;
tdovcdim=3;
tdoildim=4;
tdotdim=5;

%%-------------------------------------------------------------------------
% Experiment 1: TDO oscillation
% TDO oscillates (R=0.2 in TDO)
%%-------------------------------------------------------------------------


% Exp 1, Case 1: Monte Carlo

tdoMCNominalData = rrt2mat('rv_tdo_sim_ok.rrt');

fig1=figure(1);
drawTree2D(tdoMCNominalData, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig1, '-painters', '-dpdf', '-r600' , 'tdo_sim_ok.pdf' );

fig2=figure(2);
drawTrace(tdoMCNominalData, tdovcdim, tdotdim, 'v_c(t)', tdotmax, tdovcmin, tdovcmax);
print ( fig2, '-painters', '-dpdf', '-r600' , 'tdo_sim_ok_vc.pdf' );

fig3=figure(3);
drawTrace(tdoMCNominalData, tdoildim, tdotdim, 'i_L(t)', tdotmax, tdoilmin, tdoilmax);
print ( fig3, '-painters', '-dpdf', '-r600' , 'tdo_sim_ok_iL.pdf' );


fig_mc1_3d=figure(15);
drawTree3D(tdoMCNominalData, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig_mc1_3d, '-painters', '-dpdf', '-r600' , 'tdo_mc3d_ok.pdf' );


% Exp 1, Case 2: Monte Carlo with deviations (0.05 on vc, 005 on iL)
tdoMCVariationData = rrt2mat('rv_tdo_sim_bad.rrt');

fig4=figure(4);
drawTree2D(tdoMCVariationData, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig4, '-painters', '-dpdf', '-r600' , 'tdo_sim_bad.pdf' );

fig5=figure(5);
drawTrace(tdoMCVariationData, tdovcdim, tdotdim, 'v_c(t)', tdotmax, tdovcmin, tdovcmax);
print ( fig5, '-painters', '-dpdf', '-r600' , 'tdo_sim_bad_vc.pdf' );


fig6=figure(6);
drawTrace(tdoMCVariationData, tdoildim, tdotdim, 'i_L(t)', tdotmax, tdoilmin, tdoilmax);
print ( fig6, '-painters', '-dpdf', '-r600' , 'tdo_sim_bad_iL.pdf' );

fig_mc2_3d=figure(16);
drawTree3D(tdoMCVariationData, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig_mc2_3d, '-painters', '-dpdf', '-r600' , 'tdo_mc3d_ok.pdf' );


% Exp 1, Case 3: RRT with deviation 1 (0.05 on vc, 005 on iL), 20,000 samples
tdoRRTVariationData1 = rrt2mat('rv_tdo_rrt_bad_0.05.rrt');

fig7=figure(7);
drawTree2D(tdoRRTVariationData1, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig7, '-painters', '-dpdf', '-r600' , 'tdo_rrt_bad.pdf' );

fig8=figure(8);
drawTrace(tdoRRTVariationData1, tdovcdim, tdotdim, 'v_c(t)', -1, tdovcmin, tdovcmax);
print ( fig8, '-painters', '-dpdf', '-r600' , 'tdo_rrt_bad_vc.pdf' );

fig9=figure(9);
drawTrace(tdoRRTVariationData1, tdoildim, tdotdim, 'i_L(t)', -1, tdoilmin, tdoilmax);
print ( fig9, '-painters', '-dpdf', '-r600' , 'tdo_rrt_bad_iL.pdf' );


fig10=figure(10);
drawTree3D(tdoRRTVariationData1, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig10, '-painters', '-dpdf', '-r600' , 'tdo_rrt3d_bad.pdf' );


% Exp 1, Case 4: RRT with deviation 2 (0.01 on vc, 001 on iL), 10,000 samples
tdoRRTVariationData2 = rrt2mat('rv_tdo_experiment1.rrt');

fig11=figure(11);
drawTree2D(tdoRRTVariationData2, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig11, '-painters', '-dpdf', '-r600' , 'tdo_rrt_ok.pdf' );

fig12=figure(12);
drawTrace(tdoRRTVariationData2, tdovcdim, tdotdim, 'v_c(t)', -1, tdovcmin, tdovcmax);
print ( fig12, '-painters', '-dpdf', '-r600' , 'tdo_rrt_ok_vc.pdf' );

fig13=figure(13);
drawTrace(tdoRRTVariationData2, tdoildim, tdotdim, 'i_L(t)', -1, tdoilmin, tdoilmax);
print ( fig13, '-painters', '-dpdf', '-r600' , 'tdo_rrt_ok_iL.pdf' );

fig14=figure(14);
drawTree3D(tdoRRTVariationData2, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig14, '-painters', '-dpdf', '-r600' , 'tdo_rrt3d_ok.pdf' );




%%-------------------------------------------------------------------------
% Experiment 2: TDO No oscillation (R=0.5 in TDO)
% Using low variation, circuit doesn't oscillate in MC, but it does in RRT
% Using high variation, circuit oscillates in both cases
%%-------------------------------------------------------------------------

% Exp 2 Case 1: MC, R=0.5, Var=vc->0.005, iL->0.0005
% input file name=tdo_sim_no_osc.rrt
data2_1 = rrt2mat('tdo_sim_no_osc.rrt');

fig2_1_1=figure(20);
drawTree2D(data2_1, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig2_1_1, '-painters', '-dpdf', '-r600' , 'tdo_exp2_1_1.pdf' );

fig2_1_2=figure(21);
drawTrace(data2_1, tdovcdim, tdotdim, 'v_c(t)', -1, tdovcmin, tdovcmax);
print ( fig2_1_2, '-painters', '-dpdf', '-r600' , 'tdo_exp2_1_2.pdf' );

fig2_1_3=figure(22);
drawTrace(data2_1, tdoildim, tdotdim, 'i_L(t)', -1, tdoilmin, tdoilmax);
print ( fig2_1_3, '-painters', '-dpdf', '-r600' , 'tdo_exp2_1_3.pdf' );

fig2_1_4=figure(23);
drawTree3D(data2_1, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig2_1_4, '-painters', '-dpdf', '-r600' , 'tdo_exp2_1_4.pdf' );


% Exp 2 Case 2: MC, R=0.5, Var=vc->0.05, iL->0.005
%input file name=tdo_sim_no_osc_high_var.rrt
data2_2 = rrt2mat('tdo_sim_no_osc_high_var.rrt');

fig2_2_1=figure(30);
drawTree2D(data2_2, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig2_2_1, '-painters', '-dpdf', '-r600' , 'tdo_exp2_2_1.pdf' );

fig2_2_2=figure(31);
drawTrace(data2_2, tdovcdim, tdotdim, 'v_c(t)', -1, tdovcmin, tdovcmax);
print ( fig2_2_2, '-painters', '-dpdf', '-r600' , 'tdo_exp2_2_2.pdf' );

fig2_2_3=figure(32);
drawTrace(data2_2, tdoildim, tdotdim, 'i_L(t)', -1, tdoilmin, tdoilmax);
print ( fig2_2_3, '-painters', '-dpdf', '-r600' , 'tdo_exp2_2_3.pdf' );

fig2_2_4=figure(33);
drawTree3D(data2_2, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig2_2_4, '-painters', '-dpdf', '-r600' , 'tdo_exp2_2_4.pdf' );



% Exp 2 Case 3: RRT, R=0.5, Var=vc->0.005, iL->0.0005
%input file name=tdo_rrt_no_osc.rrt
data2_3 = rrt2mat('tdo_rrt_no_osc.rrt');

fig2_3_1=figure(40);
drawTree2D(data2_3, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig2_3_1, '-painters', '-dpdf', '-r600' , 'tdo_exp2_3_1.pdf' );

fig2_3_2=figure(41);
drawTrace(data2_3, tdovcdim, tdotdim, 'v_c(t)', -1, tdovcmin, tdovcmax);
print ( fig2_3_2, '-painters', '-dpdf', '-r600' , 'tdo_exp2_3_2.pdf' );

fig2_3_3=figure(42);
drawTrace(data2_3, tdoildim, tdotdim, 'i_L(t)', -1, tdoilmin, tdoilmax);
print ( fig2_3_3, '-painters', '-dpdf', '-r600' , 'tdo_exp2_3_3.pdf' );

fig2_3_4=figure(43);
drawTree3D(data2_3, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig2_3_4, '-painters', '-dpdf', '-r600' , 'tdo_exp2_3_4.pdf' );


% Exp 2 Case 4: RRT, R=0.5, Var=vc->0.05, iL->0.005
%input file name=tdo_rrt_no_osc_high_var.rrt
data2_4 = rrt2mat('tdo_rrt_no_osc_high_var.rrt');

fig2_4_1=figure(50);
drawTree2D(data2_4, tdovcdim, tdoildim, 'v_c', 'i_L', tdovcmin, tdovcmax, tdoilmin, tdoilmax);
print ( fig2_4_1, '-painters', '-dpdf', '-r600' , 'tdo_exp2_4_1.pdf' );

fig2_4_2=figure(51);
drawTrace(data2_4, tdovcdim, tdotdim, 'v_c(t)', -1, tdovcmin, tdovcmax);
print ( fig2_4_2, '-painters', '-dpdf', '-r600' , 'tdo_exp2_4_2.pdf' );

fig2_4_3=figure(52);
drawTrace(data2_4, tdoildim, tdotdim, 'i_L(t)', -1, tdoilmin, tdoilmax);
print ( fig2_4_3, '-painters', '-dpdf', '-r600' , 'tdo_exp2_4_3.pdf' );

fig2_4_4=figure(53);
drawTree3D(data2_4, tdovcdim, tdoildim, tdotdim, 'v_c(t)', 'i_L(t)', tdovcmin, tdovcmax, tdoilmin, tdoilmax, -1);
view(-45, 45)
drawnow
print ( fig2_4_4, '-painters', '-dpdf', '-r600' , 'tdo_exp2_4_4.pdf' );


%%-------------------------------------------------------------------------
% Experiment 3: Ring modulator results
%%-------------------------------------------------------------------------




%%-------------------------------------------------------------------------
% Experiment 4: PLL
%%-------------------------------------------------------------------------
pll_e = 3;
pll_eb = 4;
pll_in = 5;
pll_inb = 6;
pll_mout = 7;
pll_moutb = 8;
pll_osc = 9;
pll_oscb = 10;
pll_out = 11;
pll_outb = 12;
pll_xvco_c = 13;
pll_xvco_s = 14;
pll_xvco_s_clip	= 15;
pll_xpd_clip1 = 16;
pll_xpd_clip2 = 17;
pll_xpd_n1 = 18;
pll_time = 19;
% Case i: pll MC, low var
fig_4_1_1=figure(411);
data41=rrt2mat('pll_sim_ok_10000.rrt');
drawTrace(data41, pll_e, pll_time, 'pll_e', -1, -1, -1);





%%-------------------------------------------------------------------------
% Experiment 5: INV
%%-------------------------------------------------------------------------


inv_vdd= 3;
inv_vin= 4;
inv_vout=5;
inv_gnd=6;
inv_time=31; 
%Case 1: MC
data51=rrt2mat('inverter_mc_1000.rrt');
Fig_5_1_1 = figure(511);
drawEye(data51, inv_vout, inv_time, 100e-12, 'Voltage(v)', -0.2, 1.2)

%Case 2: RRT, with variation on input and supply
data52=rrt2mat('inv_rrt_2000.rrt');
Fig_5_2_1 = figure(531);
drawTrace(data52, inv_vout, inv_time, 'v_out', -1, -0.2, 1.2);

Fig_5_2_2 = figure(532);
drawTree2D( data52 , inv_vin, inv_vout, 'vin', 'vout', -0.2, 1.2, -0.2, 1.2)

%Case 3: RRT, with variation on input and supply
%data53=rrt2mat('inv_rrt_10000_fp.rrt');
data53=rrt2mat('inverter10.rrt');
Fig_5_3_1 = figure(531);
drawTrace(data53, inv_vout, inv_time, 'v_{out}', -1, -0.2, 1.2);

Fig_5_3_2 = figure(532);
drawTree2D( data53 , inv_vin, inv_vout, 'vin', 'vout', -0.2, 1.2, -0.2, 1.2)

Fig_5_3_3 = figure(533);
drawEye(data53, inv_vout, inv_time, 100e-12, 'Voltage(v)', -0.2, 1.2);



for i=3:size(data51, 1),
    x=min(data51(i, :));
    y=max(data51(i, :));

    x = x - 0.2*x; 
    y = y + 0.2*y;
    disp(sprintf('edu.uiuc.csl.system.var.ic[%d]=0', i-3));
    disp(sprintf('edu.uiuc.csl.system.var.min[%d]=%f', i-3, x));
    disp(sprintf('edu.uiuc.csl.system.var.max[%d]=%f', i-3, y));
    disp(sprintf('edu.uiuc.csl.system.var.name[%d]=v%d', i-3, i-3));

end