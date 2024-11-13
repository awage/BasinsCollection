using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes


function circadian_cell_cycle!(dX, X, p, t)
Mp, Mc, Mbmal, Pc, Cc, Pcp, Ccp, PCc, PCn, PCcp, PCnp, Bc, Bcp, Bn, Bnp, In, Mr, Rc, Rn, AP1, pRB, pRBc1, pRBp, pRBc2, pRBpp, E2F, E2Fp, Cd, Mdi, Md, Mdp27, Mce, Ce, Mei, Me, Skp2, Mep27, Pei, Pe, Ca, Mai, Map27, p27, p27p, Cdh1i, Cdh1a, Pai, Pa, Cb, Mbi, Mb, Mbp27, Cdc20i, Cdc20a, Pbi, Pb, Mw, Wee1, Wee1p = X

#parameters
k1_clock=0.8;k2_clock=0.4;k3_clock=0.8;k4_clock=0.4;k5=0.8;k6=0.4;k7=1;k8=0.2;k9=0.8;k10=0.4;
K_AP=0.6;K_AC=0.6;K_AR=0.6;K_IB=1;
k_dmb=0.02;k_dmc=0.02;k_dmp=0.02;k_dmr=0.02;k_dn=0.02;k_dnc=0.02;
K_d=0.3;K_dp=0.1;K_p=1.006;K_mB=0.4;K_mC=0.4;K_mP=0.3;K_mR=0.4;k_sB=0.32;
k_sC=3.2;k_sP=1.2;k_sR=1.7;m=2;h=2;n=2;
V_1B=1.4;V_1C=1.2;V_1P=9.6;V_1PC=2.4;V_2B=0.2;V_2C=0.2;V_2P=0.6;V_2PC=0.2;
V_3B=1.4;V_3PC=2.4;V_4B=0.4;V_4PC=0.2;V_phos=0.4;v_dBC=3;v_dBN=3;v_dCC=1.4;
v_dIN=1.6;v_dPC=3.4;v_dPCC=1.4;v_dPCN=1.4;v_dRC=4.4;v_dRN=0.8;v_mB=1.3;
v_mC=2.0;v_mP=2.2;v_mR=1.6;v_sB=1.8;                                        
v_sC=2.2;v_sP=2.4;v_sR=1.6

#Cell cycle
Chk1=0;u=0.04
GF=1;K_agf=0.1;k_dap1=0.15;eps=21.4;v_sap1=1
k_de2f=0.002;k_de2fp=1.1;k_dprb=0.01;k_dprbp=0.06;k_dprbpp=0.04
k_pc1=0.05;k_pc2=0.5;k_pc3=0.025;k_pc4=0.5;K1=0.1;K2=0.1;K3=0.1
K4=0.1;V1=2.2;V2=2;V3=1;V4=2;K_1e2f=5;K_2e2f=5;V_1e2f=4
V_2e2f=0.75;v_se2f=0.15;v_sprb=0.8
Cdk4_tot=1.5;K_i7=0.1;K_i8=2;k_cd1=0.4;k_cd2=0.005;k_decom1=0.1
k_com1=0.175;k_c1=0.15;k_c2=0.05;k_ddd=0.005;K_dd=0.1;K_1d=0.1;K_2d=0.1
V_dd=5;V_m1d=1;V_m2d=0.2
a_e=0.25;Cdk2_tot=2;i_b1=0.5;K_i9=0.1;K_i10=2;k_ce=0.29;k_c3=0.2
k_c4=0.1;k_decom2=0.1;k_com2=0.2;k_dde=0.005;k_ddskp2=0.005;k_dpe=0.075
k_dpei=0.15;K_de=0.1;K_dceskp2=2;K_dskp2=0.5;K_cdh1=0.4;K_1e=0.1
K_2e=0.1;K_5e=0.1;K_6e=0.1;V_de=3;V_dskp2=1.1;V_m1e=2;V_m2e=1.4;V_m5e=5
V_6e=0.8;v_spei=0.13;v_sskp2=0.15;x_e1=1;x_e2=1
a_a=0.2;i_b2=0.5;K_i11=0.1;K_i12=2;K_i13=0.1;K_i14=2;k_ca=0.0375
k_decom3=0.1;k_com3=0.2;k_c5=0.15;k_c6=0.125;k_dda=0.005;k_ddp27=0.06
k_ddp27p=0.01;k_dcdh1a=0.1;k_dcdh1i=0.2;k_dpa=0.075;k_dpai=0.15;K_da=1.1
K_dp27p=0.1;K_dp27skp2=0.1;K_acdc20=2;K_1a=0.1;K_2a=0.1;K_1cdh1=0.01
K_2cdh1=0.01;K_5a=0.1;K_6a=0.1;K_1p27=0.5;K_2p27=0.5;V_dp27p=5;V_da=2.5
V_m1a=2;V_m2a=1.85;V_m5a=4;V_6a=1;v_scdh1a=0.11;v_spai=0.105
v_s1p27=0.8;v_s2p27=0.1;V_1cdh1=1.25;V_2cdh1=8;V_1p27=100;V_2p27=0.1
x_a1=1;x_a2=1
a_b=0.11;Cdk1_tot=0.5;i_b=0.75;i_b3=0.5;k_c7=0.12;k_c8=0.2
k_decom4=0.1;k_com4=0.25;k_dcdc20a=0.05;k_dcdc20i=0.14;k_ddb=0.005
k_dpb=0.1;k_dpbi=0.2;k_dwee1=0.1;k_dwee1p=0.2;K_db=0.005;K_dbcdc20=0.2;K_dbcdh1=0.1
k_sw=5;K_1b=0.1;K_2b=0.1;K_3b=0.1;K_4b=0.1;K_5b=0.1;K_6b=0.1;K_7b=0.1
K_8b=0.1;v_cb=0.055;V_db=0.06;V_m1b=3.9;V_m2b=2.1;v_scdc20i=0.1;V_m3b=8;V_m4b=0.7
V_m5b=5;V_6b=1;V_m7b=1.2;V_m8b=1;v_spbi=0.12;x_b1=1;x_b2=1

#coupling 
v_sw=1;v_swee1=0.0117;nmw=4;K_aw=2;V_dmw=0.5;K_dmw=0.5
v_sce = 1.0;
K_ice=1;V_dmce=0.5;K_dmce=0.5;nce=4;k_ce2=5

g(x,y) = x/(x+y)
f(x,y,n) = x^n/(x^n+y^n)

dMp = v_sP*f(Bn,K_AP,n)-v_mP*g(Mp,K_mP)-k_dmp*Mp   
dMc = v_sC*f(Bn,K_AC,n)-v_mC*g(Mc,K_mC)-k_dmc*Mc      
dMbmal = v_sB*f(K_IB,Rn,m)-v_mB*g(Mbmal,K_mB)-k_dmb*Mbmal     
dPc = k_sP*Mp-V_1P*g(Pc,K_p)+V_2P*g(Pcp,K_dp)+k4_clock*PCc-k3_clock*Pc*Cc-k_dn*Pc  
dCc = k_sC*Mc-V_1C*g(Cc,K_p)+V_2C*g(Ccp,K_dp)+k4_clock*PCc-k3_clock*Pc*Cc-k_dnc*Cc  
dPcp = V_1P*g(Pc,K_p)-V_2P*g(Pcp,K_dp)-v_dPC*g(Pcp,K_d)-k_dn*Pcp    
dCcp = V_1C*g(Cc,K_p)-V_2C*g(Ccp,K_dp)-v_dCC*g(Ccp,K_d)-k_dn*Ccp    
dPCc = -V_1PC*g(PCc,K_p)+V_2PC*g(PCcp,K_dp)-k4_clock*PCc+k3_clock*Pc*Cc+k2_clock*PCn-k1_clock*PCc-k_dn*PCc    
dPCn = -V_3PC*g(PCn,K_p)+V_4PC*g(PCnp,K_dp)-k2_clock*PCn+k1_clock*PCc-k7*Bn*PCn+k8*In-k_dn*PCn    
dPCcp = V_1PC*g(PCc,K_p)-V_2PC*g(PCcp,K_dp)-v_dPCC*g(PCcp,K_d)-k_dn*PCcp           
dPCnp = V_3PC*g(PCn,K_p)-V_4PC*g(PCnp,K_dp)-v_dPCN*g(PCnp,K_d)-k_dn*PCnp           
dBc = k_sB*Mbmal-V_1B*g(Bc,K_p)+V_2B*g(Bcp,K_dp)-k5*Bc+k6*Bn-k_dn*Bc               
dBcp = V_1B*g(Bc,K_p)-V_2B*g(Bcp,K_dp)-v_dBC*g(Bcp,K_d)-k_dn*Bcp              
dBn = -V_3B*g(Bn,K_p)+V_4B*g(Bnp,K_dp)+k5*Bc-k6*Bn-k7*Bn*PCn+k8*In-k_dn*Bn        
dBnp = V_3B*g(Bn,K_p)-V_4B*g(Bnp,K_dp)-v_dBN*g(Bnp,K_d)-k_dn*Bnp               
dIn = -k8*In+k7*Bn*PCn-v_dIN*g(In,K_d)-k_dn*In
dMr = v_sR*f(Bn,K_AR,h)-v_mR*g(Mr,K_mR)-k_dmr*Mr
dRc = k_sR*Mr-k9*Rc+k10*Rn-v_dRC*g(Rc,K_d)-k_dn*Rc
dRn = k9*Rc-k10*Rn-v_dRN*g(Rn,K_d)-k_dn*Rn

#cell cycle
dAP1 = (v_sap1*g(GF,K_agf)-k_dap1*AP1)*eps
dpRB = (v_sprb-k_pc1*pRB*E2F+k_pc2*pRBc1-V1*g(pRB,K1)*(Md+Mdp27)+V2*g(pRBp,K2)-k_dprb*pRB)*eps
dpRBc1 = (k_pc1*pRB*E2F-k_pc2*pRBc1)*eps
dpRBp = (V1*g(pRB,K1)*(Md+Mdp27)-V2*g(pRBp,K2)-V3*g(pRBp,K3)*Me+V4*g(pRBpp,K4)-k_pc3*pRBp*E2F+k_pc4*pRBc2-k_dprbp*pRBp)*eps
dpRBc2 = (k_pc3*pRBp*E2F-k_pc4*pRBc2)*eps
dpRBpp = (V3*g(pRBp,K3)*Me-V4*g(pRBpp,K4)-k_dprbpp*pRBpp)*eps
dE2F = (v_se2f-k_pc1*pRB*E2F+k_pc2*pRBc1-k_pc3*pRBp*E2F+k_pc4*pRBc2-V_1e2f*Ma*g(E2F,K_1e2f)+V_2e2f*g(E2Fp,K_2e2f)-k_de2f*E2F)*eps
dE2Fp = (V_1e2f*Ma*g(E2F,K_1e2f)-V_2e2f*g(E2Fp,K_2e2f)-k_de2fp*E2Fp)*eps

# Module Cyclin D/Cdk4-6 : G1 phase
dCd = (k_cd1*AP1+k_cd2*E2F*g(K_i7,pRB)*g(K_i8,pRBp)-k_com1*Cd*(Cdk4_tot-(Mdi+Md+Mdp27))+k_decom1*Mdi-V_dd*g(Cd,K_dd)-k_ddd*Cd)*eps
dMdi = (k_com1*Cd*(Cdk4_tot-(Mdi+Md+Mdp27))-k_decom1*Mdi+V_m2d*g(Md,K_2d)-V_m1d*g(Mdi,K_1d))*eps
dMd = (V_m1d*g(Mdi,K_1d)-V_m2d*g(Md,K_2d)-k_c1*Md*p27+k_c2*Mdp27)*eps
dMdp27 = (k_c1*Md*p27-k_c2*Mdp27)*eps

#Module Cyclin E/Cdk2: G1 phase and transition G1/S
dMce = u*v_sce*f(K_ice,Bn,nce)-V_dmce*g(Mce,K_dmce) 
dCe = (k_ce*E2F*(K_i9/(K_i9+pRB))*(K_i10/(K_i10+pRBp))+k_ce2*Mce-k_com2*Ce*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))+k_decom2*Mei-V_de*g(Skp2,K_dceskp2)*g(Ce,K_de)-k_dde*Ce)*eps
dMei = (k_com2*Ce*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))-k_decom2*Mei+V_m2e*(Wee1+i_b1)*g(Me,K_2e)-V_m1e*Pe*g(Mei,K_1e))*eps
dMe = (V_m1e*Pe*g(Mei,K_1e)-V_m2e*(Wee1+i_b1)*g(Me,K_2e)-k_c3*Me*p27+k_c4*Mep27)*eps
dSkp2 = (v_sskp2-V_dskp2*g(Skp2,K_dskp2)*g(Cdh1a,K_cdh1)-k_ddskp2*Skp2)*eps
dMep27 = (k_c3*Me*p27-k_c4*Mep27)*eps
dPei = (v_spei+V_6e*(x_e1+x_e2*Chk1)*g(Pe,K_6e)-V_m5e*(Me+a_e)*g(Pei,K_5e)-k_dpei*Pei)*eps
dPe = (V_m5e*(Me+a_e)*g(Pei,K_5e)-V_6e*(x_e1+x_e2*Chk1)*g(Pe,K_6e)-k_dpe*Pe)*eps

# Module Cyclin A/Cdk2 : S phase and transition S/G2
dCa = (k_ca*E2F*K_i11/(K_i11+pRB)*K_i12/(K_i12+pRBp)-k_com3*Ca*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))+
k_decom3*Mai-V_da*g(Ca,K_da)*g(Cdc20a,K_acdc20)-k_dda*Ca)*eps

dMai = (k_com3*Ca*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))-k_decom3*Mai+V_m2a*(Wee1+i_b2)*g(Ma,K_2a)-V_m1a*Pa*g(Mai,K_1a))*eps

dMa = (V_m1a*Pa*g(Mai,K_1a)-V_m2a*(Wee1+i_b2)*g(Ma,K_2a)-k_c5*Ma*p27+k_c6*Map27)*eps
dMap27 = (k_c5*Ma*p27-k_c6*Map27)*eps
dp27 = (v_s1p27+v_s2p27*E2F*g(K_i13,pRB)*g(K_i14,pRBp)-k_c1*Md*p27+k_c2*Mdp27-k_c3*Me*p27+k_c4*Mep27-k_c5*Ma*p27+k_c6*Map27-k_c7*Mb*p27+k_c8*Mbp27-V_1p27*Me*g(p27,K_1p27)+V_2p27*g(p27p,K_2p27)-k_ddp27*p27)*eps
dp27p = (V_1p27*Me*g(p27,K_1p27)-V_2p27*g(p27p,K_2p27)-V_dp27p*g(Skp2,K_dp27skp2)*g(p27p,K_dp27p)-k_ddp27p*p27p)*eps
dCdh1i = (V_2cdh1*g(Cdh1a,K_2cdh1)*(Ma+Mb)-V_1cdh1*g(Cdh1i,K_1cdh1)-k_dcdh1i*Cdh1i)*eps
dCdh1a = (v_scdh1a+V_1cdh1*g(Cdh1i,K_1cdh1)-V_2cdh1*g(Cdh1a,K_2cdh1)*(Ma+Mb)-k_dcdh1a*Cdh1a)*eps
dPai = (v_spai+V_6a*(x_a1+x_a2*Chk1)*g(Pa,K_6a)-V_m5a*(Ma+a_a)*g(Pai,K_5a)-k_dpai*Pai)*eps
dPa = (V_m5a*(Ma+a_a)*g(Pai,K_5a)-V_6a*(x_a1+x_a2*Chk1)*g(Pa,K_6a)-k_dpa*Pa)*eps

# Module Cyclin B/Cdk1 : G2 phase and transition G2/M

dCb = (v_cb-k_com4*Cb*(Cdk1_tot-(Mbi+Mb+Mbp27))+k_decom4*Mbi-V_db*g(Cb,K_db)*(g(Cdc20a,K_dbcdc20)+g(Cdh1a,K_dbcdh1))-k_ddb*Cb)*eps
dMbi = (k_com4*Cb*(Cdk1_tot-(Mbi+Mb+Mbp27))-k_decom4*Mbi+V_m2b*(Wee1+i_b3)*g(Mb,K_2b)-V_m1b*Pb*g(Mbi,K_1b))*eps
dMb = (V_m1b*Pb*g(Mbi,K_1b)-V_m2b*(Wee1+i_b3)*g(Mb,K_2b)-k_c7*Mb*p27+k_c8*Mbp27)*eps
dMbp27 = (k_c7*Mb*p27-k_c8*Mbp27)*eps
dCdc20i = (v_scdc20i-V_m3b*Mb*g(Cdc20i,K_3b)+V_m4b*g(Cdc20a,K_4b)-k_dcdc20i*Cdc20i)*eps
dCdc20a = (V_m3b*Mb*g(Cdc20i,K_3b)-V_m4b*g(Cdc20a,K_4b)-k_dcdc20a*Cdc20a)*eps
dPbi = (v_spbi+V_6b*(x_b1+x_b2*Chk1)*g(Pb,K_6b)-V_m5b*(Mb+a_b)*g(Pbi,K_5b)-k_dpbi*Pbi)*eps
dPb = (V_m5b*(Mb+a_b)*g(Pbi,K_5b)-V_6b*(x_b1+x_b2*Chk1)*g(Pb,K_6b)-k_dpb*Pb)*eps

#coupling

dMw = v_swee1+u*v_sw*f(Bn,K_aw,nmw)-V_dmw*g(Mw,K_dmw)
dWee1 = (k_sw*Mw-V_m7b*(Mb+i_b)*g(Wee1,K_7b)+V_m8b*g(Wee1p,K_8b)-k_dwee1*Wee1)*eps
dWee1p = (V_m7b*(Mb+i_b)*g(Wee1,K_7b)-V_m8b*g(Wee1p,K_8b)-k_dwee1p*Wee1p)*eps

# @show X[5]

dX .= dMp, dMc, dMbmal, dPc, dCc, dPcp, dCcp, dPCc, dPCn, dPCcp, dPCnp, dBc, dBcp, dBn, dBnp, dIn, dMr, dRc, dRn, dAP1, dpRB, dpRBc1, dpRBp, dpRBc2, dpRBpp, dE2F, dE2Fp, dCd, dMdi, dMd, dMdp27, dMce, dCe, dMei, dMe, dSkp2, dMep27, dPei, dPe, dCa, dMai, dMap27, dp27, dp27p, dCdh1i, dCdh1a, dPai, dPa, dCb, dMbi, dMb, dMbp27, dCdc20i, dCdc20a, dPbi, dPb, dMw, dWee1, dWee1p 

    return nothing
end


Mp=1.1569999; Mc=1.7789; Mbmal=7.2195001; Pc=0.0062671001; Cc=325.28311; Pcp=0.0034558999; Ccp=0.76152998; PCc=0.70789999; PCn=0.19575; PCcp=5.2571998; PCnp=0.12734; Bc=2.1372001; Bcp=0.11593; Bn=1.3151; Bnp=0.076318003;  In=0.055312; Mr=0.96881002; Rc=0.1548; Rn=0.041343; AP1=6.0605998; pRB=1.5329; pRBc1=0.48471999; pRBp=12.9828; pRBc2=2.0545001; pRBpp=0.13319001; E2F=3.2801001; E2Fp=0.026215; Cd=0.094091997; Mdi=0.022845; Md=1.325; Mdp27=0.013466; Mce=0.024901001; Ce=0.005537; Mei=0.031500999; Me=1.1559; Skp2=11.8677; Mep27=0.0081121; Pei=0.014376; Pe=1.7049;  Ca=0.0046791998; Mai=0.046822; Ma=0.019317999; Map27=0.0002039; p27=0.0035139001; p27p=0.019300999; Cdh1i=0.54710001; Cdh1a=0.0055264998; Pai=0.63407999; Pa=0.29585001; Cb=0.95919001; Mbi=0.034823999; Mb=0.45019999; Mbp27=0.00094055; Cdc20i=0.026786; Cdc20a=1.8816; Pbi=0.057624999; Pb=1.0843; Mw=0.016368; Wee1=0.13579001; Wee1p=0.32049

u0 = [Mp, Mc, Mbmal, Pc, Cc, Pcp, Ccp, PCc, PCn, PCcp, PCnp, Bc, Bcp, Bn, Bnp, In, Mr, Rc, Rn, AP1, pRB, pRBc1, pRBp, pRBc2, pRBpp, E2F, E2Fp, Cd, Mdi, Md, Mdp27, Mce, Ce, Mei, Me, Skp2, Mep27, Pei, Pe, Ca, Mai, Map27, p27, p27p, Cdh1i, Cdh1a, Pai, Pa, Cb, Mbi, Mb, Mbp27, Cdc20i, Cdc20a, Pbi, Pb, Mw, Wee1, Wee1p]

# 3 -> Mbmal 
# 33 -> Ce : CyclinE/Cdk2
# 49 -> Cb : CyclinB/Cdk1 
diffeq = (reltol = 1e-6,  alg = Tsit5())
df = CoupledODEs(circadian_cell_cycle!, rand(59); diffeq)
# _complete(y) = (length(y) == 2) ? rand(59) : y; 
# proj = ProjectedDynamicalSystem(df, [33, 49], u0[[1:32; 34:48; 50:59]])
# yg = xg   = range(0., 3., length = 2001)
# mapper = AttractorsViaRecurrences(proj, grid; Δt = 1.0)

Mp=1.0568; Mc=1.8480999; Mbmal=6.3362999; Pc=0.0063677998; Cc=331.87189; Pcp=0.0036138999; Ccp=0.76160997; PCc=1.1508; PCn=0.53987998; PCcp=6.4289999; PCnp=0.9738; Bc=1.7359; Bcp=0.10505; Bn=0.78301001; Bnp=0.054529;  In=0.13778999; Mr=0.80400997; Rc=0.12732001; Rn=0.056669999; AP1=6.0605998; pRB=1.6623; pRBc1=0.064092003; pRBp=14.2593; pRBc2=0.27438; pRBpp=0.073599003; E2F=0.38652; E2Fp=0.13256; Cd=0.094083004; Mdi=0.022838; Md=1.3206; Mdp27=0.017968001; Mce=0.043490998; Ce=0.0096067004; Mei=0.029376; Me=0.85652; Skp2=10.3268; Mep27=0.0081107998; Pei=0.01898; Pe=1.6959;  Ca=0.0052708001; Mai=0.058660999; Ma=0.57453001; Map27=0.0034302; p27=0.0047358; p27p=0.019261001; Cdh1i=0.54900002; Cdh1a=0.001937; Pai=0.049509; Pa=1.3054; Cb=1.5501; Mbi=0.034116; Mb=0.45583999; Mbp27=0.0012757001; Cdc20i=0.026214; Cdc20a=1.6991; Pbi=0.056650002; Pb=1.0693001; Mw=0.012853; Wee1=0.11596; Wee1p=0.26273999; 

u1 = [Mp, Mc, Mbmal, Pc, Cc, Pcp, Ccp, PCc, PCn, PCcp, PCnp, Bc, Bcp, Bn, Bnp, In, Mr, Rc, Rn, AP1, pRB, pRBc1, pRBp, pRBc2, pRBpp, E2F, E2Fp, Cd, Mdi, Md, Mdp27, Mce, Ce, Mei, Me, Skp2, Mep27, Pei, Pe, Ca, Mai, Map27, p27, p27p, Cdh1i, Cdh1a, Pai, Pa, Cb, Mbi, Mb, Mbp27, Cdc20i, Cdc20a, Pbi, Pb, Mw, Wee1, Wee1p]

diffeq = (reltol = 1e-6,  alg = Tsit5())
df = CoupledODEs(circadian_cell_cycle!, rand(59); diffeq)
y,t = trajectory(df, 40., u1)

lines(y[:,3],y[:,33])
