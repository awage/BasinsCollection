using DrWatson
@quickactivate
using OrdinaryDiffEqVerner
using Attractors
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))


function circadian_cell_cycle!(du, u, p, t)
# Variables
    Mp, Mc, Mbmal, Pc, Cc, Pcp, Ccp, PCc, PCn, PCcp, PCnp, Bc, Bcp, Bn, Bnp, In, Mr, Rc, Rcp, Rn, Rnp, AP1, pRB, pRBc1, pRBp, pRBc2, pRBpp, E2F, E2Fp, Cd, Mdi, Md, Mdp27, Mce, Ce, Mei, Me, Skp2, Mep27, Pei, Pe, Ca, Mai, Ma, Map27, p27, p27p, Cdh1i, Cdh1a, Pai, Pa, Cb, Mbi, Mb, Mbp27, Cdc20i, Cdc20a, Pbi, Pb, Mw, Wee1, Wee1p = u

### CONSTANTS
V_cdk = 0.; KI_cdk1 = 0.8; ncdk = 1.; v_sw = p[1]
    k1_clock = 0.8; k2_clock = 0.4; k3_clock = 0.8; k4_clock = 0.4; k5 = 0.8; k6 = 0.4; k7 = 1; k8 = 0.2; k9 = 0.63; k10 = 0.4; 
    K_AP = 0.6; K_AC = 0.6; K_AR = 0.6; K_IB = 1; 
    k_dmb = 0.02; k_dmc = 0.02; k_dmp = 0.02; k_dmr = 0.02; k_dn = 0.02; k_dnc = 0.02; 
    K_d = 0.3; K_dp = 0.1; K_p = 1.006; K_mB = 0.4; K_mC = 0.4; K_mP = 0.3; K_mR = 0.4; k_sB = 0.32; 
    k_sC = 3.2; k_sP = 1.2; k_sR = 1.7; m = 2; h = 2; n = 2; 
    V_1B = 1.4; V_1C = 1.2; V_1P = 9.6; V_1PC = 2.4; V_2B = 0.2; V_2C = 0.2; V_2P = 0.6; V_2PC = 0.2; 
    V_3B = 1.4; V_3PC = 2.4; V_4B = 0.4; V_4PC = 0.2; V_phos = 0.4; v_dBC = 3; v_dBN = 3; v_dCC = 1.4; 
    v_dIN = 1.6; v_dPC = 3.4; v_dPCC = 1.4; v_dPCN = 1.4; v_dRC = 4.4; v_dRN = 0.8; v_mB = 1.3; 
    v_mC = 2.0; v_mP = 2.2; v_mR = 1.6; v_sB = 1.8; 
    v_sC = 2.2; v_sP = 2.4; v_sR = 1.6; 
    V_1R = 4; V_2R = 8; V_3R = 8; V_4R = 4; v_in = 0.7
    Chk1 = 0

    delta = 1.
    eps = 22.5

    GF = 1; K_agf = 0.1; k_dap1 = 0.15;  v_sap1 = 1
    k_de2f = 0.002; k_de2fp = 1.1; k_dprb = 0.01; k_dprbp = 0.06; k_dprbpp = 0.04
    k_pc1 = 0.05; k_pc2 = 0.5; k_pc3 = 0.025; k_pc4 = 0.5; K1 = 0.1; K2 = 0.1; K3 = 0.1
    K4 = 0.1; V1 = 2.2; V2 = 2; V3 = 1; V4 = 2; K_1e2f = 5; K_2e2f = 5; V_1e2f = 4
    V_2e2f = 0.75; v_se2f = 0.15; v_sprb = 0.8
    Cdk4_tot = 1.5; K_i7 = 0.1; K_i8 = 2; k_cd1 = 0.4; k_cd2 = 0.005; k_decom1 = 0.1
    k_com1 = 0.175; k_c1 = 0.15; k_c2 = 0.05; k_ddd = 0.005; K_dd = 0.1; K_1d = 0.1; K_2d = 0.1
    V_dd = 5; V_m1d = 1; V_m2d = 0.2
    a_e = 0.25; Cdk2_tot = 2; i_b1 = 0.5; K_i9 = 0.1; K_i10 = 2; k_ce = 0.29; k_c3 = 0.2
    k_c4 = 0.1; k_decom2 = 0.1; k_com2 = 0.2; k_dde = 0.005; k_ddskp2 = 0.005; k_dpe = 0.075
    k_dpei = 0.15; K_de = 0.1; K_dceskp2 = 2; K_dskp2 = 0.5; K_cdh1 = 0.4; K_1e = 0.1
    K_2e = 0.1; K_5e = 0.1; K_6e = 0.1; V_de = 3; V_dskp2 = 1.1; V_m1e = 2; V_m2e = 1.4; V_m5e = 5
    V_6e = 0.8; v_spei = 0.13; v_sskp2 = 0.15; x_e1 = 1; x_e2 = 1
    a_a = 0.2; i_b2 = 0.5; K_i11 = 0.1; K_i12 = 2; K_i13 = 0.1; K_i14 = 2; k_ca = 0.0375
    k_decom3 = 0.1; k_com3 = 0.2; k_c5 = 0.15; k_c6 = 0.125; k_dda = 0.005; k_ddp27 = 0.06
    k_ddp27p = 0.01; k_dcdh1a = 0.1; k_dcdh1i = 0.2; k_dpa = 0.075; k_dpai = 0.15; K_da = 1.1
    K_dp27p = 0.1; K_dp27skp2 = 0.1; K_acdc20 = 2; K_1a = 0.1; K_2a = 0.1; K_1cdh1 = 0.01
    K_2cdh1 = 0.01; K_5a = 0.1; K_6a = 0.1; K_1p27 = 0.5; K_2p27 = 0.5; V_dp27p = 5; V_da = 2.5
    V_m1a = 2; V_m2a = 1.85; V_m5a = 4; V_6a = 1; v_scdh1a = 0.11; v_spai = 0.105
    v_s1p27 = 0.8; v_s2p27 = 0.1; V_1cdh1 = 1.25; V_2cdh1 = 8; V_1p27 = 100; V_2p27 = 0.1
    x_a1 = 1; x_a2 = 1
    a_b = 0.11; Cdk1_tot = 0.5; i_b = 0.75; i_b3 = 0.5; k_c7 = 0.12; k_c8 = 0.2
    k_decom4 = 0.1; k_com4 = 0.25; k_dcdc20a = 0.05; k_dcdc20i = 0.14; k_ddb = 0.005
    k_dpb = 0.1; k_dpbi = 0.2; k_dwee1 = 0.1; k_dwee1p = 0.2; K_db = 0.005; K_dbcdc20 = 0.2; K_dbcdh1 = 0.1
    k_sw = 5; K_1b = 0.1; K_2b = 0.1; K_3b = 0.1; K_4b = 0.1; K_5b = 0.1; K_6b = 0.1; K_7b = 0.1
    K_8b = 0.1; v_cb = 0.055; V_db = 0.06; V_m1b = 3.9; V_m2b = 2.1; v_scdc20i = 0.1; V_m3b = 8; V_m4b = 0.7
    V_m5b = 5; V_6b = 1; V_m7b = 1.2; V_m8b = 1; v_spbi = 0.12; x_b1 = 1; x_b2 = 1
    v_swee1 = 0.0117; nmw = 4; K_aw = 2; V_dmw = 0.5; K_dmw = 0.5
    v_sce = 0.000; K_ice = 1; V_dmce = 0.5; K_dmce = 0.5; nce = 4; k_ce2 = 5

    f(x,y,p) = x^p/(x^p + y^p) 
    g(x,y) = x/(x+y)
#circadian clock
    Mp_d = (v_sP*f(Bn,K_AP,n) + v_in*f(KI_cdk1,Mb,ncdk) - v_mP*g(Mp,K_mP) - k_dmp*Mp)*delta
    Mc_d = (v_sC*f(Bn,K_AC,n) + v_in*f(KI_cdk1,Mb,ncdk) - v_mC*g(Mc,K_mC) - k_dmc*Mc)*delta
    Mbmal_d = (v_sB*f(K_IB, Rn, m) + v_in*f(KI_cdk1, Mb, ncdk) - v_mB*g(Mbmal,K_mB) -k_dmb*Mbmal)*delta
    Pc_d = (k_sP*Mp - V_1P*g(Pc,K_p) + V_2P*g(Pcp,K_dp)  + k4_clock*PCc - k3_clock*Pc*Cc - k_dn*Pc)*delta
    Cc_d = (k_sC*Mc - V_1C*g(Cc,K_p) + V_2C*g(Ccp,K_dp) + k4_clock*PCc - k3_clock*Pc*Cc - k_dnc*Cc)*delta
    Pcp_d = (V_1P*g(Pc, K_p) - V_2P*g(Pcp,K_dp) - v_dPC*g(Pcp,K_d) - k_dn*Pcp)*delta
    Ccp_d = (V_1C*g(Cc,K_p) - V_2C*g(Ccp, K_dp) - v_dCC*g(Ccp, K_d) - k_dn*Ccp)*delta
    PCc_d = (-V_1PC*g(PCc,K_p) + V_2PC*g(PCcp,K_dp) - k4_clock*PCc + k3_clock*Pc*Cc + k2_clock*PCn - k1_clock*PCc - k_dn*PCc)*delta
    PCn_d = (-V_3PC*g(PCn,K_p) +V_4PC*g(PCnp,K_dp) -k2_clock*PCn + k1_clock*PCc - k7*Bn*PCn + k8*In - k_dn*PCn)*delta
    PCcp_d = (V_1PC*g(PCc,K_p) -V_2PC*g(PCcp,K_dp) -v_dPCC*g(PCcp,K_d) -k_dn*PCcp)*delta
    PCnp_d = (V_3PC*g(PCn,K_p) -V_4PC*g(PCnp,K_dp) - v_dPCN*g(PCnp,K_d) - k_dn*PCnp)*delta
    Bc_d = (k_sB*Mbmal - (V_1B)*g(Bc,K_p) +V_2B*g(Bcp,K_dp) - k5*Bc + k6*Bn - k_dn*Bc)*delta
    Bcp_d = ((V_1B)*g(Bc,K_p) - V_2B*g(Bcp,K_dp) - v_dBC*g(Bcp,K_d) - k_dn*Bcp)*delta
    Bn_d = (-V_3B*g(Bn,K_p) +V_4B*g(Bnp,K_dp) + k5*Bc - k6*Bn - k7*Bn*PCn + k8*In - k_dn*Bn)*delta
    Bnp_d = (V_3B*g(Bn,K_p) - V_4B*g(Bnp,K_dp) - v_dBN*g(Bnp,K_d) - k_dn*Bnp)*delta
    In_d = (-k8*In + k7*Bn*PCn - v_dIN*g(In,K_d) - k_dn*In)*delta
    Mr_d = (v_sR*f(Bn,K_AR,h) + v_in*f(KI_cdk1,Mb,ncdk) - v_mR*g(Mr,K_mR) - k_dmr*Mr)*delta
    Rc_d = (k_sR*Mr - k9*Rc + k10*Rn - (V_1R+V_cdk*Mb)*g(Rc,K_p) + V_2R*g(Rcp,K_dp) - k_dn*Rc)*delta
    Rcp_d = ((V_1R+V_cdk*Mb)*g(Rc,K_p) -V_2R*g(Rcp,K_dp) -v_dRC*g(Rcp,K_d) - k_dn*Rcp)*delta
    Rn_d = (k9*Rc - k10*Rn - (V_3R+V_cdk*Mb)*g(Rn,K_p) +V_4R*g(Rnp,K_dp) - k_dn*Rn)*delta
    Rnp_d = ((V_3R+V_cdk*Mb)*g(Rn,K_p) -V_4R*g(Rnp,K_dp) -v_dRN*g(Rnp,K_d) - k_dn*Rnp)*delta


#cell cycle
# Mitotic stimulation by growth factor, GF
    AP1_d = (v_sap1*g(GF,K_agf) -k_dap1*AP1)*eps
# Antagonistic regulation exerted by pRB and E2F
    pRB_d = (v_sprb-k_pc1*pRB*E2F+k_pc2*pRBc1-V1*g(pRB,K1)*(Md+Mdp27)+V2*g(pRBp,K2) -k_dprb*pRB)*eps
    pRBc1_d = (k_pc1*pRB*E2F-k_pc2*pRBc1)*eps
    pRBp_d = (V1*g(pRB,K1)*(Md+Mdp27)-V2*g(pRBp,K2) - V3*g(pRBp,K3) *Me+V4*g(pRBpp,K4) -k_pc3*pRBp*E2F+k_pc4*pRBc2-k_dprbp*pRBp)*eps
    pRBc2_d = (k_pc3*pRBp*E2F-k_pc4*pRBc2)*eps
    pRBpp_d = (V3*g(pRBp,K3)*Me-V4*g(pRBpp,K4) -k_dprbpp*pRBpp)*eps
    E2F_d = (v_se2f-k_pc1*pRB*E2F+k_pc2*pRBc1-k_pc3*pRBp*E2F+k_pc4*pRBc2*V_1e2f*Ma*g(E2F,K_1e2f) +V_2e2f*g(E2Fp,K_2e2f) - k_de2f*E2F)*eps
    E2Fp_d = (V_1e2f*Ma*g(E2F,K_1e2f) -V_2e2f*g(E2Fp,K_2e2f) -k_de2fp*E2Fp)*eps
# Module Cyclin D/Cdk4-6 : G1 phase
    Cd_d = (k_cd1*AP1+k_cd2*E2F*g(K_i7,K_i7)*g(K_i8,K_i8) -k_com1*Cd*(Cdk4_tot-(Mdi+Md+Mdp27))+k_decom1*Mdi-V_dd*g(Cd,K_dd) -k_ddd*Cd)*eps
    Mdi_d = (k_com1*Cd*(Cdk4_tot-(Mdi+Md+Mdp27))-k_decom1*Mdi+V_m2d*g(Md,K_2d) -V_m1d*g(Mdi,K_1d) )*eps
    Md_d = (V_m1d*g(Mdi,K_1d) -V_m2d*g(Md,K_2d) -k_c1*Md*p27+k_c2*Mdp27)*eps
    Mdp27_d = (k_c1*Md*p27-k_c2*Mdp27)*eps
#Module Cyclin E/Cdk2: G1 phase and transition G1/S
    Mce_d = v_sce*f(K_ice,Bn,nce) - V_dmce*g(Mce,K_dmce) 
    Ce_d = (k_ce*E2F*g(K_i9,K_i9)*g(K_i10,K_i10) + k_ce2*Mce-k_com2*Ce*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))+k_decom2*Mei-V_de*g(Skp2,K_dceskp2)*g(Ce,K_de) -k_dde*Ce)*eps
    Mei_d = (k_com2*Ce*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))-k_decom2*Mei+V_m2e*(Wee1+i_b1)*g(Me,K_2e) - V_m1e*Pe*g(Mei,K_1e))*eps
    Me_d = (V_m1e*Pe*g(Mei,K_1e) -V_m2e*(Wee1+i_b1)*g(Me,K_2e) -k_c3*Me*p27 + k_c4*Mep27)*eps
    Skp2_d = (v_sskp2-V_dskp2*g(Skp2,K_dskp2)*g(Cdh1a,K_cdh1) -k_ddskp2*Skp2)*eps
    Mep27_d = (k_c3*Me*p27-k_c4*Mep27)*eps
    Pei_d = (v_spei+V_6e*(x_e1+x_e2*Chk1)*g(Pe,K_6e) -V_m5e*(Me+a_e)*g(Pei,K_5e) -k_dpei*Pei)*eps
    Pe_d = (V_m5e*(Me+a_e)*g(Pei,K_5e) -V_6e*(x_e1+x_e2*Chk1)*g(Pe,K_6e) -k_dpe*Pe)*eps

# Module Cyclin A/Cdk2 : S phase and transition S/G2
    Ca_d = (k_ca*E2F*g(K_i11,K_i11)*g(K_i12,K_i12) -k_com3*Ca*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))+k_decom3*Mai-V_da*g(Ca,K_da)*g(Cdc20a,K_acdc20) -k_dda*Ca)*eps
    Mai_d = (k_com3*Ca*(Cdk2_tot-(Mei+Me+Mep27+Mai+Ma+Map27))-k_decom3*Mai+V_m2a*(Wee1+i_b2)*g(Ma,K_2a) -V_m1a*Pa*g(Mai,K_1a))*eps
    Ma_d = (V_m1a*Pa*g(Mai,K_1a) -V_m2a*(Wee1+i_b2)*g(Ma,K_2a) -k_c5*Ma*p27+k_c6*Map27)*eps
    Map27_d = (k_c5*Ma*p27-k_c6*Map27)*eps
    p27_d = (v_s1p27+v_s2p27*E2F*g(K_i13,K_i13) *g(K_i14,K_i14) -k_c1*Md*p27+k_c2*Mdp27-k_c3*Me*p27+k_c4*Mep27-k_c5*Ma*p27+k_c6*Map27-k_c7*Mb*p27+k_c8*Mbp27-V_1p27*Me*g(p27,K_1p27) +V_2p27*g(p27p,K_2p27) -k_ddp27*p27)*eps
    p27p_d = (V_1p27*Me*g(p27,K_1p27) -V_2p27*g(p27p,K_2p27) -V_dp27p*g(Skp2,K_dp27skp2)*g(p27p,K_dp27p) - k_ddp27p*p27p)*eps
    Cdh1i_d = (V_2cdh1*g(Cdh1a,K_2cdh1)*(Ma+Mb)-V_1cdh1*g(Cdh1i,K_1cdh1) -k_dcdh1i*Cdh1i)*eps
    Cdh1a_d = (v_scdh1a+V_1cdh1*g(Cdh1i,K_1cdh1) -V_2cdh1*g(Cdh1a,K_2cdh1) *(Ma+Mb)-k_dcdh1a*Cdh1a)*eps
    Pai_d = (v_spai+V_6a*(x_a1+x_a2*Chk1)*g(Pa,K_6a) -V_m5a*(Ma+a_a)*g(Pai,K_5a) -k_dpai*Pai)*eps
    Pa_d = (V_m5a*(Ma+a_a)*g(Pai,K_5a) -V_6a*(x_a1+x_a2*Chk1)*g(Pa,K_6a) -k_dpa*Pa)*eps

# Module Cyclin B/Cdk1 : G2 phase and transition G2/M
    Cb_d = (v_cb-k_com4*Cb*(Cdk1_tot-(Mbi+Mb+Mbp27))+k_decom4*Mbi-V_db*g(Cb,K_db) *(g(Cdc20a,K_dbcdc20) + g(Cdh1a,K_dbcdh1) )-k_ddb*Cb)*eps
    Mbi_d = (k_com4*Cb*(Cdk1_tot-(Mbi+Mb+Mbp27))-k_decom4*Mbi+V_m2b*(Wee1+i_b3)*g(Mb,K_2b) -V_m1b*Pb*g(Mbi,K_1b))*eps
    Mb_d = (V_m1b*Pb*g(Mbi,K_1b) -V_m2b*(Wee1+i_b3)*g(Mb,K_2b) -k_c7*Mb*p27+k_c8*Mbp27)*eps
    Mbp27_d = (k_c7*Mb*p27-k_c8*Mbp27)*eps
    Cdc20i_d = (v_scdc20i-V_m3b*Mb*g(Cdc20i,K_3b) +V_m4b*g(Cdc20a,K_4b)  - k_dcdc20i*Cdc20i)*eps
    Cdc20a_d = (V_m3b*Mb*g(Cdc20i,K_3b) -V_m4b*g(Cdc20a,K_4b) -k_dcdc20a*Cdc20a)*eps
    Pbi_d = (v_spbi+V_6b*(x_b1+x_b2*Chk1)*g(Pb,K_6b) -V_m5b*(Mb+a_b)*g(Pbi,K_5b) -k_dpbi*Pbi)*eps
    Pb_d = (V_m5b*(Mb+a_b)*g(Pbi,K_5b) -V_6b*(x_b1+x_b2*Chk1)*g(Pb,K_6b) -k_dpb*Pb)*eps

#coupling via wee1
    Mw_d = v_swee1+ v_sw*f(Bn,K_aw,nmw) - V_dmw*g(Mw,K_dmw) 
    Wee1_d = (k_sw*Mw - V_m7b*(Mb+i_b)*g(Wee1,K_7b) +V_m8b*g(Wee1p,K_8b) -k_dwee1*Wee1)*eps
    Wee1p_d = (V_m7b*(Mb+i_b)*g(Wee1,K_7b) -V_m8b*g(Wee1p,K_8b) -k_dwee1p*Wee1p)*eps


    du .= Mp_d, Mc_d, Mbmal_d, Pc_d, Cc_d, Pcp_d, Ccp_d, PCc_d, PCn_d, PCcp_d, PCnp_d, Bc_d, Bcp_d, Bn_d, Bnp_d, In_d, Mr_d, Rc_d, Rcp_d, Rn_d, Rnp_d, AP1_d, pRB_d, pRBc1_d, pRBp_d, pRBc2_d, pRBpp_d, E2F_d, E2Fp_d, Cd_d, Mdi_d, Md_d, Mdp27_d, Mce_d, Ce_d, Mei_d, Me_d, Skp2_d, Mep27_d, Pei_d, Pe_d, Ca_d, Mai_d, Ma_d, Map27_d, p27_d, p27p_d, Cdh1i_d, Cdh1a_d, Pai_d, Pa_d, Cb_d, Mbi_d, Mb_d, Mbp27_d, Cdc20i_d, Cdc20a_d, Pbi_d, Pb_d, Mw_d, Wee1_d, Wee1p_d

# @show length(du)
    return nothing
end

Mp=0.1; Mc=0.1; Mbmal=0.1; Pc=0.1; Cc=0.1; 
Pcp=0.1; Ccp=0.1; PCc=0.1; PCn=0.1; PCcp=0.1
PCnp=0.1; Bc=0.1; Bcp=0.1; Bn=0.1; Bnp=0.1; 
In=0.1; Mr=0.1; Rc=0.1; Rcp=0.1; Rn=0.1; Rnp=0.1
AP1=0.01; pRB=1; pRBc1=0.25; pRBp=0.1; pRBc2=0.01; pRBpp=0.01; E2F=0.1; E2Fp=0.05
Cd=0.01; Mdi=0.01; Md=0.01; Mdp27=0.01; 
Mce=0.1; Ce=0.01; Mei=0.01; Me=0.01; Skp2=0.01; Mep27=0.01; Pei=0.01; Pe=0.01
Ca=0.01; Mai=0.01; Ma=0.01; Map27=0.01; p27=0.25; p27p=0.01; Cdh1i=0.01; Cdh1a=0.01; Pai=0.01; Pa=0.01
Cb=0.01; Mbi=0.01; Mb=0.01; Mbp27=0.01; Cdc20i=0.01; Cdc20a=0.01; Pbi=0.01; Pb=0.01; Mw=0; Wee1=0.1
Wee1p=0.01


u0 = [Mp, Mc, Mbmal, Pc, Cc, Pcp, Ccp, PCc, PCn, PCcp, PCnp, Bc, Bcp, Bn, Bnp, In, Mr, Rc, Rcp, Rn, Rnp, AP1, pRB, pRBc1, pRBp, pRBc2, pRBpp, E2F, E2Fp, Cd, Mdi, Md, Mdp27, Mce, Ce, Mei, Me, Skp2, Mep27, Pei, Pe, Ca, Mai, Ma, Map27, p27, p27p, Cdh1i, Cdh1a, Pai, Pa, Cb, Mbi, Mb, Mbp27, Cdc20i, Cdc20a, Pbi, Pb, Mw, Wee1, Wee1p]

# 14 -> Bn
# 54 -> Mb 
# 37 -> Me
diffeq = (reltol = 1e-6,  alg = Tsit5())
df = CoupledODEs(circadian_cell_cycle!, rand(62), [2.2]; diffeq)
_complete(y) = (length(y) == 3) ? rand(62) : y; 
proj = ProjectedDynamicalSystem(df, [14, 37, 54], u0[[1:13; 15:36; 38:53; 55:62]])
# proj = ProjectedDynamicalSystem(df, [14, 37, 54], _complete)
# yg = range(-1, 3; length = 1001)
# grid = ntuple(x -> yg, dimension(df))
pow = 3; xg = range(0, 2^(1/pow); length = 20001).^pow
# yg = xg = zg = range(0., 3., length = 2001)
mapper = AttractorsViaRecurrences(proj, (xg,xg,xg); Δt = 1.0, consecutive_recurrences = 150)
# mapper = AttractorsViaRecurrences(df, grid; Δt = 1.0, consecutive_recurrences = 500)    

# for _ in 1:100
#     @show mapper(rand(3))
#     # @show mapper(rand(62)*0.4)
# end
#
# prob = ODEProblem(circadian_cell_cycle!, u0, (0,1.))
# sol = solve(prob, DP8(), abstol = 1e-12; dt = 0.01, adaptive = false)

# y,t = trajectory(df, 10., rand(62)*100; Δt = 0.01)
# y,t = trajectory(df, 400., rand(3); Ttr = 200, Δt = 0.01)
# plot(t, y[:,62])
# plot(y[:,1], y[:,2])

prange = 2:0.5:5
pidx = 1 # index of the parameter
xg = range(0,2,length =10)
grid = (xg,xg,xg)
sampler, = statespace_sampler(grid)
ascm = AttractorSeedContinueMatch(mapper)
fractions_cont, attractors_cont = global_continuation(
	ascm, prange, pidx, sampler; samples_per_parameter = 10
)
