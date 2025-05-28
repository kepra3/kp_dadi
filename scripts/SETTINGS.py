import demo_models_kp

# GLOBAL VARIABLES
SET_PTS = [100, 120, 130]


def get_settings(model, ALL=False):
    if model == "snm.1d":
        # standard neutral model 1d
        num = 1
        p_labels = "nu"
        upper = [200]
        lower = [0.001]
        model_fun = demo_models_kp.no_divergence_1d
    elif model == "size_change":
        num = 2
        p_labels = "nu, T"
        upper = [200, 200]
        lower = [0.001, 0.001]
        model_fun = demo_models_kp.instant_change
    elif model == "bottle":
        num = 3
        p_labels = "nuB, nuF, T"
        upper = [200, 200, 200]
        lower = [0.001, 0.001, 0.001]
        model_fun = demo_models_kp.bottlegrowth
    elif model == "bottle_neck":
        num = 4
        p_labels = "nuB, nuF, TB, TF"
        upper = [200, 200, 200, 200]
        lower = [0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.bottleneck
    elif model == "no_mig":
        # divergence with no migration
        num = 3
        p_labels = "nu1, nu2, T"
        upper = [150, 150, 15]
        lower = [0.001, 0.001, 0.001]
        model_fun = demo_models_kp.no_migration
    elif model == "snm":
        # standard neutral model, no divergence
        num = 1
        p_labels = "nu"
        upper = [150]
        lower = [0.001]
        model_fun = demo_models_kp.no_divergence
    elif model == "sym_mig":
        # divergence with symmetrical migration
        num = 4
        p_labels = "nu1, nu2, m, T"
        upper = [150, 150, 10, 15]
        lower = [0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sym_migration
    elif model == "asym_mig":
        # divergence with asymmetrical migration
        num = 5
        p_labels = "nu1, nu2, m12, m21, T"
        upper = [150, 150, 10, 10, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.asym_migration
    elif model == "anc_sym_mig":
        num = 5
        p_labels = "nu1, nu2, m, T1, T2"
        upper = [150, 150, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_sym_migration
    elif model == "anc_asym_mig":
        num = 6
        p_labels = "nu1, nu2, m12, m21, T1, T2"
        upper = [150, 150, 10, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_asym_migration
    elif model == "sec_cont_sym_mig":
        num = 5
        p_labels = "nu1, nu2, m, T1, T2"
        upper = [150, 150, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_sym_migration
    elif model == "sec_cont_asym_mig":
        num = 6
        p_labels = "nu1, nu2, m12, m21, T1, T2"
        upper = [150, 150, 10, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_asym_migration
    elif model == "iso_inbred":
        # divergence with inbreeding
        num = 5
        p_labels = "nu1, nu2, F1, F2, T"
        upper = [150, 150, 1, 1, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001]
        model_fun = demo_models_kp.iso_inbreeding
    elif model == "mig_inbred":
        # divergence with migration with inbreeding
        num = 6
        p_labels = "nu1, nu2, F1, F2, m, T"
        upper = [150, 150, 1, 1, 10, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001]
        model_fun = demo_models_kp.mig_inbreeding
    elif model == "anc_mig_inbred":
        # ancient migration with inbreeding
        num = 7
        p_labels = "nu1, nu2, F1, F2, m, T1, T2"
        upper = [150, 150, 1, 1, 10, 15, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_sym_mig_inbred
    elif model == "sec_cont_inbred":
        # secondary contact with inbreeding
        num = 7
        p_labels = "nu1, nu2, F1, F2, m, T1, T2"
        upper = [150, 150, 1, 1, 10, 15, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_sym_mig_inbred
    elif model == "mig_be_inbred":
        # population size changes with divergence with migration and inbreeding
        num = 10
        p_labels = "nu1, nu2, nu1a, nu2a, F1, F2, m1, m2,  T1, T2"
        upper = [150, 150, 150, 150, 1, 1, 10, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.mig_be_inbred
    elif model == "split_bottle":
        num = 6
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2"
        upper = [200, 200, 200, 200, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_bottlegrowth
    elif model == "split_bottle_asym_mig":
        num = 10
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2, m12T1, m21T1, m12T2, m21T2"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_bottlegrowth_asym_mig
    elif model == "split_bottle_sec_asym_mig":
        num = 8
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2, m12, m21"
        upper = [200, 200, 200, 200, 15, 15, 10, 10]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_bottlegrowth_second_asym_mig
    elif model == "split_bottle_anc_asym_mig":
        num = 8
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2, m12, m21"
        upper = [200, 200, 200, 200, 15, 15, 10, 10]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_bottlegrowth_ancient_asym_mig
    elif model == "split_sizechange":
        num = 6
        p_labels = "nu1T1, nu2T2, nu1T1, nu2T2, T1, T2"
        upper = [200, 200, 200, 200, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_sizechange
    elif model == "split_sizechange_asym_mig":
        num = 10
        p_labels = "nu1T1, nu2T2, nu1T1, nu2T2, T1, T2, m12T1, m21T1, m12T2, m21T2"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_sizechange_asym_mig
    elif model == "split_sizechange_anc_asym_mig":
        num = 8
        p_labels = "nu1T1, nu2T2, nu1T1, nu2T2, T1, T2, m12, m21"
        upper = [200, 200, 200, 200, 15, 15, 10, 10]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_sizechange_ancient_asym_mig
    elif model == "split_sizechange_sec_asym_mig":
        num = 8
        p_labels = "nu1T1, nu2T2, nu1T1, nu2T2, T1, T2, m12, m21"
        upper = [200, 200, 200, 200, 15, 15, 10, 10]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_sizechange_second_asym_mig
    elif model == "het_asym_mig":
        num = 8
        p_labels = "nu1, nu2, m12, m21, me12, me21, T, P"
        upper = [200, 200, 10, 10, 10, 10, 15, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.hetero_asym_migration
    elif model == "anc_het_asym_mig":
        num = 9
        p_labels = "nu1, nu2, m12, m21, me12, me21, T1, T2, P"
        upper = [200, 200, 10, 10, 10, 10, 15, 15, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_hetero_asym_migration
    elif model == "sec_het_asym_mig":
        num = 9
        p_labels = "nu1, nu2, m12, m21, me12, me21, T1, T2, P"
        upper = [200, 200, 10, 10, 10, 10, 15, 15, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_hetero_asym_migration
    elif model == "split_bottle_het_asym_mig":
        num = 15
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2, m12T1, m21T1, me12T1, me21T1, m12T2, m21T2, me12T2, me21T2, P"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
                 0.001]
        model_fun = demo_models_kp.split_bottlegrowth_hetero_asym_mig
    elif model == "split_bottle_anc_het_asym_mig":
        num = 11
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2, m12, m21, me12, me21, P"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_bottlegrowth_ancient_hetero_asym_mig
    elif model == "split_bottle_sec_het_asym_mig":
        num = 11
        p_labels = "nu1, nu2, nu1F, nu2F, T1, T2, m12, m21, me12, me21, P"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_bottlegrowth_second_hetero_asym_mig
    elif model == "split_sizechange_het_asym_mig":
        num = 15
        p_labels = "nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12T1, m21T1, me12T1, me21T1, m12T2, m21T2, me12T2, me21T2, P"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
                 0.001]
        model_fun = demo_models_kp.split_sizechange_hetero_asym_mig
    elif model == "split_sizechange_anc_het_asym_mig":
        num = 11
        p_labels = "nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12, m21, me12, me21, P"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_sizechange_ancient_hetero_asym_mig
    elif model == "split_sizechange_sec_het_asym_mig":
        num = 11
        p_labels = "nu1T1, nu2T1, nu1T2, nu2T2, T1, T2, m12, m21, me12, me21, P"
        upper = [200, 200, 200, 200, 15, 15, 10, 10, 10, 10, 1]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.split_sizechange_second_hetero_asym_mig
    elif model == "twoepoch_het_mig":
        num = 17
        p_labels = "nu_1, nu_2, t1, nu11, nu12, m1_12, me1_12, m1_21, me1_21, t2, nu21, nu22, m2_12, m2_21, me2_12, me2_21, P"
        model_fun = demo_models_kp.twoepoch_het_mig
        upper = [200, 200, 10, 200, 200, 10, 10, 10, 10, 10, 200, 200, 10, 10, 10, 10, 1]
        lower = [1e-3, 1e-3, 0, 1e-3, 1e-3, 0, 0, 0, 0, 0, 1e-3, 1e-3, 0, 0, 0, 0, 0]
    elif model == "model_split1_then_23":
        num = 5
        p_labels = "nu1, nu2, nu3, T1_23, T23"
        model_fun = demo_models_kp.model_split1_then_23
        upper = [100, 100, 100, 10, 10]
        lower = [0.001, 0.001, 0.001, 0, 0]
    elif model == "model_split2_then_13":
        num = 5
        p_labels = "nu1, nu2, nu3, T1_23, T23"
        model_fun = demo_models_kp.model_split2_then_13
        upper = [100, 100, 100, 10, 10]
        lower = [0.001, 0.001, 0.001, 0, 0]
    elif model == "model_split3_then_12":
        num = 5
        p_labels = "nu1, nu2, nu3, T1_23, T23"
        model_fun = demo_models_kp.model_split3_then_12
        upper = [100, 100, 100, 10, 10]
        lower = [0.001, 0.001, 0.001, 0, 0]
    else:
        raise ValueError(
            'Model nickname undefined please check in SETTINGS.py you are using the correct model nickname!')
    if ALL:
        return model_fun, num, p_labels, upper, lower
    else:
        return model_fun
