def matching_ms_hb_choice_q_model(data_holder_path):
    """
    ms = multi_session
    hb = hierarchical bayesian
        < - find hyperparameters using simulation
            & MCMC(Hamiltonian MC) sampling
            from hyperparameters to get prior
    Matching Analysis Function
    Developed by Antonio Lee @ BCCN Berlin
    Version 1.0 ~ May 2026
    Model iteration see the end of script
    """

    # import libraries
    import os
    import scipy.io as sio
    import math
    import numpy as np
    import matlab.engine
    from cmdstanpy import CmdStanModel

    # leverage existing matlab scripts
    eng = matlab.engine.start_matlab()

    # import data_holder
    data_folder_path = eng.OttLabDataServerFolderPath()
    # data_holder_path = '71\\bpod_session_pooled_analysis\\71_TwoArmBanditVariant_20241028_20250321\\Selected_Data'
    full_data_holder_path = os.path.join(data_folder_path, data_holder_path)
    data_holder = sio.loadmat(full_data_holder_path)
    data_holder = data_holder["DataHolder"]

    # convert the data in python dictionary format
    # how missing data is modelled or treated?
    data_code = """
        
        """

    """
        data {
            int
        }
        transformed data {}
        parameters {}
        transformed parameters {}
        model {}
        """

    model_code = """
        data {
            int n_session;
        }
        parameters {
            real logit_mu_alpha;
            real log_kappa_alpha;
            real logit_mu_gamma_v;
            real log_kappa_gamma_v;
            real logit_mu_gamma_m;
            real log_kappa_gamma_m;
            real mu_beta;
            real log_precision_beta;
            real mu_phi;
            real log_precision_phi;
            real mu_bias;
            real log_precision_bias;
            array[n_sessions] real logit_alpha;
            array[n_sessions] real logit_gamma_v;
            array[n_sessions] real logit_gamma_m;
            array[n_sessions] real beta;
            array[n_sessions] real phi;
            array[n_sessions] real bias;
        }
        transformed parameters {
            // convert parameters from real number space to working space
            // hyper-parameters
            real mu_alpha = inv_logit(logit_mu_alpha);
            real kappa_alpha = exp(log_kappa_alpha);
            real alpha_alpha = mu_alpha * kappa_alpha;
            real mu_gamma_v = inv_logit(logit_mu_gamma_v);
            real kappa_gamma_v = exp(log_kappa_gamma_v);
            real alpha_gamma_v = mu_gamma_v * kappa_gamma_v;
            real mu_gamma_m = inv_logit(logit_mu_gamma_m);
            real kappa_gamma_m = exp(log_kappa_gamma_m);
            real alpha_gamma_m = mu_gamma_m * kappa_gamma_m;
            real precision_beta = inv_logit(log_precision_beta);
            real sigma_beta = sqrt(precision_beta);
            real precision_phi = inv_logit(log_precision_phi);
            real sigma_phi = sqrt(precision_phi);
            real precision_bias = inv_logit(log_precision_bias);
            real sigma_bias = sqrt(precision_bias);
            
            // session parameters
            array[n_sessions] real alpha = inv_logit(logit_alpha);
            array[n_sessions] real gamma_v = inv_logit(logit_gamma_v);
            array[n_sessions] real gamma_m = inv_logit(logit_gamma_m);
            
            // 
            array[] choice_left;
        }
        model {
            // pdf of hyper-parameters with hyper-prior
            mu_alpha ~ beta(4, 12);
            kappa_alpha ~ gamma(16, 1);
            mu_gamma_v ~ beta(2, 10);
            kappa_gamma_v ~ gamma(12, 1);
            mu_gamma_m ~ beta(10, 6);
            kappa_gamma_m ~ gamma(16, 1);
            mu_beta ~ normal(8, 1);
            precision_beta ~ gamma(1, 1);
            mu_phi ~ normal(-1, 1);
            precision_phi ~ gamma(1, 1);
            mu_bias ~ normal(0, 1);
            precision_bias ~ gamma(1, 1);
            
            // pdf of parameters with prior
            alpha ~ beta(alpha_alpha, beta_alpha);
            gamma_v ~ beta(alpha_gamma_v, beta_gamma_v);
            gamma_m ~ beta(alpha_gamma_m, beta_gamma_m);
            beta ~ normal(mu_beta, sigma_beta);
            phi ~ normal(mu_phi, sigma_phi);
            bias ~ normal(mu_bias, sigma_bias);
            
            // pdf of session nll
            for i_session = 1:n_sessions{
                choice_left ~ bernoulli_logit(b);
            }      
        }
        """
    # hyper - prior
    # a set initial point of MCMC around simulation results
    sampler_initial_hyper_parameters = [
        0.25, 8,  # LearningRate: Mu, Kappa
        8, 1,     # InverseTemperature: Mean, precision
        1/6, 8,   # ForgettingRate: Mu, Kappa
        -1, 1,    # ChoiceStickiness: Mean, precision
        5/8, 8,   # ChoiceForgettingRate: Mu, Kappa
        0, 1,     # Bias: Mean, precision
        ]

    # convert Parameters to real number space
    for idx in range(len(sampler_initial_hyper_parameters)):
        if idx % 4 == 0:  # matlab idx 1, 5, 9
            sampler_initial_hyper_parameters[idx] = math.log(
                sampler_initial_hyper_parameters[idx] / (1 - sampler_initial_hyper_parameters[idx]))

        if idx % 2 == 1:  # matlab idx 2, 4, 6, 8, 10, 12
            sampler_initial_hyper_parameters[idx] = math.log(sampler_initial_hyper_parameters[idx])

    # nSession x nParameter
    session_initial_parameters = [
        0.25, 8,
        1/6, -1,
        5/8, 0
        ]

    n_sessions = data_holder.size
    session_initial_parameters = repmat(SesseionInitialParameters, [1, nSessions])
    session_initial_parameters = reshape(SamplerInitialParameters, [], 1)



    lin_reg_dat = {
        'n': n,
        'x': x,
        'y': y
        }

    model = CmdStanModel(
        modelname="matching_ms_hb_choice_q_model",
        exe_file="")  # or stan_file
    fit = model.sample(
        data="",
        chains=2,
        parallel_chains=1,
        show_console=True,
        iter_warmup=100000,
        iter_sampling=100000,
        show_progress=True,
        )  # in json format

    """
    for i_session in range(data_holder.size):
        print(i_session)
        session_data = data_holder[0, i_session]
        
        trial_data = session_data[0, 0]["Custom"][0, 0]["TrialData"]
        
    """

    # model = session_data_path
    print("")
    return

if __name__ == "__main__":
    matching_ms_hb_choice_q_model("71\\bpod_session_pooled_analysis\\71_TwoArmBanditVariant_20241028_20250321\\Selected_Data")