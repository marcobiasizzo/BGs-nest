# -*- coding: utf-8 -*-
"""
"""

__author__ = 'marco'

from BGs_nest.manager import manager_class

class BGs_class:
    def __init__(self, nest, number_of_neurons, cortical_mode, parameters_file_name, dopa_depl, cortex_type,
                 in_vitro=False, n_spike_generators='n_neurons'):
        # create Basal Ganglia neurons and connections
        self.N = number_of_neurons  # total BGs pop neurons
        self.mng = manager_class(parameters_file_name, dopa_depl)   # network building manager
        # save BGs populations and indexes
        self.BGs_pops, self.BGs_pop_ids = self.create_bgmodel(nest, in_vitro=in_vitro)

        if not in_vitro:
            # create cortical and external inputs and connect to BGs
            assert cortical_mode == 'active' or cortical_mode == 'slow', 'cortical mode should be defined as active or slow'
            self.CTX_pops = self.create_ctxinput(nest, cortical_mode, squared=True, in_spikes=cortex_type, n_spike_generators=n_spike_generators)
            self.EXT_pops = self.create_extinput(nest, cortical_mode)
            self.connect_cortex_BGs(nest, in_spikes=cortex_type)
            self.connect_external_BGs(nest)


    def create_bgmodel(self, nest_, in_vitro=False):
        """ Function to create a BG connected population
            N is the total number of neurons                """
        pop_dims = self.mng.derive_population_dimentions_from_csv(self.N)  # get the population dim according to N

        izhi_model = 'izhik_cond_exp'  # the two types of neuronal models
        aeif_model = 'my_aeif_cond_exp'
        rec = {izhi_model: nest_.GetDefaults(izhi_model)['receptor_types'],  # save receptors values for later
               aeif_model: nest_.GetDefaults(aeif_model)['receptor_types']}

        stc_syn = 'static_synapse'  # the two types of synapses
        tsod_syn = 'tsodyks_synapse'

        FSN_pop = self.mng.create_population(nest_, izhi_model, 'FS', pop_dims)
        MSND1_pop = self.mng.create_population(nest_, izhi_model, 'M1', pop_dims)
        MSND2_pop = self.mng.create_population(nest_, izhi_model, 'M2', pop_dims)
        if not in_vitro:
            self.mng.create_connection(nest_, tsod_syn, FSN_pop, FSN_pop, 'FS_FS_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_1'])
            self.mng.create_connection(nest_, tsod_syn, FSN_pop, MSND1_pop, 'FS_M1_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_1'])
            self.mng.create_connection(nest_, tsod_syn, FSN_pop, MSND2_pop, 'FS_M2_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_1'])
            self.mng.create_connection(nest_, stc_syn, MSND1_pop, MSND1_pop, 'M1_M1_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_2'])
            self.mng.create_connection(nest_, stc_syn, MSND1_pop, MSND2_pop, 'M1_M2_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_2'])
            self.mng.create_connection(nest_, stc_syn, MSND2_pop, MSND1_pop, 'M2_M1_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_2'])
            self.mng.create_connection(nest_, stc_syn, MSND2_pop, MSND2_pop, 'M2_M2_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_2'])

        GPeTA_pop = self.mng.create_population(nest_, aeif_model, 'GA', pop_dims)
        GPeTI_pop = self.mng.create_population(nest_, aeif_model, 'GI', pop_dims)
        if not in_vitro:
            self.mng.create_connection(nest_, tsod_syn, GPeTA_pop, FSN_pop, 'GA_FS_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_2'])
            self.mng.create_connection(nest_, tsod_syn, GPeTI_pop, FSN_pop, 'GI_FS_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_3'])
            self.mng.create_connection(nest_, stc_syn, GPeTA_pop, MSND1_pop, 'GA_M1_gaba',
                                       receptor_true=rec[izhi_model]['GABAA_3'])
            self.mng.create_connection(nest_, stc_syn, GPeTA_pop, MSND2_pop, 'GA_M2_gaba',  # check, it was 5
                                       receptor_true=rec[izhi_model]['GABAA_3'])
            self.mng.create_connection(nest_, stc_syn, GPeTA_pop, GPeTA_pop, 'GA_GA_gaba',
                                       receptor_true=rec[aeif_model]['GABAA_2'])
            self.mng.create_connection(nest_, stc_syn, GPeTA_pop, GPeTI_pop, 'GA_GI_gaba',
                                       receptor_true=rec[aeif_model]['GABAA_2'])
            self.mng.create_connection(nest_, stc_syn, GPeTI_pop, GPeTA_pop, 'GI_GA_gaba',
                                       receptor_true=rec[aeif_model]['GABAA_2'])
            self.mng.create_connection(nest_, stc_syn, GPeTI_pop, GPeTI_pop, 'GI_GI_gaba',
                                       receptor_true=rec[aeif_model]['GABAA_2'])
            self.mng.create_connection(nest_, tsod_syn, MSND2_pop, GPeTI_pop, 'M2_GI_gaba',
                                       receptor_true=rec[aeif_model]['GABAA_1'])

        STN_pop = self.mng.create_population(nest_, aeif_model, 'ST', pop_dims)
        if not in_vitro:
            self.mng.create_connection(nest_, stc_syn, GPeTI_pop, STN_pop, 'GI_ST_gaba',
                                       receptor_true=rec[aeif_model]['GABAA_1'])
            self.mng.create_connection(nest_, stc_syn, STN_pop, GPeTA_pop, 'ST_GA_ampa',
                                       receptor_true=rec[aeif_model]['AMPA_1'])
            self.mng.create_connection(nest_, stc_syn, STN_pop, GPeTI_pop, 'ST_GI_ampa',
                                       receptor_true=rec[aeif_model]['AMPA_1'])

        SNr_pop = self.mng.create_population(nest_, aeif_model, 'SN', pop_dims)
        if not in_vitro:
            self.mng.create_connection(nest_, tsod_syn, MSND1_pop, SNr_pop, 'M1_SN_gaba',  # direct pathway
                                       receptor_true=rec[aeif_model]['GABAA_1'])
            self.mng.create_connection(nest_, tsod_syn, GPeTI_pop, SNr_pop, 'GI_SN_gaba',  # indirect pathway
                                       receptor_true=rec[aeif_model]['GABAA_2'])
            self.mng.create_connection(nest_, tsod_syn, STN_pop, SNr_pop, 'ST_SN_ampa',  # hyperdirect pathway
                                       receptor_true=rec[aeif_model]['AMPA_1'])

        BGs_pops = {'FSN': FSN_pop, 'MSND1': MSND1_pop, 'MSND2': MSND2_pop, 'GPeTA': GPeTA_pop,
                    'GPeTI': GPeTI_pop, 'STN': STN_pop, 'SNr': SNr_pop}
        save_ids = lambda pop: (min(pop), max(pop))
        pop_ids = {'FSN': save_ids(FSN_pop), 'MSND1': save_ids(MSND1_pop), 'MSND2': save_ids(MSND2_pop),
                   'GPeTA': save_ids(GPeTA_pop), 'GPeTI': save_ids(GPeTI_pop), 'STN': save_ids(STN_pop),
                   'SNr': save_ids(SNr_pop)}
        return BGs_pops, pop_ids

    def create_ctxinput(self, nest_, cortical_mode_, squared, in_spikes='poisson', n_spike_generators='n_neurons'):
        rate_FSN = {'active': 787.0, 'slow': 646.0}
        a_FSN = {'active': 0.11, 'slow': 0.11}
        rate_MSN1 = {'active': 546.0, 'slow': 448.0}
        a_MSN1 = {'active': 0.11, 'slow': 0.11}
        rate_MSN2 = {'active': 722.0, 'slow': 592.0}
        a_MSN2 = {'active': 0.11, 'slow': 0.11}
        rate_STN = {'active': 250.0, 'slow': 170.0}
        a_STN = {'active': 0.35, 'slow': 0.9}
        oscillations_freq = {'active': 20., 'slow': 1.}
        if squared:
            poisson_type = 'poisson_generator_periodic'  # generates a squared poisson input
            poisson_periodic_params = lambda rate_, a_: {
                'period_first': (1. / oscillations_freq[cortical_mode_] * 1000.) / 2.,
                'period_second': (1. / oscillations_freq[cortical_mode_] * 1000.) / 2.,
                'rate_first': rate_[cortical_mode_] * (1. + a_[cortical_mode_]),
                'rate_second': rate_[cortical_mode_] * (1. - a_[cortical_mode_])}
        else:   # then, it is sinusoidal
            poisson_type = 'sinusoidal_poisson_generator'  # generates a poisson input modulated by a sinusoid
            poisson_periodic_params = lambda rate_, a_: {
                'rate': rate_[cortical_mode_],
                'amplitude': rate_[cortical_mode_] * a_[cortical_mode_],
                'frequency': oscillations_freq[cortical_mode_]}

        # input is poisson generator with squared or sinusoidal input
        if in_spikes == 'poisson':
            CF_pop = nest_.Create(poisson_type, params=poisson_periodic_params(rate_FSN, a_FSN))
            # 1 will transmit ampa, 2 will transmit nmda
            CM1_pop1 = nest_.Create(poisson_type, params=poisson_periodic_params(rate_MSN1, a_MSN1))
            CM1_pop2 = nest_.Create(poisson_type, params=poisson_periodic_params(rate_MSN1, a_MSN1))
            CM2_pop1 = nest_.Create(poisson_type, params=poisson_periodic_params(rate_MSN2, a_MSN2))
            CM2_pop2 = nest_.Create(poisson_type, params=poisson_periodic_params(rate_MSN2, a_MSN2))
            CS_pop1 = nest_.Create(poisson_type, params=poisson_periodic_params(rate_STN, a_STN))
            CS_pop2 = nest_.Create(poisson_type, params=poisson_periodic_params(rate_STN, a_STN))

        # input is a spike generator, it will be summed to a poisson generator, see connect_cortex_BGs
        # useful for the multiarea model with BGs and Cerebellum
        elif in_spikes == 'spike_generator':
            # one spike generator for each neuron
            if n_spike_generators == 'n_neurons':
                pop_dims = self.mng.derive_population_dimentions_from_csv(self.N)  # get the population dim according to N
                n_s_g = pop_dims
            # user-defined number of spike generator
            else:
                n_s_g = n_spike_generators

            CF_pop = nest_.Create("spike_generator", n_s_g['FS'])
            CM1_pop1 = nest_.Create("spike_generator", n_s_g['M1'])
            CM1_pop2 = nest_.Create("spike_generator", n_s_g['M1'])
            CM2_pop1 = nest_.Create("spike_generator", n_s_g['M2'])
            CM2_pop2 = nest_.Create("spike_generator", n_s_g['M2'])
            CS_pop1 = nest_.Create("spike_generator", n_s_g['ST'])
            CS_pop2 = nest_.Create("spike_generator", n_s_g['ST'])

        # classical poisson generator with constant (and updatable) rate
        elif in_spikes == 'dynamic_poisson':
            CF_pop = nest_.Create("poisson_generator", params={'rate': rate_FSN[cortical_mode_]})
            CM1_pop1 = nest_.Create("poisson_generator", params={'rate': rate_MSN1[cortical_mode_]})
            CM1_pop2 = nest_.Create("poisson_generator", params={'rate': rate_MSN1[cortical_mode_]})
            CM2_pop1 = nest_.Create("poisson_generator", params={'rate': rate_MSN2[cortical_mode_]})
            CM2_pop2 = nest_.Create("poisson_generator", params={'rate': rate_MSN2[cortical_mode_]})
            CS_pop1 = nest_.Create("poisson_generator", params={'rate': rate_STN[cortical_mode_]})
            CS_pop2 = nest_.Create("poisson_generator", params={'rate': rate_STN[cortical_mode_]})

        CTX_pops = {'CF': CF_pop, 'CM1_1': CM1_pop1, 'CM1_2': CM1_pop2, 'CM2_1': CM2_pop1, 'CM2_2': CM2_pop2,
                    'CS_1': CS_pop1, 'CS_2': CS_pop2}
        return CTX_pops

    def create_extinput(self, nest_, cortical_mode_):
        rate_EA = {'active': 200.0, 'slow': 200.0}
        rate_EI = {'active': 1530.0, 'slow': 620.0}
        rate_ES = {'active': 1800.0, 'slow': 2000.0}

        EA_pop = nest_.Create('poisson_generator', params={'rate': rate_EA[cortical_mode_]})
        EI_pop = nest_.Create('poisson_generator', params={'rate': rate_EI[cortical_mode_]})
        ES_pop = nest_.Create('poisson_generator', params={'rate': rate_ES[cortical_mode_]})

        EXT_pops = {'EA': EA_pop, 'EI': EI_pop, 'ES': ES_pop}
        return EXT_pops


    def connect_cortex_BGs(self, nest_, in_spikes='poisson'):
        CF_pop = self.CTX_pops['CF']
        CM1_pop1 = self.CTX_pops['CM1_1']
        CM1_pop2 = self.CTX_pops['CM1_2']
        CM2_pop1 = self.CTX_pops['CM2_1']
        CM2_pop2 = self.CTX_pops['CM2_2']
        CS_pop1 = self.CTX_pops['CS_1']
        CS_pop2 = self.CTX_pops['CS_2']

        FSN_pop = self.BGs_pops['FSN']
        MSND1_pop = self.BGs_pops['MSND1']
        MSND2_pop = self.BGs_pops['MSND2']
        STN_pop = self.BGs_pops['STN']

        izhi_model = 'izhik_cond_exp'
        aeif_model = 'my_aeif_cond_exp'
        stc_syn = 'static_synapse'
        rec = {izhi_model: nest_.GetDefaults(izhi_model)['receptor_types'],  # save receptors values for later
               aeif_model: nest_.GetDefaults(aeif_model)['receptor_types']}

        self.mng.create_connection(nest_, stc_syn, CF_pop, FSN_pop, 'CF_FS_ampa',
                                   receptor_true=rec[izhi_model]['AMPA_1'], input_type=in_spikes)
        self.mng.create_connection(nest_, stc_syn, CM1_pop1, MSND1_pop, 'C1_M1_ampa',
                                   receptor_true=rec[izhi_model]['AMPA_1'], input_type=in_spikes)
        self.mng.create_connection(nest_, stc_syn, CM1_pop2, MSND1_pop, 'C1_M1_nmda',
                                   receptor_true=rec[izhi_model]['NMDA_1'], input_type=in_spikes)
        self.mng.create_connection(nest_, stc_syn, CM2_pop1, MSND2_pop, 'C2_M2_ampa',
                                   receptor_true=rec[izhi_model]['AMPA_1'], input_type=in_spikes)
        self.mng.create_connection(nest_, stc_syn, CM2_pop2, MSND2_pop, 'C2_M2_nmda',
                                   receptor_true=rec[izhi_model]['NMDA_1'], input_type=in_spikes)
        self.mng.create_connection(nest_, stc_syn, CS_pop1, STN_pop, 'CS_ST_ampa',
                                   receptor_true=rec[aeif_model]['AMPA_1'], input_type=in_spikes)
        self.mng.create_connection(nest_, stc_syn, CS_pop2, STN_pop, 'CS_ST_nmda',
                                   receptor_true=rec[aeif_model]['NMDA_1'], input_type=in_spikes)


        if in_spikes == 'spike_generator':
            # create also a background activity
            CF_pop_pois = nest_.Create("poisson_generator", params={'rate': 700.4})
            CM1_pop1_pois = nest_.Create("poisson_generator", params={'rate': 486.})
            CM1_pop2_pois = nest_.Create("poisson_generator", params={'rate': 486.})
            CM2_pop1_pois = nest_.Create("poisson_generator", params={'rate': 642.6})
            CM2_pop2_pois = nest_.Create("poisson_generator", params={'rate': 642.6})
            CS_pop1_pois = nest_.Create("poisson_generator", params={'rate': 162.5})
            CS_pop2_pois = nest_.Create("poisson_generator", params={'rate': 162.5})

            self.mng.create_connection(nest_, stc_syn, CF_pop_pois, FSN_pop, 'CF_FS_ampa',
                                       receptor_true=rec[izhi_model]['AMPA_1'], input_type='second_poisson')
            self.mng.create_connection(nest_, stc_syn, CM1_pop1_pois, MSND1_pop, 'C1_M1_ampa',
                                       receptor_true=rec[izhi_model]['AMPA_1'], input_type='second_poisson')
            self.mng.create_connection(nest_, stc_syn, CM1_pop2_pois, MSND1_pop, 'C1_M1_nmda',
                                       receptor_true=rec[izhi_model]['NMDA_1'], input_type='second_poisson')
            self.mng.create_connection(nest_, stc_syn, CM2_pop1_pois, MSND2_pop, 'C2_M2_ampa',
                                       receptor_true=rec[izhi_model]['AMPA_1'], input_type='second_poisson')
            self.mng.create_connection(nest_, stc_syn, CM2_pop2_pois, MSND2_pop, 'C2_M2_nmda',
                                       receptor_true=rec[izhi_model]['NMDA_1'], input_type='second_poisson')
            self.mng.create_connection(nest_, stc_syn, CS_pop1_pois, STN_pop, 'CS_ST_ampa',
                                       receptor_true=rec[aeif_model]['AMPA_1'], input_type='second_poisson')
            self.mng.create_connection(nest_, stc_syn, CS_pop2_pois, STN_pop, 'CS_ST_nmda',
                                       receptor_true=rec[aeif_model]['NMDA_1'], input_type='second_poisson')

    def connect_external_BGs(self, nest_):
        EA_pop = self.EXT_pops['EA']
        EI_pop = self.EXT_pops['EI']
        ES_pop = self.EXT_pops['ES']

        GPeTA_pop = self.BGs_pops['GPeTA']
        GPeTI_pop = self.BGs_pops['GPeTI']
        SNr_pop = self.BGs_pops['SNr']

        izhi_model = 'izhik_cond_exp'
        aeif_model = 'my_aeif_cond_exp'
        stc_syn = 'static_synapse'
        rec = {izhi_model: nest_.GetDefaults(izhi_model)['receptor_types'],  # save receptors values for later
               aeif_model: nest_.GetDefaults(aeif_model)['receptor_types']}

        self.mng.create_connection(nest_, stc_syn, EA_pop, GPeTA_pop, 'EA_GA_ampa',
                                   receptor_true=rec[aeif_model]['AMPA_2'], input_type='poisson')
        self.mng.create_connection(nest_, stc_syn, EI_pop, GPeTI_pop, 'EI_GI_ampa',
                                   receptor_true=rec[aeif_model]['AMPA_2'], input_type='poisson')
        self.mng.create_connection(nest_, stc_syn, ES_pop, SNr_pop, 'ES_SN_ampa',
                                   receptor_true=rec[aeif_model]['AMPA_2'], input_type='poisson')
