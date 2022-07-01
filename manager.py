import csv
from numpy import random

class manager_class:
    def __init__(self, parameters_file_name, dopa_depl_level):
        self.params_file_name = parameters_file_name
        self.dopa_depl = dopa_depl_level
        assert self.dopa_depl >= -0.8 and self.dopa_depl <=0., f'The level of dopamine depletion should be between 0 (healthy) and -0.8 (severe Parkinson)'

    def get_csv_value(self, p0, p1, p2, convert_to_type):
        """ Function to get a single value from csv file """
        with open(self.params_file_name, mode='r') as csv_file:
            value = None
            csv_reader = csv.DictReader(csv_file, delimiter=',')
            for idx, row in enumerate(csv_reader):  # read all csv rows
                if row['p0'] == p0:
                    if row['p1'] == p1:
                        if row['p2'] == p2:
                            try:
                                value = convert_to_type(row['p3'])
                            except ValueError:
                                print(f'Could not convert to {convert_to_type}: {row["p3"]}')
        return value

    def copy_params_from_csv(self, area, element, model):
        """ Function to copy parameters of  from csv file """
        params = {}
        with open(self.params_file_name, mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file, delimiter=',')
            for idx, row in enumerate(csv_reader):  # read all csv rows
                if row['p0'] == area:
                    if row['p1'] == element:
                        if row['p2'] != 'template':
                            params[row['p2']] = float(
                                    row['p3'])  # if possible, convert element to float and save in dic
                        # except ValueError:  # if it is a literal string
                        #    print(f'Could not convert to float: {row["p3"]}')
                        elif row["p2"] == 'template':
                            assert row["p3"] == model, 'Model different from csv declaration'
        return params

    def derive_population_dimentions_from_csv(self, N):
        """ Function to derive pop dimentions from csv file and save to dictionary """
        pop_dims = {}
        with open(self.params_file_name, mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file, delimiter=',')
            for idx, row in enumerate(csv_reader):  # read all csv rows
                if row['p2'] == 'prop':
                    try:
                        prop = float(row['p3'])  # if possible, convert element to float and save in dic
                        pop_dims[row['p1']] = round(prop * N)
                    except ValueError:  # if it is a literal string
                        print(f'Could not convert to float: {row["p3"]}')
        assert abs(sum(pop_dims.values()) - N) <= 1, f'total number of neuron is not N: {sum(pop_dims.values())}'
        return pop_dims

    def create_population(self, nest, model, population, pop_dims, variations=True):
        """ Function to create a population and set the proper params """
        params = self.copy_params_from_csv('nest', population, model)
        params['tata_dop'] = self.dopa_depl  # save dopa depl level
        pop = nest.Create(model, pop_dims[population], params)
        if variations:
            dCms = [{"C_m": random.normal(params['C_m'], params['C_m'] * 0.1)} for _ in range(pop_dims[population])]
            nest.SetStatus(pop, dCms)
            dVts = [{"V_th": random.normal(params['V_th'], 1)} for _ in range(pop_dims[population])]
            nest.SetStatus(pop, dVts)
            dVms = [{"V_m": random.uniform(params['V_m'] - 10, params['V_m'] + 10)} for _ in range(pop_dims[population])]
            nest.SetStatus(pop, dVms)

        return pop


    def create_connection(self, nest, template, pre_pop, post_pop, connection, receptor_true, input_type='neurons'):
        """ Function to connect 2 pops and set the proper params """

        params = self.copy_params_from_csv('nest', connection, template)

        receptor_type = int(params.pop('receptor_type'))        # get and remove item from dic
        if receptor_true is not None:
            assert receptor_type == receptor_true, 'Receptor number different from csv declaration'

        if input_type != 'second_poisson':
            nest.CopyModel(template, f'{connection}_syn', params)
            # connect depending on the input population/generator
            if input_type == 'poisson' or input_type == 'dynamic_poisson':  # one poisson generator connected to all of them
                conn = nest.Connect(pre_pop, post_pop, {'rule': 'all_to_all', "multapses": False, "autapses": False},
                                    syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type})
            elif input_type == 'spike_generator':
                assert len(pre_pop) <= len(post_pop), 'Input generators should be less or equal to target neurons'
                if len(pre_pop) == len(post_pop):   # if we have 1 spike generator for each target neuron
                    conn = nest.Connect(pre_pop, post_pop, {'rule': 'one_to_one', "multapses": False, "autapses": False},
                                        syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type})
                else:   # some neurons will receive the same input
                    post_pop = list(post_pop)
                    random.shuffle(post_pop)
                    n_s_g = len(pre_pop)
                    n_targets = len(post_pop)/n_s_g
                    for i in range(n_s_g - 1):
                        post = post_pop[round(i*n_targets):round((i+1)*n_targets)]
                        # connect one source neuron to a group of target neurons
                        conn = nest.Connect([pre_pop[i]], post,
                                            {'rule': 'all_to_all', "multapses": False, "autapses": False},
                                            syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type})
                    post = post_pop[round((n_s_g - 1)*n_targets):]
                    conn = nest.Connect([pre_pop[n_s_g - 1]], post,
                                        {'rule': 'all_to_all', "multapses": False, "autapses": False},
                                        syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type})

            elif input_type == 'parrots':  # for multiple external inputs use parrots neurons
                conn = nest.Connect(pre_pop, post_pop, {'rule': 'one_to_one', "multapses": False, "autapses": False},
                                    syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type})

            elif input_type == 'neurons':
                in_deg = round(self.get_csv_value('conn', connection, 'fan_in', int))   # read the indegree from csv
                weight_dic = {'distribution': 'uniform', 'low': params['weight'] - params['weight'] * 0.5,
                                                         'high': params['weight'] + params['weight'] * 0.5}
                if connection == 'GI_FS_gaba':
                    # only 10% of GPeTI connects to FS
                    pre = pre_pop[0:int(0.1 * (max(pre_pop) - min(pre_pop)))]  # select the 10% of GI, to connect to FSN
                    conn = nest.Connect(pre, post_pop,
                                        {'rule': 'fixed_indegree', 'indegree': in_deg, "multapses": False, "autapses": False},
                                        syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type,
                                                  'weight': weight_dic})
                elif connection == 'FS_M1_gaba' or connection == 'FS_M2_gaba':
                    # groups of FS could connect to just 540 MSN
                    conn = restrict_connect_MSN(540, pre_pop, post_pop, nest, in_deg, connection,
                                                receptor_type, weight_dic)
                elif connection == 'M1_M1_gaba' or connection == 'M1_M2_gaba' or \
                        connection == 'M2_M1_gaba' or connection == 'M2_M2_gaba':
                    # groups of MSN could connect to just 2800 MSN
                    conn = restrict_connect_MSN(2800, pre_pop, post_pop, nest, in_deg, connection,
                                                receptor_type, weight_dic)
                else:   # other connections are free
                    conn = nest.Connect(pre_pop, post_pop,
                                        {'rule': 'fixed_indegree', 'indegree': in_deg, "multapses": False, "autapses": False},
                                        syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type,
                                                  'weight': weight_dic})
            else:
                print('WARNING: input_type connection not found')
        elif input_type == 'second_poisson':
            nest.CopyModel(template, f'{connection}_syn_bground', params)
            conn = nest.Connect(pre_pop, post_pop, {'rule': 'all_to_all', "multapses": False, "autapses": False},
                                syn_spec={'model': f'{connection}_syn', 'receptor_type': receptor_type})
        return conn


def restrict_connect_MSN(post_blocks, pre_pop_, post_pop_, nest_, in_deg_, connection_, receptor_type_, weight_dic_):
    """ Create a spatial connectivity in FSN/MSN """

    N_in = max(pre_pop_) - min(pre_pop_)        # input pop dim
    N_out = max(post_pop_) - min(post_pop_)     # output pop dim
    ratio_o_i = int(N_out/N_in)                 # output/input ratio
    in_rec_field = int((post_blocks-1) / ratio_o_i)     # input cells connecting to one post block
    for i_out in range(N_out):
        i_in = int(i_out / ratio_o_i)
        i_in1 = i_in - int(in_rec_field/2)      # inferior range of inputs
        i_in2 = i_in1 + in_rec_field            # superior range of inputs
        if i_in1 < 0:
            pre = pre_pop_[i_in1+N_in:-1] + pre_pop_[0:i_in2]
        elif i_in2 >= N_in:
            pre = pre_pop_[i_in1:-1] + pre_pop_[0:i_in2-N_in]
        else:
            pre = pre_pop_[i_in1:i_in2]
        post = [post_pop_[i_out]]
        conn = nest_.Connect(pre, post,
                            {'rule': 'fixed_indegree', 'indegree': in_deg_, "multapses": False, "autapses": False},
                            syn_spec={'model': f'{connection_}_syn', 'receptor_type': receptor_type_,
                                      'weight': weight_dic_})
    return conn