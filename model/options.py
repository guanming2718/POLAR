import argparse
import configparser
import sys
import os.path


def read_parameters():
    parser = argparse.ArgumentParser(description='Active Polar Dynamics')
    parser.add_argument('-in',help='The location of the input file')
    parser.add_argument('-out',help='The location of the output folder')
    parser.add_argument('-v',help='1D or 2D version')
    args = vars(parser.parse_args())
    input_file = args['in']
    output_dir = args['out']
    # check if the input and output are specified
    if input_file is None:
        raise Exception('The input file is None,use -in to specify it')
    if output_dir is None:
        raise Exception('The output folder is None,use -out to specify it')
    if args['v'] is None:
        raise Exception('The version is None,use -v to specify it')
    #check if the input file exists
    if os.path.isfile(input_file):
        config = configparser.ConfigParser()
        config.read(input_file)
    else:
        raise Exception('The input file does not exist')
    settings = dict()
    settings['input'] = input_file
    settings['output'] = output_dir
    settings['version'] = args['v']
    
    config_params = dict()
    config_params['nx'] = config.getint('Configuration','nx')
    config_params['ny'] = config.getint('Configuration','ny')
    config_params['Lx'] = config.getint('Configuration','Lx')
    config_params['Ly'] = config.getint('Configuration','Ly')
    config_params['init-config'] = config['Configuration']['init-config']
    config_params['BC_density'] = config['Configuration']['BC_density']
    config_params['BC_polar'] = config['Configuration']['BC_polar']
    config_params['BC_velocity'] = config['Configuration']['BC_velocity']
    if (config_params['init-config'] =='cluster') or (config_params['init-config'] == 'circular-wound'):
            config_params['R'] = config.getfloat('Configuration','R')
    elif (config_params['init-config'] == 'two-side-wound') or (config_params['init-config'] == 'one-side-wound'):
            config_params['wound-ratio'] = config.getfloat('Configuration','wound-ratio')
    # model parameters
    model_params = dict()
    parameters = ['Ac','Bc','Kc','Ap','Bp','Kp','w','gamma','nu','zeta',
                       'alpha','Ks','Kv','xi']
    for p in parameters:
        model_params[p] = config.getfloat('Model',p)
    return [settings,config_params,model_params]

def print_settings(settings,config_params,model_params):
    print(" ╔═╗┌─┐┌┐┌┌┬┐┬┌┐┌┬ ┬┬ ┬┌┬┐  ╔╦╗┌─┐┌┬┐┌─┐┬    ┌─┐┌─┐┬─┐  ┌┬┐┬┌─┐┌─┐┬ ┬┌─┐")
    print(" ║  │ ││││ │ │││││ ││ ││││  ║║║│ │ ││├┤ │    ├┤ │ │├┬┘   │ │└─┐└─┐│ │├┤ ")
    print(" ╚═╝└─┘┘└┘ ┴ ┴┘└┘└─┘└─┘┴ ┴  ╩ ╩└─┘─┴┘└─┘┴─┘  └  └─┘┴└─   ┴ ┴└─┘└─┘└─┘└─┘\n")
    print("-"*65)
    print("Settings")
    for s in settings:
        print("|{:<30}|{:>20}".format(s,settings[s]))
    print("-"*65)
    print("Confiuration")
    for p in config_params:
        print("|{:<30}|{:>20}".format(p,config_params[p]))
    print("-"*65)
    print("Model")
    for p in model_params:
        print("|{:<30}|{:>20}".format(p,model_params[p]))
    print("-"*65)


    




