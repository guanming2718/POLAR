import argparse
import configparser
import sys
import os.path
sys.path.append('./model/*')
from model import options
from model import polar2D

[settings,config_params,model_params] = options.read_parameters()
options.print_settings(settings,config_params,model_params)


    




