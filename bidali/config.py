# Bidali configuration
import configparser, os
configFileOptions = [
    'bidali.cfg', # in current working dir
    os.path.expanduser('~/.bidali.cfg'),
    '/etc/bidali.cfg'
]

# Default configuration
config = configparser.ConfigParser()
config['bidali'] = {
    'plotting_backend': 'TkAgg'
}

config['LSD'] = {
    'cache': '4w', #supports w[eeks], d[ays] or h[ours]
    'cachedir': os.path.expanduser('~/LSData/')
}

# Read configuration file
for configFile in configFileOptions:
    if os.path.exists(configFile):
        config.read(configFile)
        break #only reads the first config file found
