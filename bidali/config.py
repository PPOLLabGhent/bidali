# -*- coding: utf-8 -*-
"""Bidali configuration

This module organizes the configuration for the general bidali package and
the LSD subpackage.

Todo:
    * Give example symlinking privatedir ln -s ~/Dropbiz/Lab/z_archive/Datasets ~/LSData/private
"""

import configparser, os
configFileOptions = [
    'bidali.cfg', # in current working dir
    os.path.expanduser('~/.bidali.cfg'),
    '/usr/local/etc/bidali.cfg'
]

# Default configuration
config = configparser.ConfigParser()
config['bidali'] = {
    'plotting_backend': 'TkAgg'
}

config['LSD'] = {
    'cachetime': '4w', #supports w[eeks], d[ays] or h[ours]
    'cachedir': os.path.expanduser('~/LSData/cache/'),
    'privatedir': os.path.expanduser('~/LSData/private/')
}

config['plotting'] = {
    'mpl_backend': 'TkAgg',
    'interactive': 'yes',
    'use_seaborn_bug_workaround': 'no'
}

# Read configuration file
for configFile in configFileOptions:
    if os.path.exists(configFile):
        config.read(configFile)
        break #only reads the first config file found
