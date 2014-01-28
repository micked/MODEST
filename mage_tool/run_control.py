
"""
Configuration/run control module for mage_tool

This module contains a dictionary containing the module-wide configuration.
"""

import yaml
import os.path


CONF = {#Oligo options
        "oligo_length": 90,
        "min_homology": 15,
        "long_homology": 50,

        #Program options
        "processes": None,

        #Operations
        "operations": ['mage_tool.operations.manual',
                       'mage_tool.operations.translation',
                       'mage_tool.operations.promoter'],

        #Gene options
        "leader_length": 35,
        "promoter_length": 250,
        }


def load_rcfile(filepath):
    """Load a mage_tool configuration.

    Configuration must be in YAML format.

    """
    if os.path.isfile(filepath):
        with open(filepath, "r") as f:
            CONF.update(yaml.safe_load(f))


#DEFAULT_RC_FILES = ["~/.conf/mage_tool.rc", "./mage_tool.rc",]

#for rcfile in DEFAULT_RC_FILES:
    #load_rcfile(rcfile)

