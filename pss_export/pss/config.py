'''
Class to hold config for PSS and the exports.
'''

from pathlib import Path

import yaml
import re

#-------------------------------------
#  Settings class
#-------------------------------------

class Config:
    """
    A class to hold config for PSS and the exports.
    """

    def __init__(self, filename=None, **kwargs):
        """
        Initializes the settings with yaml file and/or provided keyword arguments.

        Parameters
        ----------
        filename (str, optional):
            The path to the YAML file to load settings from.

        **kwargs:
            Arbitrary keyword arguments to set as attributes.
            (Overrides attributes in YAML file if provided.)
        """
        settings = self.load_settings(filename=filename, **kwargs)
        for key, value in settings.items():
            # print(key, value)
            setattr(self, key, value)

    def load_settings(self, filename=None, **kwargs):
        """
        Load settings from a YAML file and/or keyword arguments.
        Parameters
        ----------
        filename (str, optional):
            The path to the YAML file to load settings from.
            If not provided, defaults to None.
        **kwargs:
            Arbitrary keyword arguments to set as attributes.
            These will override attributes loaded from the YAML file if provided.
        """

        if filename:
            with open(filename, 'r', encoding="utf-8") as file:
                settings = yaml.safe_load(file) if filename else None
        else:
            settings = {}

        for key, value in kwargs.items():
            settings[key] = value

        return settings

pss_export_config = Config(filename=Path(__file__).with_name("pss_export_config.yaml"))
pss_schema_config = Config(filename=Path(__file__).with_name("pss_schema_config.yaml"))


def species_codes(species):
    ''' The species to collect for as a tuple of codes of pss_schema_config.species, from a code, a list of codes
    or None (no species filter, returned as None). Raises ValueError for an unknown or repeated code. '''
    if species is None:
        return None
    codes = (species,) if isinstance(species, str) else tuple(species)
    if not codes:
        raise ValueError("No species given (None for no species filter)")
    unknown = [code for code in codes if code not in pss_schema_config.species]
    if unknown:
        raise ValueError(f"Unknown species '{', '.join(unknown)}', one of: "
                         f"{', '.join(pss_schema_config.species)} (or None)")
    if len(set(codes)) < len(codes):
        raise ValueError(f"Species given more than once: {', '.join(codes)}")
    return codes
