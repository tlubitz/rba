# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
# package imports
from rba.core.constraint_blocks import ConstraintBlocks
from rbastructure.element_block import ElementBlock


class ModuleBlock(ElementBlock):
    """
    Class holding information on the modules in the model.

  Attributes
  ----------
  Elements : dict
      Each model-enzyme is represented by a key.
      The values, holding information on each enzyme, are dicts with predefined keys:
          'ID' : So far dummy (type str)
          'Name' : So far dummy (type str)
    """

    def fromFiles(self, model, Info):
        self.Elements = {}
        self.Elements['DummyModule'] = {'ID': 'Dummy',
                                        'Name': 'Dummy'}

    def overview(self):
        """
        Derive statistics on modules.

        Returns
        -------
        Dictionary with general numbers on modules.

        """
        nT = len(self.Elements.keys())
        out = {'nModulesTotal': nT}
        return(out)
