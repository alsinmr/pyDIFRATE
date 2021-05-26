# __init__.py


from pyDIFRATE.r_class import all_sens as sens
#import pyDIFRATE.r_class as sens

from .Struct.structure import molecule
from .data.data_class import data
from .data import fitting
from .data import in_out as io
from .plots import plotting_funs as plotting
from .tools import DRtools as tools
from .iRED import Ct_fast
from .Struct import eval_fr as frames
from .chimera import chimeraX_funs as chimeraX
from .chimera import cmx_3D_plots as cmx_plots

#import Sensitivities
#"Temporary hack to recover old data objects that have been saved"
#import os
#os.chdir('pyDIFRATE')
#import DR_old
#os.chdir('..')



