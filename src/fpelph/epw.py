#region modules
import os
import h5py
import numpy as np
from ase.units import Hartree, Bohr, Rydberg, _amu, _me
import xmltodict
import jmespath
from fpflow.inputs.inputyaml import InputYaml
from fpflow.structure.qe.qe_struct import QeStruct
from benedict import benedict
from ase import Atoms
from ase.data import atomic_numbers, atomic_masses
from fpflow.structure.kpts import Kpts

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class Epw:
    def __init__(self, src_folder: str = './save'):
        self.src_folder: str = src_folder
        self.inputdict: dict = InputYaml.from_yaml_file().inputdict
        self.qestruct: QeStruct = QeStruct.from_yaml_file()
        self.atoms: Atoms = self.qestruct.atoms[self.qestruct.struct_idx]

        # Some data from input yaml. 
        self.read_input_yaml()

    def read_input_yaml(self):
        # Read epw_nk for number of pools. 
        self.npool: int = jmespath.search('epw.job_info.nk', self.inputdict)
        self.npool = 1 if self.npool is None else self.npool

        self.nocc: int = self.qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        self.nc: int = jmespath.search('epw.cond_bands', self.inputdict)
        self.nv: int = jmespath.search('epw.val_bands', self.inputdict)
        self.kgrid: int = jmespath.search('dfpt.qgrid', self.inputdict)

    def read_structure_info(self):
        self.prefix: str = 'struct'
        self.nat: int = len(self.atoms)
        self.ntyp: int = len(set(self.atoms.get_chemical_symbols()))
        self.positions: np.ndarray = self.atoms.get_positions() / Bohr
        self.lat: np.ndarray = self.atoms.get_cell().array / Bohr
        self.scaled_positions: np.ndarray = self.atoms.get_scaled_positions()
        self.species: list = list(set(self.atoms.get_chemical_symbols()))
        self.species_atomic_numbers: list = [atomic_numbers[sym] for sym in self.species]
        self.masses: np.ndarray = np.array([atomic_masses[atomic_numbers[sym]] for sym in self.atoms.get_chemical_symbols()])
        self.recip_lat: np.ndarray = 2 * np.pi * np.linalg.inv(self.lat).T

    def read_band_kpt_info(self):
        self.nbnd: int = self.nocc + self.nc
        self.nbnd_begin: int = self.nocc - self.nv
        self.nbnd_end: int = self.nocc + self.nc 
        self.nbnd_red: int = self.nv + self.nc
        self.nk: int = int(np.prod(self.kgrid))
        self.nq: int = self.nk
        self.kpts: Kpts = Kpts.from_kgrid(
            kgrid = self.kgrid,
            is_reduced=False,
        ).epw_kpts
        self.nmodes: int = 3 * self.nat
        self.nk_per_pool: int = int(self.nk / self.npool)

        # Some error checking.
        if self.nbnd_begin > self.nbnd_end:
            raise Exception('Keyword "nbnd_begin" cannot be greater than "nbnd_end".')
        elif self.nbnd_end > self.nbnd:
            raise Exception('Keyword "nbnd_begin" cannot be greater than "nbnd in QE".')
        if self.nbnd_begin < 0:
            raise Exception('Keyword "nbnd_begin" cannot be less than 0.')

    def read_epw(self):
        pass

    def calc_dyn(self):
        pass

    def rotate_to_mode_basis(self):
        pass

    def read(self):
        self.read_structure_info()
        self.read_band_kpt_info()
        self.read_epw()
        self.calc_dyn()
        self.rotate_to_mode_basis()

    def write(self):
        dsets_dict = {
            'key': 'value',
        }

        # Write all the datasets to h5 file.
        with h5py.File('elph.h5', 'w') as f:
            for key, value in dsets_dict.items():
                f.create_dataset(key, data=value)
    
#endregion