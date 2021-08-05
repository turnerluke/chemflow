import warnings
import numpy as np
import numpy.ma as ma
from scipy.optimize import fsolve
import pandas as pd
import pint
# Self-Made Imports
from ChemFormula import formula_to_MW_atoms
import ChemFormula

# Additional settings
np.seterr(divide='ignore')

# Initialize a unit registry
u = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)

# TODO: Mass Yield Reactor
# TODO: Heater
# TODO: EOS, Activity Coeff
# TODO: Atomic Balances


class process(object):

    def __init__(self, comps:list, streams:list, **props):
        """

        :param comps: type:list
        A list of strings corresponding to the components in the system.
        :param streams: type:list
        A list of strings corresponding to the process stream names.
        :param props: type:dict
        Dictionary of properties and additional specifications.
        Supported keys: MW, Cp, T, P, formulas.
        MW:
        Cp:
        T:
        P:
        formulas: Bool, If comps are chemical formuas, collects MW and atom count info from ChemFormula
        """
        # Unpack components and streams inputs
        self.comps = comps
        self.streams = streams

        # Unpack Properties
        self.MW = np.array(props.get("MW")) * u('g / mol') if props.get("MW") is not None else None
        self.Cp = props.get("Cp")  # TODO: Units & etc
        self.T = dict(zip(streams, [u(Ti) for Ti in props.get("T")])) if props.get("T") is not None else None
        self.P = dict(zip(streams, [u(Pi) for Pi in props.get("P")])) if props.get("P") is not None else None
        if props.get("formulas"):
            self.MW, self.atoms = zip(*[formula_to_MW_atoms(comp_i) for comp_i in self.comps])
            self.MW = np.array(self.MW) * u('g / mol')
        else: self.atoms = None


        # TODO: Throw Errors

        # Initialize the mass & mol flow dictionaries of vectors
        empty_vector = np.zeros(len(comps))
        self.mass_flows = {}
        self.mol_flows = {}
        for stream in streams:
            self.mass_flows[stream] = empty_vector * u('kg / hr')
            self.mol_flows[stream] = empty_vector * u('mol / hr')


    def specify_flow(self, stream:str, comp:str, val:float, units:str):
        if u(units).check("[substance] / [time]"):
            self.mol_flows[stream][[comp == comp_i for comp_i in self.comps]] = val * u(units)
            self.update_mass(stream)
        elif u(units).check("[mass] / [time]"):
            self.mass_flows[stream][[comp == comp_i for comp_i in self.comps]] = val * u(units)
            self.update_mol(stream)
        else:
            raise Exception("Units must be substance or mass per time.")

    # TODO: Specify T&P
    def specify_stream(self, stream:str, vals:list, units:str):
        """
        Allows the user to specify all flows in a single stream. Can be molar or mass, depending on the units.

        :param stream: type:str
        Stream name with the specified flows.
        :param vals: type:list
        List of values corresponding to flowrates.
        :param units: type:str
        Units of flows. All flows must be the same units. The units must be supported by the pint package.
        :return:
        """
        if len(vals) != len(self.comps):
            raise ValueError("Flowrate values must be the same length as the component list.")
        if u(units).check("[substance] / [time]"):
            self.mol_flows[stream][:] = np.array(vals) * u(units)
            self.update_mass(stream)
        elif u(units).check("[mass] / [time]"):
            self.mass_flows[stream][:] = np.array(vals) * u(units)
            self.update_mol(stream)
        else:
            raise Exception("Units must be substance or mass per time.")

    def update_mol(self, stream):
        """
        Updates molar flows of the process, after a mass flow change.  Does nothing if molecular weights are not specified.
        :param stream: Type:str
        Stream name to be updated.
        """
        if self.MW is not None:
            self.mol_flows[stream] = self.mass_flows[stream] / self.MW
        else:
            pass

    def update_mass(self, stream):
        """
        Updates mass flows of the process, after a molar flow change.  Does nothing if molecular weights are not specified.
        :param stream: Type:str
        Stream name to be updated.
        """
        if self.MW is not None:
            self.mass_flows[stream] = self.mol_flows[stream] * self.MW
        else:
            pass

    # TODO: Add T&P
    def add_streams(self, new_streams:list):
        """
        Adds stream names to the total streams of the process, initializing their flow vectors as zero.
        :param new_streams: Type:list
        A list of strings to add to the stream list.
        :return:
        """
        if new_streams == []:
            return
        # Iterate through new specified streams and append to streams list
        empty_vector = np.zeros(len(self.comps))
        for stream in new_streams:
            if stream in self.streams:
                raise Exception("Stream name {} has already been specified.".format(stream))
            else:
                self.streams.append(stream)
                self.mol_flows[stream] = empty_vector * u("mol / hr")
                self.mass_flows[stream] = empty_vector * u("kg / hr")


    def stoich_reaction(self, stream_in:str, stream_out:str, nu, spec=None):
        # TODO: Fix the requirement to adjust units to base units
        """
        Calculates a theoretical outlet stream, given an inlet stream, and all reaction specifications.
        :param stream_in: type:str
        The name of the stream to be the inlet of the reactor.
        :param stream_out: type:str
        The name of the stream to adjust to being the outlet of the reactor.
        :param nu: type: np.ndarray, list, list of lists
        The stoichiometric array.
        :param basis: type:dict
        A dictionary contining the keys: Basis, Vals, Units, Conversion Components
        Spec - Possible values: Extent, Conversion
        Vals - np.ndarray vector or list of values corresponding to the reaction basis extents
        Units - Only required if extent is specified. String of units of type [substance] / [time].
        Components - Optional, when a basis of conversion is specified, a list of component names of the same
        size of reactions to base the fractional conversion upon.
        If nothing is specified, it is assumed to be a baisis of fractional conversions equal to one for each limiting
        reagent for each reaction.

        :return:
        """
        # Unpack nu into a numpy array
        nu = self.unpack_nu(nu)

        # Unpack the specified basis
        if spec == None:  # If no specifications are provided, the reaction is assumed to go to 100% conversion of the limiting reagent in each reaction.
            frac_conv = np.ones(nu.shape[0])
            frac_basis = True
        elif type(spec) == dict:
            if spec["Basis"] == "Extent":
                if type(spec["Vals"]) == list:
                    ext = np.array(basis["Vals"])
                elif type(spec["Vals"]) == np.ndarray:
                    ext = spec["Vals"]
                else:
                    raise TypeError("Basis value must be a list or numpy array.")

                if u(spec["Units"]).check("[substance] / [time]"):
                    ext *= u(spec["Units"])
                else:
                    raise ValueError(
                        "Basis units must be of type [substance] / [time] for extent of reaction specification.")
                frac_basis = False

            elif spec["Basis"] == "Conversion":
                try:
                    if type(spec["Vals"]) == list:
                        frac_conv = np.array(spec["Vals"])
                    elif type(spec["Vals"]) == np.ndarray:
                        frac_conv = spec["Vals"]
                except:
                    frac_conv = np.ones(nu.shape[0])
                frac_basis = True
            else:  # Basis is not specified to be extent or fractional conversion
                raise ValueError("Basis must be Extent or Conversion.")

        # Check for proper shapes
        if frac_basis:
            if frac_conv.shape[0] != nu.shape[0]:
                raise ValueError("Fractional conversion must have the same first dimension as the stoichiometric array.")
        else:
            if ext.shape[0] != nu.shape[0]:
                raise ValueError("Extent of reaction must have the same first dimension as the stoichiometric array.")

        if nu.shape[1] != len(self.comps):
            raise ValueError("The second dimension of the stoichiometric matrix must be the same as the number of components.")

        # Unpack flow in for ease of calculation
        n0 = self.mol_flows[stream_in]

        # Calculate extents if the basis is fractional conversion:
        if frac_basis:
            try:  # If component conversion is specified
                lim = spec["Components"]
                lim_idx = [[i if lim_reag == comp_i else np.NaN for i, comp_i in enumerate(comps)] for lim_reag in lim]
                lim_idx = np.array(lim_idx)
                lim_idx = lim_idx(np.logical_not(np.isnan(x)))
            except:  # If no component conversion is specified, pick the limiting reagents in each reaction
                reac_contribution = ma.masked_array(ma.masked_invalid(-n0 / nu), mask= (nu>=0))  #
                lim_idx = np.argmin(reac_contribution, axis=1)

            ext = np.array([(-frac_conv[i] * n0[lim_comp_i] / nu[i, lim_comp_i]).to_base_units().magnitude for i, lim_comp_i in enumerate(lim_idx)]) * u("mol / s")
            #print(ext)
        # With calculated extent: calculate change in moles of each species
        dn = np.sum(nu * ext[np.newaxis].T, axis=0)
        # Calc moles out from moles in and change
        nf = n0 + dn

        # Check for errors in reaction specifications
        if np.any(nf<0):
            raise ValueError("Reaction specifications are impossible, review specifications.")

        self.mol_flows[stream_out] = nf
        self.update_mass(stream_out)

        # Check the mass balance in & out of reactor
        if not np.isclose(np.sum(self.mass_flows[stream_in]), np.sum(self.mass_flows[stream_in])):
            warnings.warn("The mass balance around this reactor does not close.")

    def equilib_reaction(self, K, nu, stream_in:str, stream_out:str):
        """

        :param K:
        :param nu:
        :return:
        """
        # Unpack nu into numpy array
        nu = self.unpack_nu(nu)
        # Unpack K in to a numpy vector
        if type(K) == list:
            K = np.array(K)
        if type(K) == np.ndarray:
            if K.ndim != 1:
                raise ValueError("K must be a 1 dimensional vector.")
        else:
            raise TypeError("K must be a list, or numpy vector.")

        # Solve for extents of reation
        # Unpack arguments
        n0 = self.mol_flows[stream_in]
        args = (K, n0, nu, self.P[stream_out], n0.units)
        ext_g = np.ones_like(K)
        theta = fsolve(equilib_residuals, ext_g, args=args, full_output=True)
        if theta[2] == 1:
            ext = theta[0]
            ext *= n0.units
        else:
            raise ValueError("Extent of reaction could not be solved in the equilibrium reactor.")

        # Calculate delta n and flows out from extent
        dn = np.sum(nu * ext[np.newaxis].T, axis=0)
        nf = n0 + dn

        # Check for errors in reaction specifications
        if np.any(nf < 0):
            raise ValueError("Reaction specifications are impossible, review specifications.")

        self.mol_flows[stream_out] = nf
        self.update_mass(stream_out)

        # Check the mass balance in & out of reactor
        if not np.isclose(np.sum(self.mass_flows[stream_in]), np.sum(self.mass_flows[stream_in])):
            warnings.warn("The mass balance around this reactor does not close.")

    #def mass_yield_reaction(self, nu):

    def splitter(self, stream_in:str, streams_out:list, mass_frac:list):
        """
        Performs a splitting operation, from the stream_in to streams_out with mass fractions specified.
        :param stream_in: type:str
        String specifing the inlet stream to be split.
        :param streams_out: type:list
        A list of streams out to recieve fractions of the inlet stream. Adds the streams if they have not yet been added.
        :param mass_frac: type:list
        A list of fractional split each stream out recieves. If they do not sum to one, they are normalized to do so.
        :return:
        """
        # Ensure stream_in is an already specified stream and has mass flows
        if not stream_in in self.streams:
            raise ValueError("{} is not a specified stream.".format(stream_in))
        if np.all(self.mass_flows[stream_in] == 0):
            raise ValueError("{} does not have a mass flow specified.".format(stream_in))
        # Add any unspecified streams
        new_streams = []
        for stream_out in streams_out:
            if not stream_out in self.streams:
                new_streams.append(stream_out)
        self.add_streams(new_streams)
        # Ensure streams_out is the same size as mass_frac
        if len(stream_out) != len(mass_frac):
            raise ValueError("Mass fraction specifications must be the same size as the streams out.")
        # Normalize mass fractions
        if sum(mass_frac) != 1:
            warnings.warn("The mass fractions do not add to 1. Values are assumed to be ratios and normalized.")
            mass_frac = [frac_i / sum(mass_frac) for frac_i in mass_frac]

        # Loop through the streams out, calculating their mass flows
        for stream_out, frac in zip(streams_out, mass_frac):
            self.mass_flows[stream_out] = frac * self.mass_flows[stream_in]
            self.update_mol(stream_out)




    def report_flows(self, units:str="mol / hr", display:bool=False):
        table = {}
        for stream in self.streams:
            if u(units).check("[substance] / [time]"):
                table[stream] = self.mol_flows[stream].to(units).magnitude
            elif u(units).check("[mass] / [time]"):
                table[stream] = self.mass_flows[stream].to(units).magnitude

        table = pd.DataFrame(table, index=self.comps, columns=self.streams)
        table.index.name = units
        if display:
            print(table)
        return table

    def check_atom_bal(self, streams_in:list, streams_out:list, tol:float=1E-2, tol_units:str="mol / s"):
        """

        :param streams_in: type: list
        A list of stream names into a definable system.
        :param streams_out: type: list
        A list of stream names out of a definable system.
        :param tol: type: float
        Tolerance of atom flow difference, with units mol/s.
        :return:
        """
        # Errors: wrong stream name, components don't have atom counts
        for stream in [*streams_in, *streams_out]:
            if stream not in self.streams:
                raise ValueError("{} is not an initialized stream.".format(stream))
        if self.atoms is None:
            raise ValueError("Atom counts of components are not known, please initalize components with formulas.")

        # Initialize dictionaries of atomic flows
        atoms_in = dict.fromkeys(ChemFormula.ATOMS.keys(), 0 * u('mol / s'))
        atoms_out = dict.fromkeys(ChemFormula.ATOMS.keys(), 0 * u('mol / s'))

        n_ins = [self.mol_flows[stream_in] for stream_in in streams_in]
        n_outs = [self.mol_flows[stream_out] for stream_out in streams_out]
        for i, comp in enumerate(self.comps):
            for atom in ChemFormula.ATOMS.keys():
                for n_in in n_ins:
                    atoms_in[atom] += n_in[i] * self.atoms[i][atom]
                for n_out in n_outs:
                    atoms_out[atom] += n_out[i] * self.atoms[i][atom]



        status = True
        unbalanced_atoms = []
        total_discrepancy = 0
        for atom in ChemFormula.ATOMS.keys():
            discrepancy = abs(atoms_in[atom] - atoms_out[atom])
            total_discrepancy += discrepancy
            if discrepancy.to(tol_units).magnitude >= tol:
                status = False
                unbalanced_atoms.append(atom)

        res = {
            "status": "success" if status else "failure",
            "atoms_in": {atom:num for atom,num in atoms_in.items() if num!=0},
            "atoms_out": {atom:num for atom,num in atoms_out.items() if num!=0},
            "mesg": "Atom balance " + ("succeeds" if status else "fails") + " with total discrepancy of {}.".format(
                total_discrepancy)
        }
        return res

    def check_mass_bal(self, streams_in:list, streams_out:list, tol:float=1E-2, tol_units='kg/s'):
        """

        :param streams_in: type:list
        List of stream names into a definable system, of which to check the mass balance of.
        :param streams_out: type:list
        List of stream names out of a definable system, of which to check the mass balance of.
        :param tol: type:float
        Tolerance of mass flow difference, with units kg/s.
        :return:
        """
        # Initialize scalars of masses in & out to check.
        m_in = 0 * u("kg / s")
        m_out = 0 * u("kg / s")
        tol *= u(tol_units)
        for stream_in in streams_in:
            if stream_in in self.streams:
                m_in += sum(self.mass_flows[stream_in])
            else:
                raise ValueError('"{}" is not a defined stream.'.format(stream_in))
        for stream_out in streams_out:
            if stream_out in self.streams:
                m_out += sum(self.mass_flows[stream_out])
            else:
                raise ValueError('"{}" is not a defined stream.'.format(stream_out))

        discrepancy = abs(m_in - m_out)
        status = discrepancy <= tol

        res = {
            "status": "success" if status else "failure",
            "m_in": m_in,
            "m_out": m_out,
            "mesg": "Mass balance " + ("succeeds" if status else "fails") + " with discrepancy of {}.".format(discrepancy)
        }
        return res

    # Assistant Static Methods
    @staticmethod
    def unpack_nu(nu):
        """
        Unpacks a user input stoichiometric matrix into the desired numpy array.
        :param nu:
        :return:
        """
        if type(nu) == list:
            nu = np.array(nu, ndmin=2)
        elif type(nu) == np.ndarray:
            pass
        else:
            raise TypeError("The stoichiometric array must be a list, list of lists, or numpy.ndarray.")
        return nu


# Functions
def equilib_residuals(ext, *args):
    """
    A function to use with scipy.optimize.fsolve to iteratively solve for the extents of reaction, provided the inlet molar flows,
    stoichiometric array, and pressure (if specified).
    :param ext:
    :param args:
    :return:
    """
    K = args[0]
    n0 = args[1]
    nu = args[2]
    P = args[3]
    ext_units = args[4]

    ext *= ext_units

    # Calc moles out with the extent
    dn = np.sum(nu * ext[np.newaxis].T, axis=0)
    # Calc moles out from moles in and change
    nf = n0 + dn
    y = nf / np.sum(nf)
    if P is not None:
        P_adjustment = P / (101.325 * u("kPa"))
    else:
        P_adjustment = 1

    return K - np.prod(y**nu, axis=1) * P_adjustment**np.sum(nu, axis=1)  # TODO: Thermodynamic models and activity coeff
