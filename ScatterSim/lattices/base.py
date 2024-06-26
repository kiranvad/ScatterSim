from copy import deepcopy
import numpy as np

class Lattice:
    """Defines a lattice type and provides methods for adding objects to the
    lattice. This is the starting point of all the crystalline samples.  It is
    recommended here to set the basis vectors (lattice spacing, and alpha,
    beta, gamma) and then have objects inherit this, such as FCC BCC here.
    Finally, you need to add a list of NanoObjects with form factors for this
    to work.

    sigma_D : the DW factor. Can be one number or a two tuple:
            (unit normal, DW factor) pair. The latter is for more complex
            anisotropic DW factors.
    """

    def __init__(
            self,
            objects,
            lattice_spacing_a=1.0,
            lattice_spacing_b=None,
            lattice_spacing_c=None,
            alpha=90,
            beta=None,
            gamma=None,
            sigma_D=0.01,
            lattice_positions=None,
            lattice_coordinates=None,
            symmetry=None,
            lattice_types=None):
        ''' Initialize Lattice object.
                objects : list of NanoObjects
        '''
        # accept just one object too
        if not isinstance(objects, list):
            objects = [objects]

        self.lattice_spacing_a = lattice_spacing_a
        if lattice_spacing_b is None:
            self.lattice_spacing_b = lattice_spacing_a
        else:
            self.lattice_spacing_b = lattice_spacing_b

        if lattice_spacing_c is None:
            self.lattice_spacing_c = lattice_spacing_a
        else:
            self.lattice_spacing_c = lattice_spacing_c

        self.lattice_spacings = np.array([self.lattice_spacing_a,
                                          self.lattice_spacing_b,
                                          self.lattice_spacing_c])

        self.alpha = np.radians(alpha)
        if beta is None:
            self.beta = np.radians(alpha)
        else:
            self.beta = np.radians(beta)
        if gamma is None:
            self.gamma = np.radians(alpha)
        else:
            self.gamma = np.radians(gamma)

        if lattice_positions is None:
            lattice_positions = ['center']
        if lattice_coordinates is None:
            lattice_coordinates = np.array([[0, 0, 0]])
        if lattice_types is None:
            lattice_types = np.ones(len(lattice_coordinates))
        if symmetry is None:
            symmetry = {
                'crystal family': 'N/A',
                'crystal system': 'N/A',
                'Bravais Lattice': 'N/A',
                'crystal class': 'N/A',
                'point group': 'N/A',
                'space group': 'N/A',
            }

        # now set Lattice disorder
        self.sigma_D = np.array(sigma_D)          
        self.symmetry = symmetry
        self.lattice_positions = lattice_positions
        self.lattice_coordinates = lattice_coordinates
        self.lattice_types = lattice_types
        order = np.argsort(self.lattice_types)
        # re-order types
        self.lattice_types = np.array(self.lattice_types)[order]
        self.lattice_coordinates = np.array(self.lattice_coordinates)[order]
        self.lattice_positions = np.array(self.lattice_positions)[order]
        # number_types is number of *unique* objects
        self.number_types = len(np.unique(lattice_types))
        # number_objects is total number objects
        self.number_objects = len(lattice_types)

        unique_types = np.unique(lattice_types)

        # finally make objects
        self.lattice_objects = list()
        if len(objects) == 1:
            # make multiple copies
            obj = objects[0]
            for i in range(self.number_objects):
                self.lattice_objects.append(deepcopy(obj))
        elif len(objects) == self.number_types:
            # TODO : right now assumes all types are ordered, should fix for
            # unordered later
            for i, typ in enumerate(unique_types):
                w = np.where(self.lattice_types == typ)[0]
                if len(w) > 0:
                    for k in range(len(w)):
                        self.lattice_objects.append(deepcopy(objects[i]))
            # move objects into list
        else:
            raise ValueError("Can only handle one or " +
                             "{} objects".format(self.number_types) +
                             ", but received {}.".format(len(objects))
                             )

        # now update the positions
        for i in range(len(self.lattice_objects)):
            pos = self.lattice_coordinates[i] * self.lattice_spacings
            self.lattice_objects[i].set_origin(*pos)
            self.lattice_objects[i].rebuild()

    def update_sigma_D(self, sigma_D):
        '''sigma_D : the DW factor. Can be one number or a two tuple:
        (unit normal, DW factor) pair '''
        self.sigma_D = np.array(sigma_D)

    def unit_cell_volume(self):
        V = np.sqrt(1 -
                    (np.cos(self.alpha))**2 -
                    (np.cos(self.beta))**2 -
                    (np.cos(self.gamma))**2 +
                    2 *
                    np.cos(self.alpha) *
                    np.cos(self.beta) *
                    np.cos(self.gamma))
        V *= self.lattice_spacing_a\
            * self.lattice_spacing_b\
            * self.lattice_spacing_c

        return V

    def sum_over_objects(self, funcname, shape, dtype, *args, **kwargs):
        ''' Sum the function with name 'funcname' over variable 'vec'.
        Forwards other keyword arguments to function

        Parameters
        ---------
        funcname : the function name
        shape : the shape of the result
        dtype : the data type
        args : arguments to the function
        kwargs : keyword arguments to function

        '''
        res = np.zeros(shape, dtype=dtype)
        # cts = 0.

        for obj in self.lattice_objects:
            restmp = getattr(obj, funcname)(*args, **kwargs)
            res += restmp

        return res

    # Components of intensity calculation

    def multiplicity_lookup(self, h, k, l):
        """Returns the peak multiplicity for the given reflection."""
        return 1

    def symmetry_factor(self, h, k, l):
        """Returns the symmetry factor (0 for forbidden)."""
        return 1

    def q_hkl(self, h, k, l):
        """Determines the position in reciprocal space for the given
        reflection."""

        # NOTE: This is assuming cubic/rectangular only!
        qhkl_vector = np.array([2 * np.pi * h / (self.lattice_spacing_a),
                                2 * np.pi * k / (self.lattice_spacing_b),
                                2 * np.pi * l / (self.lattice_spacing_c)])
        
        qhkl = np.sqrt(qhkl_vector[0]**2 + qhkl_vector[1]**2 + qhkl_vector[2]**2)

        return (qhkl, qhkl_vector)


    def iterate_over_hkl_compute(self, max_hkl=6):
        """Returns a sequence of hkl lattice peaks (reflections)."""

        # r will contain the return value, an array with rows that contain:
        # h, k, l, peak multiplicity, symmetry factor, qhkl, qhkl_vector
        r = []

        for h in range(-max_hkl, max_hkl + 1):
            for k in range(-max_hkl, max_hkl + 1):
                for l in range(-max_hkl, max_hkl + 1):

                    # Don't put a reflection at origin
                    if not (h == 0 and k == 0 and l == 0):
                        m = self.multiplicity_lookup(h, k, l)
                        f = self.symmetry_factor(h, k, l)

                        if m != 0 and f != 0:
                            qhkl, qhkl_vector = self.q_hkl(h, k, l)
                            r.append([h, k, l, m, f, qhkl, qhkl_vector])

        return r

    def form_factor(self, qvec):
        """Returns the sum of the form factor of particles in the unit cell."""
        return self.sum_over_objects(
            'form_factor', qvec[0].shape, complex, qvec)

    def projections(self, length, npoints=100):
        """Returns the sum of the form factor of particles in the unit cell."""
        V1 = np.zeros((npoints, npoints))
        V2 = np.zeros((npoints, npoints))
        V3 = np.zeros((npoints, npoints))
        # cts = 0.

        for obj in self.lattice_objects:
            V1tmp, V2tmp, V3tmp = obj.projections(length, npoints=npoints)
            V1 += V1tmp
            V2 += V2tmp
            V3 += V3tmp

        return V1, V2, V3

    def V(self, rvec):
        """Returns the sum of the form factor of particles in the unit cell."""
        return self.sum_over_objects('V', rvec[0].shape, float, rvec)
    
    def Pq(self, q):
        """Returns the sum of the form factor of particles in the unit cell."""
        return self.sum_over_objects('form_factor_squared_isotropic', q.shape, float, q)
    
    def Gq(self, qx, qy=None, qz=None):
        ''' the G(q) from the Debye-Waller factor.
            If sigma_D is a two tuple, then compute the DW factor per unit
            normal.
            if qy and qz none, qx treated as magnitude.

            Note : normalized to lattice_spacing for now for backwards
            compatibility
            so sigma = 1/(2sd ld)
        '''
        # just 1D Debye Waller
        if self.sigma_D.ndim == 0:
            if qy is not None:
                qnorm = np.sqrt(qx**2 + qy**2 + qz**2)
            else:
                qnorm = qx
            res = np.exp(-(self.sigma_D**2) * ((qnorm) **
                                               2 * self.lattice_spacing_a**2))
        else:
            # q must be 3 vector
            if qy is None:
                raise ValueError(
                    "Error: if DW factor is anisotropic, need to "
                    "supply qy and qz")
            res = 1.
            for sd, unrmx, unrmy, unrmz in self.sigma_D:
                # need to fix later, should be array form
                u = np.array([unrmx, unrmy, unrmz])
                res *= np.exp(-(sd**2) * (np.tensordot(
                    np.array([qx, qy, qz]), u, axes=(0, 0))**2) *
                    (self.lattice_spacing_a**2))
        return res

    def beta_ratio(self, q):
        """Returns the beta ratio: |<U(q)>|^2 / <|U(q)|^2>
        for the lattice."""
        # numerator over denominator
        beta_num = self.sum_over_objects('beta_numerator', q.shape, float, q)
        beta_denom = self.Pq(q)
        beta = beta_num / beta_denom

        return beta

    def Z0(self, q, peak, background=None, max_hkl=6):
        """Return Z0(q)
        """
        summation = 0
        hkl_list = self.iterate_over_hkl_compute(max_hkl=max_hkl)
        for h, k, l, m, f, qhkl, qhkl_vector in hkl_list:

            fs = self.form_factor(qhkl_vector)
            fs_squared = (fs * fs.conjugate()).real
            L = peak(q - qhkl)

            summation += (m * (f**2)) * fs_squared * L
    
        out = (1.0 / (q**2)) * summation

        if background is None:
            return out
        else:
            return background(q) + out

    # Outputs
    def to_string(self):
        """Returns a string describing the lattice."""
        s = "Lattice of type: " + self.__class__.__name__ + "\n"

        s += "    Family: {}   ".format(self.symmetry['crystal family'])
        s += "System: {}   ".format(self.symmetry['crystal system'])
        s += "Bravais: {}  ".format(self.symmetry['Bravais lattice'])
        s += " Class: {}   ".format(self.symmetry['crystal class'])
        s += "Space Group: {}\n".format(self.symmetry['space group'])

        s += "    (a, b, c) = ({:.3f},".format(self.lattice_spacing_a)
        s += " {:.3f},".format(self.lattice_spacing_b)
        s += " {:.3f}) in nm\n".format(self.lattice_spacing_c)

        s += "    (alpha, beta, gamma) = (%.3f,%.3f,%.3f) in radians\n" % (
            self.alpha, self.beta, self.gamma)
        s += "                        "
        s += "({:.2f}, ".format(np.degrees(self.alpha))
        s += "{:.2f}, ".format(np.degrees(self.beta))
        s += "{:.2f} in degrees".format(np.degrees(self.gamma))
        s += "    volume = %.4f nm^3\n\n" % self.unit_cell_volume()
        s += "    Objects:\n"
        for i, obj in enumerate(self.objects):
            if i < len(self.positions):
                pos = self.positions[i]
            else:
                pos = '---'
                s += "        {:d} ({})\t".format(i, pos)
                s += "{} ({})\n".format(obj.__class__.__name__,
                                        obj.to_short_string())
        s += "    Unit cell:\n"
        for i, pos, xi, yi, zi, obj in self.iterate_over_objects():
            objstr = obj.__class__.__name__
            s += "        %d (%s)\n" % (i, pos)
            s += "           \t%s (%s)\n" % (objstr, obj.to_short_string())
            s += "           \t  at pos = (%.3f,%.3f,%.3f)\n" % (xi, yi, zi)
        return s

class ExtendedLattice(Lattice):
    ''' An extended lattice of objects, can choose random fillings.
        fill_probs : list of probabilities
    '''

    def __init__(
            self,
            objects,
            lattice_spacing_a=1.0,
            fill_probs=None,
            lattice_spacing_b=None,
            lattice_spacing_c=None,
            alpha=90,
            beta=None,
            gamma=None,
            sigma_D=0.01,
            lattice_positions=None,
            lattice_coordinates=None,
            symmetry=None,
            lattice_types=None,
            repeat=3):
        ''' For now only supports one type.'''
        if fill_probs is None:
            fill_probs = [1]

        # sub_lattice_positions = lattice_positions
        sub_lattice_coordinates = lattice_coordinates
        # sub_lattice_types = lattice_types
        if sub_lattice_coordinates is None:
            sub_lattice_coordinates = [[0, 0, 0]]
        lattice_positions = []
        lattice_coordinates = []
        lattice_objects = []
        lattice_types = []

        # fill_prob = fill_probs[0]

        subunit_x, subunit_y, subunit_z = 1., 1., 1.
        # TODO also add positions and types
        for i in range(len(sub_lattice_coordinates)):
            xp, yp, zp = sub_lattice_coordinates[i]

            for ix in range(repeat):
                xi = ix + xp * subunit_x
                for iy in range(repeat):
                    yi = iy + yp * subunit_y
                    for iz in range(repeat):
                        zi = iz + zp * subunit_z

                        if(np.random.uniform(0.0, 1.0) <= fill_probs[i]):
                            lattice_positions.append('N/A')
                            lattice_types.append(i)
                            lattice_coordinates.append((xi, yi, zi))
                            lattice_objects.append(objects[i])

        super().__init__(
            objects,
            lattice_coordinates=lattice_coordinates,
            lattice_types=lattice_types,
            lattice_positions=lattice_positions,
            lattice_spacing_a=lattice_spacing_a,
            lattice_spacing_b=lattice_spacing_b,
            lattice_spacing_c=lattice_spacing_c,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            sigma_D=sigma_D)

    def symmetry_factor(self, h, k, l):
        """Returns the symmetry factor (0 for forbidden)."""
        return 1

    def unit_cell_volume(self):
        return self.lattice_spacing_a**3