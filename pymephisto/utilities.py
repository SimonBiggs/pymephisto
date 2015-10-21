# Copyright (C) 2015 Simon Biggs and Riverina Cancer Care Centre
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public
# License along with this program. If not, see
# http://www.gnu.org/licenses/.


import numpy as np

from scipy.interpolate import interp1d


def normalise_pdd(relative_dose):
    """Normalise pdds to maximum dose
    """
    normalisation = 100 / np.max(relative_dose)
    return relative_dose * normalisation


def normalise_profile_to_cra(distance, relative_dose):
    """Normalise profiles to a distance of 0
    """
    # Linear interpolation function
    interpolation = interp1d(distance, relative_dose)
    normalisation = 100 / interpolation(0)

    return relative_dose * normalisation


def normalise_profile_to_cm(distance, relative_dose):
    """Normalise profiles to the position of centre of mass
    """
    threshold = 0.5 * np.max(relative_dose)
    weights = relative_dose.copy()
    weights[weights < threshold] = 0

    centre_of_mass = np.average(distance, weights=weights)
    normalisation = 100 / interp1d(distance, relative_dose)(centre_of_mass)

    return relative_dose * normalisation


def normalise_profile_to_pdd(distance, relative_dose, scan_depth,
                             pdd_distance, pdd_relative_dose):
    """Normalise profiles to a distance of 0
    """
    pdd_interpolation = interp1d(pdd_distance, pdd_relative_dose)
    cra_normalisation = pdd_interpolation(scan_depth)

    # Linear interpolation function
    profile_interpolation = interp1d(distance, relative_dose)
    normalisation = cra_normalisation / profile_interpolation(0)

    return relative_dose * normalisation
