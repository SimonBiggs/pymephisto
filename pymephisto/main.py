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

import os

import numpy as np
import pandas as pd
import bokeh.plotting as bkh

from .mccread import (
    pull_mephisto_item, pull_mephisto_number, pull_mephisto_data)
from .utilities import (
    normalise_pdd, normalise_profile_to_cra, normalise_profile_to_pdd,
    normalise_profile_to_cm)


def load_mephisto(filepath, output_to_file=True, output_directory=None,
                  normalise_to=None, display=True, sort=True):
    """Input the filepath of a mephisto .mcc file and return the data of the
    scans in four lists, distance, relative_dose, scan_curvetype, and
    scan_depth. Each respective element in these lists corresponds to an
    individual scan.
    """
    # Open the file and store the contents in file_contents
    with open(filepath) as file_pointer:
            file_contents = np.array(file_pointer.readlines())

    # Use the functions defined within mccread.py to pull the desired data
    distance, relative_dose = pull_mephisto_data(file_contents)
    scan_curvetype = pull_mephisto_item('SCAN_CURVETYPE', file_contents)
    scan_depth = pull_mephisto_number('SCAN_DEPTH', file_contents)

    # If user has selected normalise (which is default) normalise the data
    if normalise_to == 'CRA':
        relative_dose = normalise_to_cra(
            distance, relative_dose, scan_curvetype)
    elif normalise_to == 'PDD':
        relative_dose = normalise_to_pdd(
            distance, relative_dose, scan_curvetype, scan_depth)
    elif normalise_to == 'CM':
        relative_dose = normalise_to_cm(
            distance, relative_dose, scan_curvetype)

    if sort:
        sort_ref = np.hstack([
            np.where(scan_curvetype == 'PDD')[0],
            np.where(scan_curvetype == 'INPLANE_PROFILE')[0],
            np.where(scan_curvetype == 'CROSSPLANE_PROFILE')[0]
        ])

        assert len(sort_ref) == len(scan_curvetype)

        distance = list(np.array(distance)[sort_ref])
        relative_dose = list(np.array(relative_dose)[sort_ref])
        scan_curvetype = list(np.array(scan_curvetype)[sort_ref])
        scan_depth = list(np.array(scan_depth)[sort_ref])

    # Output csv's if "output_to_file" is True
    if output_to_file:

        # If user didn't define an output_directory use a default one
        if output_directory is None:
            # Define output directory as a mephisto folder
            filepath_directory = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            output_directory = os.path.join(
                filepath_directory, filename)

        # If the output directory does not exist create it
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        file_output(
            output_directory, distance, relative_dose,
            scan_curvetype, scan_depth)

    # If the parameter display is set to True (default) then display all curves
    if display:
        bokeh_display(distance, relative_dose)

    return distance, relative_dose, scan_curvetype, scan_depth


def bokeh_display(distance, relative_dose):
    """Plot the distance and relative_dose of each scan. Used to display the
    results of the mephisto import.
    """
    # Initialise the plot
    fig = bkh.figure(plot_width=600, plot_height=400)

    # Loop through all data and add a line to the plot
    for i in range(len(relative_dose)):
        # Create a random colour for plotting the line
        random_colour = tuple(np.random.uniform(high=256, size=3).astype(int))

        # Add line to the figure
        fig.line(
            distance[i], relative_dose[i], alpha=0.7, line_width=2,
            line_color=random_colour)

    # Show the figure to the user
    bkh.show(fig)


def normalise_to_cra(distance, relative_dose, scan_curvetype):
    """Normalise the relative_dose of profiles to the CRA and PDDs to dmax.
    """
    normalised_relative_dose = []
    # Loop through all curve types normalising the data
    for i, curvetype in enumerate(scan_curvetype):
        # If curvetype is PDD then use normalise_pdd to normalise
        if curvetype == 'PDD':
            normalised_relative_dose.append(
                normalise_pdd(relative_dose[i]))
        # Otherwise if curvetype is a profile normalise using normalise_profile
        elif (
              (curvetype == 'INPLANE_PROFILE') |
              (curvetype == 'CROSSPLANE_PROFILE')):
            normalised_relative_dose.append(
                normalise_profile_to_cra(distance[i], relative_dose[i]))
        else:
            # Raise an error if the curve type was not as expected
            raise Exception("Unexpected scan_curvetype")

    return normalised_relative_dose


def normalise_to_cm(distance, relative_dose, scan_curvetype):
    """Normalise the relative_dose of profiles to the CRA and PDDs to dmax.
    """
    normalised_relative_dose = []
    # Loop through all curve types normalising the data
    for i, curvetype in enumerate(scan_curvetype):
        # If curvetype is PDD then use normalise_pdd to normalise
        if curvetype == 'PDD':
            normalised_relative_dose.append(
                normalise_pdd(relative_dose[i]))
        # Otherwise if curvetype is a profile normalise using normalise_profile
        elif (
              (curvetype == 'INPLANE_PROFILE') |
              (curvetype == 'CROSSPLANE_PROFILE')):
            normalised_relative_dose.append(
                normalise_profile_to_cm(distance[i], relative_dose[i]))
        else:
            # Raise an error if the curve type was not as expected
            raise Exception("Unexpected scan_curvetype")

    return normalised_relative_dose


def normalise_to_pdd(distance, relative_dose, scan_curvetype, scan_depth):
    """Normalise the relative_dose of PDDs to dmax and the profiles to the
    value of the PDD at the CRA.
    """
    first_pdd_index = np.where(scan_curvetype == 'PDD')[0][0]
    pdd_distance = distance[first_pdd_index]
    pdd_relative_dose = normalise_pdd(relative_dose[first_pdd_index])

    normalised_relative_dose = []
    # Loop through all curve types normalising the data
    for i, curvetype in enumerate(scan_curvetype):
        # If curvetype is PDD then use normalise_pdd to normalise
        if curvetype == 'PDD':
            normalised_relative_dose.append(
                normalise_pdd(relative_dose[i]))
        # Otherwise if curvetype is a profile normalise using normalise_profile
        elif (
              (curvetype == 'INPLANE_PROFILE') |
              (curvetype == 'CROSSPLANE_PROFILE')):
            normalised_relative_dose.append(normalise_profile_to_pdd(
                distance[i], relative_dose[i], scan_depth[i],
                pdd_distance, pdd_relative_dose))
        else:
            # Raise an error if the curve type was not as expected
            raise Exception("Unexpected scan_curvetype")

    return normalised_relative_dose


def file_output(output_directory, distance, relative_dose,
                scan_curvetype, scan_depth):
    """Store the loaded mephisto data into csv files for easy user confirmation
    and use.
    """
    # Determines the filepaths for the output
    filepaths = determine_output_filepaths(
        output_directory, scan_curvetype, scan_depth)

    columns = ['distance (mm)', 'relative dose']

    # Loop over each curvetype and save the data to csv
    for i, curvetype in enumerate(scan_curvetype):
        # Stacks the data into one array and transposes into column orientation
        data = np.vstack([distance[i], relative_dose[i]]).T

        # Use pandas to save data to csv
        df = pd.DataFrame(data, columns=columns)
        df.to_csv(filepaths[i])


def determine_output_filepaths(output_directory, scan_curvetype, scan_depth):
    """Determine a useful filepath for the saving of each mephisto scan.
    """
    filepaths = []

    # Loop over each scan curvetype creating a relevant filepath
    for i, curvetype in enumerate(scan_curvetype):
        if curvetype == 'PDD':
            # Create the filename to be pdd_[number].csv
            filepaths.append(os.path.join(
                output_directory, "pdd_[{0:d}].csv".format(i)))

        elif curvetype == 'INPLANE_PROFILE':
            # Create the filename to be inplaneprofile_depth_[number].csv
            filepaths.append(os.path.join(
                output_directory,
                "inplaneprofile_{0:d}mm_[{1:d}].csv".format(
                    int(scan_depth[i]), i)))

        elif curvetype == 'CROSSPLANE_PROFILE':
            # Create the filename to be crossplaneprofile_depth_[number].csv
            filepaths.append(os.path.join(
                output_directory,
                "crossplaneprofile_{0:d}mm_[{1:d}].csv".format(
                    int(scan_depth[i]), i)))

        else:
            # Raise an error if the curve type was not as expected
            raise Exception("Unexpected scan_curvetype")

    return filepaths
