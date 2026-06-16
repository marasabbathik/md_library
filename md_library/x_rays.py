import matplotlib as mpl
import numpy as np
import re
import os
import xraydb 
import math
from typing import Any, Optional, Union
from datetime import datetime
from datetime import timedelta
import pandas as pd

num = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

def Angles_to_Q(Omega, Two_theta, wavelength):
    k = 2*np.pi / wavelength
    Omega_rad = np.radians(Omega)
    Two_theta_rad = np.radians(Two_theta)
    Qx = k * (np.cos(Two_theta_rad - Omega_rad) - np.cos(Omega_rad))
    Qz = k * (np.sin(Two_theta_rad - Omega_rad) + np.sin(Omega_rad))
    return Qx, Qz


def rigaku_value_finder(variable_name: str, data_list: list) ->list:
    matches = re.findall(rf"{variable_name}.*", data_list)
    values = []
    if matches:
        for i in range(len(matches)):
            value = matches[i].strip().split(None, 1)[1].strip(" \"")
            values.append(value)
    return values

def huber_value_finder(variable_name: str, data_list: list) -> str:
    matches = re.findall(rf"{variable_name}.*", data_list)
    if matches:
        value = matches[0].strip().split("=")[1].strip("';\" ")
        return value
    else: return None

def rigaku_components_settings_finder(variable_name, data_list):
    matches = re.findall(rf"{variable_name}.*", data_list)
    if matches:
        components = []
        for i in range(len(matches)):
            component = matches[i].strip().split(None, 1)[1].strip(" \"")
            components.append(component)
        return components

def value_finder_cif(keyword, cif_data):
    matches = re.findall(rf"{keyword}.*", cif_data)
    if matches:
        if "(" in matches[0]:
            value = matches[0].strip().split()[1].strip().split("(")[0]
        else:
            value = matches[0].strip().split()[1].strip()
    return value

def XRD_peaks_value_finder(keyword, comment):
    matches = re.findall(rf"{keyword}.*", "".join(comment))
    values = []
    units = []
    if matches:
        for i in range(len(matches)):
            value = matches[i].strip().split(":")[1].strip().split()[0]
            value_unit = matches[i].strip().split(":")[1].strip().split()[1]
            values.append(value)
            units.append(value_unit)
    return values, units

def normalize_element_symbol(text):
    match = re.match(r"[A-Za-z]+", str(text).strip())
    if not match:
        return str(text).strip()
    symbol = match.group(0)
    return symbol[0].upper() + symbol[1:].lower()

def cif_loader(cif_file_name, file_location):
    cif_file_name = cif_file_name + ".cif"
    with open(os.path.join(file_location, cif_file_name), 'r') as file:
        data = file.read()
        lines = data.splitlines()
        atom_data = {"atom_label": [], "element": [], "x_coordinate": [], "y_coordinate": [], "z_coordinate": []}
        a = float(value_finder_cif("_cell_length_a", data))
        b = float(value_finder_cif("_cell_length_b", data))
        c = float(value_finder_cif("_cell_length_c", data))
        alpha = float(value_finder_cif("_cell_angle_alpha", data))
        beta = float(value_finder_cif("_cell_angle_beta", data))
        gamma = float(value_finder_cif("_cell_angle_gamma", data))
        lattice_parameters = [a, b, c]
        lattice_angles = [alpha, beta, gamma]        
    
    loop_indexes = [i for i, line in enumerate(lines) if line.startswith("loop_")]
    symetry_operations_start_index = loop_indexes[1]+1
    while lines[symetry_operations_start_index].startswith("_"):
        symetry_operations_start_index += 1
    symetry_operations_end_index = loop_indexes[2] 
    symetry_operations_type = lines[symetry_operations_start_index - 1].strip()
    symetry_operations = lines[symetry_operations_start_index:symetry_operations_end_index]
    atom_positions_start_index = loop_indexes[2] + 1
    if len(loop_indexes) > 3:
        atom_positions_end_index = loop_indexes[3]
    else:
        atom_positions_end_index = None 
    values_names = []
    for line in lines[atom_positions_start_index:atom_positions_end_index]:
        if line.startswith("_"):
            values_names.append(line.strip())
        else:
            line_informations = line.strip().split()
            atom_label = line_informations[values_names.index("_atom_site_label")] if "_atom_site_label" in values_names else ""
            if "_atom_site_type_symbol" in values_names:
                element = normalize_element_symbol(line_informations[values_names.index("_atom_site_type_symbol")])
            else:
                element = normalize_element_symbol(atom_label)

            atom_data["atom_label"].append(atom_label)
            atom_data["element"].append(element)
            atom_data["x_coordinate"].append(float(line_informations[values_names.index("_atom_site_fract_x")].split("(")[0]))
            atom_data["y_coordinate"].append(float(line_informations[values_names.index("_atom_site_fract_y")].split("(")[0]))
            atom_data["z_coordinate"].append(float(line_informations[values_names.index("_atom_site_fract_z")].split("(")[0]))
        
    unique_sites = set()
    expanded = {"atom_label": [], "element": [], "x_coordinate": [], "y_coordinate": [], "z_coordinate": []}
    for i in range(len(atom_data["atom_label"])):
        x0 = atom_data["x_coordinate"][i]
        y0 = atom_data["y_coordinate"][i]
        z0 = atom_data["z_coordinate"][i]
        for operation in symetry_operations:
            # s = operation.lstrip()
            # if s[0].isdigit():
                # operation = operation.split(" ",1)[1]
            quoted_operation = re.search(r"(['\"])(.*?)\1", operation)
            if quoted_operation:
                operation = quoted_operation.group(2)
            operation = operation.strip("'").split(",")
            parts = [value.strip() for value in operation]
            local = {"x": x0, "y": y0, "z": z0}
            x = eval(parts[0], {"__builtins__": {}}, local) % 1.0
            y = eval(parts[1], {"__builtins__": {}}, local) % 1.0
            z = eval(parts[2], {"__builtins__": {}}, local) % 1.0

            key = (atom_data["atom_label"][i], round(x, 6), round(y, 6), round(z, 6))
            if key in unique_sites:
                continue
            unique_sites.add(key)
            expanded["atom_label"].append(atom_data["atom_label"][i])
            expanded["element"].append(atom_data["element"][i])
            expanded["x_coordinate"].append(x)
            expanded["y_coordinate"].append(y)
            expanded["z_coordinate"].append(z)
    return expanded, lattice_parameters, lattice_angles

def interlayer_distance(lattice_constants, lattice_angles, Miller_indexes):
    a, b, c = lattice_constants
    alpha, beta, gamma = np.radians(lattice_angles)
    h, k, l = Miller_indexes
    
    V_square = (a * b * c)**2 * (1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2)
    
    A = h**2*b**2*c**2*np.sin(alpha)**2
    B = k**2*a**2*c**2*np.sin(beta)**2
    C = l**2*a**2*b**2*np.sin(gamma)**2
    D = 2*h*k*a*b*c**2*(np.cos(alpha)*np.cos(beta) - np.cos(gamma))
    E = 2*h*l*a*b**2*c*(np.cos(alpha)*np.cos(gamma) - np.cos(beta))
    F = 2*k*l*a**2*b*c*(np.cos(beta)*np.cos(gamma) - np.cos(alpha))

    d_hkl = np.sqrt(V_square / (A + B + C + D + E + F))
    return d_hkl

def atomic_form_factor(element, hkl, lattice_parameters, lattice_angles, wavelength):
    element = normalize_element_symbol(element)
    d_hkl = interlayer_distance(lattice_parameters, lattice_angles, hkl)
    q = 2 * np.pi / d_hkl

    f0 = xraydb.f0(element, q)
    fp = xraydb.f1_chantler(element, wavelength)
    fpp = xraydb.f2_chantler(element, wavelength)
    return f0 + fp + 1j * fpp

def unitcell_structural_factor(atom_data, lattice_parameters, lattice_angles, Laue_indexes, wavelength):
    h, k, l = Laue_indexes
    Fs_hkl = []
    for i in range(len(atom_data["x_coordinate"])):
        element = atom_data["element"][i]
        x = atom_data["x_coordinate"][i]
        y = atom_data["y_coordinate"][i]
        z = atom_data["z_coordinate"][i]
        f = atomic_form_factor(element, Laue_indexes, lattice_parameters, lattice_angles, wavelength)
        phase = np.exp(2j * np.pi * (h * x + k * y + l * z))
        Fs_hkl.append(f * phase)
    F_hkl_total = np.abs(sum(Fs_hkl))**2
    return F_hkl_total

def two_theta_predicted(lattice_parameters, lattice_angles, lattice_Miller_indexes, measured_two_theta, atom_data, wavelength, tolerance_factor, peaks, num):
    two_theta_predicted = []
    indexes_two_theta_predicted = []
    if lattice_Miller_indexes == "all_possible":
            lattice_Miller_indexes = []
            for i in range(num):
                for j in range(num):
                    for k in range(num):
                        if [i,j,k] == [0,0,0]:
                            continue
                        g = math.gcd(i, math.gcd(j,k))
                        if g == 1:
                            lattice_Miller_indexes.append([i, j, k])
    else:
        lattice_Miller_indexes = [lattice_Miller_indexes]
    for i in range(len(lattice_Miller_indexes)): # range_of_interations
        two_theta, order = 0, 0
        indexes_for_drop = []        
        while two_theta < max(measured_two_theta): 
            order = order + 1
            indexes = [order*index for index in lattice_Miller_indexes[i]]
            d_hkl = interlayer_distance(lattice_parameters, lattice_angles, indexes)
            two_theta = 2*np.degrees(np.arcsin(wavelength / (2 * d_hkl))) 
            peak_match = any(abs(two_theta - peak) <= tolerance_factor for peak in peaks)
            if peak_match:
                F = unitcell_structural_factor(atom_data, lattice_parameters, lattice_angles, indexes, wavelength)
                if two_theta < max(measured_two_theta):
                    if F > 1e-16:
                        two_theta_predicted.append(two_theta)
                        indexes_two_theta_predicted.append(indexes)
                    else:
                        indexes_for_drop.append(i)
                else:
                    indexes_for_drop.append(i)
    two_theta_predicted = [two_theta_predicted[i] for i in range(len(two_theta_predicted)) if i not in indexes_for_drop]
    indexes_two_theta_predicted = [indexes_two_theta_predicted[i] for i in range(len(indexes_two_theta_predicted)) if i not in indexes_for_drop]
    return two_theta_predicted, indexes_two_theta_predicted

def value_finder_XRD(keyword, comment):
    matches = re.findall(rf"{keyword}.*", "".join(comment))
    values = []
    units = []
    if matches:
        for i in range(len(matches)):
            value = matches[i].strip().split(":")[1].strip().split()[0]
            value_unit = matches[i].strip().split(":")[1].strip().split()[1]
            values.append(value)
            units.append(value_unit)
    return values, units

def load_XRD_ras_data(file_names, data_folder_location):
    measurements_components_settings = []
    # measurements = {"sample_name": [], "machine_name": [], "detector_name": [], "Cu_K_alpha1": [], "Cu_K_alpha2": [], "Cu_K_beta": [], "start_date_time": [], "stop_date_time": [], "start_value": [], "end_value": [], "step_value": [], "measurement_speed_setted": [], "measurement_speed_real": [], "measurement_speed_unit": [], "measurement_x_axis_unit": [], "measurement_y_axis_unit": [], "scanning_mode": [], "x_axis": [], "temperature": [], "omega": [], "two_theta": [], "intensity": []}
    measurements = {"sample_name": [], "machine_name": [], "detector_name": [], "Cu_K_alpha1": [], "Cu_K_alpha2": [], "Cu_K_beta": [], "omega_origin": [], "chi_origin": [], "phi_origin": [], "two_theta_origin": [], "two_theta_start": [], "two_theta_stop": [], "two_theta_step": [], "start_date_time": [], "stop_date_time": [], "omega_start": [], "omega_stop": [], "omega_step": [], "measurement_speed_setted": [], "measurement_speed_real": [], "measurement_speed_unit": [], "measurement_x_axis_unit": [], "measurement_y_axis_unit": [], "scanning_mode": [], "x_axis": [], "temperature": [], "omega": [], "two_theta": [], "intensity": []}
    for file_name in file_names:
        with open(os.path.join(data_folder_location, file_name), 'r', encoding="ISO-8859-1") as file:
            data = file.read()
            if "#h.tp=" in data:
                start_date_time = datetime.strptime(huber_value_finder("h.DateTime", data), "%d-%b-%Y %H:%M:%S")
                stop_date_time = start_date_time + timedelta(minutes=30)
                sample_name = huber_value_finder("h.Comment", data)
                machine_name = "Huber"
                detector_name = "1D detector?"
                x_axis = huber_value_finder("h.Title", data).split(" ")[0]
                scanning_mode = "step mode"
                start_value = float(huber_value_finder("h.Title", data).split(",")[1])
                end_value = float(huber_value_finder("h.Title", data).split(",")[2])*2
                step_value = "N/A"
                measurement_x_axis_unit = "deg"
                measurement_speed_setted = "N/A"
                measurement_speed_real = "N/A"
                measurement_speed_unit = "deg/min"

                measurements["sample_name"].append(sample_name)
                measurements["start_date_time"].append(start_date_time)
                measurements["stop_date_time"].append(stop_date_time)
                measurements["machine_name"].append(machine_name)
                measurements["detector_name"].append(detector_name)
                measurements["x_axis"].append(x_axis)
                measurements["scanning_mode"].append(scanning_mode)
                measurements["start_value"].append(start_value)
                measurements["end_value"].append(end_value)
                measurements["step_value"].append(step_value)
                measurements["measurement_x_axis_unit"].append(measurement_x_axis_unit) 
                measurements["measurement_speed_setted"].append(measurement_speed_setted)
                measurements["measurement_speed_real"].append(measurement_speed_real)
                measurements["measurement_speed_unit"].append(measurement_speed_unit)

                components_names = []
                components_internal_names = []
                components_offset = []
                components_values = []
                components_units = []

            elif "RAS_INT_START" in data:
                date_time_format = "%m/%d/%y %H:%M:%S"
                
                number_of_cycles = data.count("RAS_INT_END")
                print(f"Number of cycles: {number_of_cycles}")
                sample_name = rigaku_value_finder("FILE_SAMPLE", data)[0]
                machine_name = rigaku_value_finder("FILE_OPERATOR", data)[0]
                detector_name = rigaku_value_finder("HW_COUNTER_SELECT_NAME", data)[0]
                Cu_K_alpha1 = float(rigaku_value_finder("HW_XG_WAVE_LENGTH_ALPHA1", data)[0])
                Cu_K_alpha2 = float(rigaku_value_finder("HW_XG_WAVE_LENGTH_ALPHA2", data)[0])
                Cu_K_beta = float(rigaku_value_finder("HW_XG_WAVE_LENGTH_BETA", data)[0])
                
                scan_start_values = rigaku_value_finder("MEAS_SCAN_START ", data)
                scan_end_values = rigaku_value_finder("MEAS_SCAN_STOP ", data)
                scan_step_values = rigaku_value_finder("MEAS_SCAN_STEP ", data)
                
                if "*FILE_DATA_TYPE \"RAS_3DE_RSM\"" in data:
                    omega_origin = float(rigaku_value_finder("MEAS_3DE_OMEGA_ORIGIN", data)[0]) #            
                    chi_origin = float(rigaku_value_finder("MEAS_3DE_CHI_ORIGIN", data)[0])
                    phi_origin = float(rigaku_value_finder("MEAS_3DE_PHI_ORIGIN", data)[0])
                    two_theta_start = float(rigaku_value_finder("MEAS_3DE_SCAN_START", data)[0])
                    two_theta_stop = float(rigaku_value_finder("MEAS_3DE_SCAN_STOP", data)[0])
                    two_theta_step = float(rigaku_value_finder("MEAS_3DE_SCAN_STEP", data)[0])
                    number_of_two_thetas = int(np.round((two_theta_stop - two_theta_start) / two_theta_step)) 
                    two_theta = np.round(np.linspace(two_theta_start, two_theta_stop-two_theta_step, number_of_two_thetas),4)
                    omega_relative_start = float(rigaku_value_finder("MEAS_3DE_STEP_START", data)[0])
                    omega_relative_stop = float(rigaku_value_finder("MEAS_3DE_STEP_STOP", data)[0])
                    omega_step = float(rigaku_value_finder("MEAS_3DE_STEP_STEP", data)[0])
                    omega_start = omega_origin + omega_relative_start
                    omega_stop = omega_origin + omega_relative_stop
                    number_of_omegas = int(np.round((omega_stop - omega_start) / omega_step) + 1) 
                    omega = np.round(np.linspace(omega_start, omega_stop, number_of_omegas),4)
                else:
                    omega = float(rigaku_value_finder("MEAS_COND_AXIS_POSITION-2", data)[0])
                    omega_start = omega
                    omega_stop = omega
                    omega_step = 0
                    chi_origin = float(rigaku_value_finder("MEAS_COND_AXIS_POSITION-3", data)[0])
                    phi_origin = float(rigaku_value_finder("MEAS_COND_AXIS_POSITION-4", data)[0])

                    two_theta_starts = scan_start_values
                    two_theta_stops = scan_end_values
                    two_theta_steps = scan_step_values
                    two_theta = np.array([])
                    for i in range(number_of_cycles):
                        two_theta_start = float(two_theta_starts[i])
                        two_theta_stop = float(two_theta_stops[i])
                        two_theta_step = float(two_theta_steps[i])
                        number_of_two_thetas = int(np.round((two_theta_stop - two_theta_start) / two_theta_step) + 1)
                        two_theta_cycle = np.round(np.linspace(two_theta_start, two_theta_stop, number_of_two_thetas),4)
                        two_theta = np.concatenate((two_theta, two_theta_cycle))
                    start_value = two_theta[0]
                    end_value = two_theta[-1]
                    step_value = two_theta[1] - two_theta[0]
                measurement_speed_setted = float(rigaku_value_finder("MEAS_SCAN_SPEED_USER", data)[0])
                measurement_speed_real = float(rigaku_value_finder("MEAS_SCAN_SPEED", data)[0])
                measurement_speed_unit = rigaku_value_finder("MEAS_SCAN_SPEED_UNIT", data)[0]
                measurement_x_axis_unit = rigaku_value_finder("MEAS_SCAN_UNIT_X", data)[0]
                measurement_y_axis_unit = rigaku_value_finder("MEAS_SCAN_UNIT_Y", data)[0]
                scanning_mode = rigaku_value_finder("MEAS_SCAN_MODE", data)[0]
                x_axis = rigaku_value_finder("MEAS_SCAN_AXIS_X", data)[0]
                components_names = rigaku_components_settings_finder("MEAS_COND_AXIS_NAME-", data)
                components_internal_names = rigaku_components_settings_finder("MEAS_COND_AXIS_NAME_INTERNAL-", data)
                components_magicno = rigaku_components_settings_finder("MEAS_COND_AXIS_NAME_MAGICNO-", data)
                components_offset = rigaku_components_settings_finder("MEAS_COND_AXIS_OFFSET-", data)
                components_values= rigaku_components_settings_finder("MEAS_COND_AXIS_POSITION-", data)
                components_units = rigaku_components_settings_finder("MEAS_COND_AXIS_UNIT-", data)
                if "*FILE_DATA_TYPE \"RAS_3DE_RSM\"" in data:
                    print(f"Reciprocal space map was measured")
                    start_date_time = [datetime.strptime(i, date_time_format) for i in rigaku_value_finder("MEAS_SCAN_START_TIME", data)]
                    stop_date_time = [datetime.strptime(i, date_time_format) for i in rigaku_value_finder("MEAS_SCAN_END_TIME", data)]
                    omegas = [float(i) for i in rigaku_value_finder("MEAS_COND_AXIS_POSITION-2 ", data)]
                    measurements["omega"].append(omegas)
                else:    
                    
                    start_date_time = datetime.strptime(rigaku_value_finder("MEAS_SCAN_START_TIME", data)[0], date_time_format)
                    stop_date_time = datetime.strptime(rigaku_value_finder("MEAS_SCAN_END_TIME", data)[0], date_time_format)
                    
                measurements["sample_name"].append(sample_name)
                measurements["machine_name"].append(machine_name)
                measurements["detector_name"].append(detector_name)
                measurements["Cu_K_alpha1"].append(Cu_K_alpha1)
                measurements["Cu_K_alpha2"].append(Cu_K_alpha2)
                measurements["Cu_K_beta"].append(Cu_K_beta)
                measurements["start_date_time"].append(start_date_time)
                measurements["stop_date_time"].append(stop_date_time)

                measurements["phi_origin"].append(phi_origin)
                measurements["chi_origin"].append(chi_origin)
                measurements["two_theta_origin"].append(two_theta_start)
                measurements["two_theta_start"].append(two_theta_start)
                measurements["two_theta_stop"].append(two_theta_stop)
                measurements["two_theta_step"].append(two_theta_step)
                measurements["omega_start"].append(omega_start)
                measurements["omega_stop"].append(omega_stop)
                measurements["omega_step"].append(omega_step)

                #measurements["start_value"].append(start_value)
                #measurements["end_value"].append(end_value)
                #measurements["step_value"].append(step_value)
                measurements["measurement_speed_setted"].append(measurement_speed_setted)
                measurements["measurement_speed_real"].append(measurement_speed_real)
                measurements["measurement_speed_unit"].append(measurement_speed_unit)
                measurements["measurement_x_axis_unit"].append(measurement_x_axis_unit)
                measurements["measurement_y_axis_unit"].append(measurement_y_axis_unit)
                measurements["scanning_mode"].append(scanning_mode)
                measurements["x_axis"].append(x_axis)

            print("\033[31m Sample:\033[0m", sample_name, "\t \033[31m measured on: \033[0m", machine_name, "\t \033[31m with detector: \033[0m", detector_name)
            print("\033[31m Measurement started at: \033[0m", start_date_time, "\t \033[31m ended at: \033[0m", stop_date_time)
            print(x_axis, " scan in ", scanning_mode, " mode")
            print("\033[31m Scan in:\033[0m", start_value, measurement_x_axis_unit, "-", end_value, measurement_x_axis_unit, "\033[31m with step: \033[0m", step_value, measurement_x_axis_unit)
            print("\033[31m With speed: \033[0m", measurement_speed_real, measurement_speed_unit, "(\033[31m setted to: \033[0m", measurement_speed_setted, measurement_speed_unit, ")")
            
            components_settings = pd.DataFrame({"name": components_names, "internal_name": components_internal_names, "offset": components_offset, "value": components_values, "unit": components_units})
            measurements_components_settings.append(components_settings)
            print(components_settings)
            print("###########################################################################################################")
            if machine_name == "Huber":
                pattern_data = re.compile(rf"^\s*{num}\s+{num}\s*")
            else:
                pattern_data = re.compile(rf"^{num}\s{num}\s{num}$")
            if "*FILE_DATA_TYPE \"RAS_3DE_RSM\"" in data:
                # number_of_omegas = len(omegas)
                intensity = []
                # for i in range(number_of_omegas):
                #     two_theta.append([])
                #     intensity.append([])
                # counter_of_lines = 0
                for line in data.splitlines():
                    values = line.strip().split(None)
                    if pattern_data.match(line):
                        #counter_of_lines += 1
                        #omega_index = counter_of_lines % number_of_omegas
                        #two_theta[omega_index].append(float(values[0]))
                        intensity.append(float(values[1]))
                intensity = np.array(intensity).reshape((number_of_omegas, number_of_two_thetas)).tolist()
                Two_theta, Omega = np.meshgrid(two_theta, omega)
                measurements["two_theta"].append(Two_theta)
                measurements["omega"].append(Omega)
            else:
                two_theta = []
                intensity = []
                for line in data.splitlines():
                    if pattern_data.match(line):
                        values = line.strip().split(None)
                        if machine_name == "Huber":
                            two_theta.append(float(values[0])*2)
                        else:   
                            two_theta.append(float(values[0]))
                        intensity.append(float(values[1]))
            measurements["two_theta"].append(two_theta)
            measurements["omega"].append(omega)
            measurements["intensity"].append(intensity)
    if len(measurements["temperature"]) > 1:
        temperature_order = sorted(range(len(measurements["temperature"])), key=lambda k: measurements["temperature"][k])
        print(temperature_order)
        for key in measurements:
            measurements[key] = [measurements[key][i] for i in temperature_order]
        measurements_components_settings = [measurements_components_settings[i] for i in temperature_order]
    return measurements, measurements_components_settings

def load_XRD_txt_data(file_name, file_location):
    with open(os.path.join(file_location, file_name), 'r') as file:
        lines = file.readlines()
        number_of_comment_lines = sum(1 for line in lines if line.startswith('* '))
        comment = lines[0:number_of_comment_lines]
        file.close()
    data = pd.DataFrame(pd.read_csv(os.path.join(file_location, file_name), delimiter='\t', header=0, skiprows=number_of_comment_lines))
    temperatures = []
    for i in range(len(data.columns)):
        match = re.search(num, data.columns[i])
        if match:
            temperatures.append(float(match.group()))
    return data, temperatures, comment

def load_XRD_peaks_table(sample_name, dir_all_data, file_location):
    for file_name in dir_all_data:
        if sample_name in file_name and "peak_table" in file_name and file_name.endswith('.txt'):
            print(f"\033[34m Found file \033[0m {file_name} as file with peaks positions \033[34m for sample \033[0m '{sample_name}':")
            peaks_table = pd.DataFrame(pd.read_csv(os.path.join(file_location, file_name), delimiter='\t', header=0))
    return peaks_table

def load_XRD_peak_positions(file_name, file_location):
    with open(os.path.join(file_location, file_name), 'r') as file:
        lines = file.readlines()
        number_of_comment_lines = sum(1 for line in lines if line.startswith('* '))
        comment = lines[0:number_of_comment_lines]
        file.close()
        data = pd.DataFrame(pd.read_csv(os.path.join(file_location, file_name), delimiter='\t', header=0, skiprows=number_of_comment_lines))
        return data, comment