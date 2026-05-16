import matplotlib as mpl
import numpy as np
import re
import os
import xraydb 
import math

def ColorFade(c1, c2, mix, temp, T_c):
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    if T_c and temp == int(T_c):
        return "black"
    else:
        return mpl.colors.to_hex((1-mix) * c1 + mix * c2)

def GetLabel(temp, T_c):
    if T_c is not None and temp == int(T_c):
        return r"$T_\mathrm{C} = $" + str(temp)
    else:
        return str(temp)


def rigaku_value_finder(variable_name, data_list):
    matches = re.findall(rf"{variable_name}.*", data_list)
    if matches:
        values = []
        for i in range(len(matches)):
            value = matches[i].strip().split(None, 1)[1].strip(" \"")
            values.append(value)
        return values

def huber_value_finder(variable_name, data_list):
    matches = re.findall(rf"{variable_name}.*", data_list)
    if matches:
        value = matches[0].strip().split("=")[1].strip("';\" ")
        return value
    
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
    symetry_operations_start_index = loop_indexes[1] + 2
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
    #for i in range(len(atom_data["x_coordinate"])):
        #print(f"Atom {i}: x={atom_data['x_coordinate'][i]}, y={atom_data['y_coordinate'][i]}, z={atom_data['z_coordinate'][i]}")
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
                F = unitcell_structural_factor(atom_data, lattice_parameters, lattice_angles, indexes)
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
