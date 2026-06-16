import os
import re
import pandas as pd
import scipy.constants as sccs
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl
import numpy as np
from pathlib import Path

h = sccs.h
c = sccs.c
e = sccs.e

num = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

def E_to_cm(E):
    return E *e / (h * c) * 1e-2
def E_to_nm(E):
    return (h * c) / (E * e) * 1e9
def cm_to_E(cm):
    return cm * h * c / e * 1e2
def cm_to_nm(cm):
    return 1/cm
def nm_to_E(nm):
    return (h * c) / (nm * e) * 1e9
def nm_to_cm(nm):
    return 1/nm

def pseudodielectric_function(psi, delta, AOI):
    rho = np.tan(np.deg2rad(psi)) * np.exp(-1j * np.deg2rad(delta))
    pseude_varepsilon =  np.sin(np.deg2rad(AOI))**2 * (1 + np.tan(np.deg2rad(AOI))**2 * ((1 - rho) / (1 + rho))**2)
    return np.real(pseude_varepsilon), np.imag(pseude_varepsilon)

def load_data(file_names, file_location, energy_units):  
    def search_in_list(string, list_of_strings):
        for item in list_of_strings:
            if string in item:
                return item
            
    measured_datas, measurement_temperatures = [], []  
    for file_name in file_names:
        with open(os.path.join(file_location, file_name), 'r') as file:
            lines = file.readlines()
            # i = 0
            # while i < len(lines) and len(lines[i]) == 0:
            #     i += 1
            # lines = lines[i:]
            # first_non_empty_line = next((index for index, line in enumerate(lines) if line.strip()), len(lines))
            # lines = lines[first_non_empty_line:]
            if "VASEmethod" not in lines[0]:
                commands = lines[0].strip().split(", ")
                for command in commands:
                    if "Setpoint" in command or "T=" in command:
                        measurement_temperature = int(re.findall(r"\d+", command)[0])
                        measurement_temperatures.append(measurement_temperature)
                    else:
                        measurement_temperatures.append(300)
            else:
                measurement_temperatures.append(300)
            for line in lines:
                if "VASEmethod" in line:
                    settings = line.strip().removeprefix("VASEmethod[").removesuffix("]").split(", ")
                    Ellipsometer_type = int(re.findall(r"\d+", search_in_list("EllipsometerType", settings))[0])
                    value = search_in_list("AutoSlit", settings)
                    if value is not None:
                        Slit = re.search(r"\d+", value)
                    value = search_in_list("Revs", settings)
                    if value is not None:
                        Revolutions = re.search(r"\d+", value)
                    value = search_in_list("MaxRevs", settings)
                    if value is not None:
                        Max_revolutions = re.search(r"\d+", value)
                if line.strip() in energy_units:
                    energy_unit = line.strip()
                    break
            Energy, AOI, Psi, Delta, Psi_error, Delta_error, Depolarization, Depolarization_error = [], [], [], [], [], [], [], []
            Ellipsometry_data = {"eV": [], "1/cm": [], "nm": [], "AOI": [], "Psi": [], "Delta": [], "Psi_error": [], "Delta_error": [], "Depolarization": [], "Depolarization_error": [], "pseudo_e1": [], "pseudo_e2": [], "AnE_Psi": [], "AnE_Delta": [], "AnE_Psi_error": [], "AnE_Delta_error": [], "Aps_Psi": [], "Aps_Delta": [], "Aps_Psi_error": [], "Aps_Delta_error": [], "Asp_Psi": [], "Asp_Delta": [], "Asp_Psi_error": [], "Asp_Delta_error": [], "mm":  [[[] for _ in range(4)] for _ in range(4)], "mm_error":  [[[] for _ in range(4)] for _ in range(4)]}
           
            if Ellipsometer_type == 0:
                print(f"Ellipsometer type: FIR-ellipsometer")
                pattern_data = re.compile(rf"^{num}\t{num}\t{num}\t{num}\t{num}\t{num}$")
                for line in lines:
                    if pattern_data.match(line):
                        data = line.strip().split("\t")
                        Energy.append(float(data[0]))
                        AOI.append(float(data[1]))
                        Psi.append(float(data[2]))
                        Delta.append(float(data[3]))
                        Psi_error.append(float(data[4]))
                        Delta_error.append(float(data[5]))
            elif Ellipsometer_type == 3:
                print(f"Ellipsometer type: IR-VASE")
                pattern_data = re.compile(rf"^E \t{num} \t{num} \t{num} \t{num} \t{num} \t{num}$") #  \t{num} \t{num} for some measuremets only 6 columns are present, for some 8 columns are present. dont knot what the last two are for
                patern_depol = re.compile(rf"^dpolE \t{num} \t{num} \t{num} \t{num}$") #  \t{num} \t{num}
                for line in lines:
                    if pattern_data.match(line):
                        data = line.strip().split("\t")
                        Energy.append(float(data[1]))
                        AOI.append(float(data[2]))
                        Psi.append(float(data[3]))
                        Delta.append(float(data[4]))
                        Psi_error.append(float(data[5]))
                        Delta_error.append(float(data[6]))
                    if patern_depol.match(line):
                        depol_data = line.strip().split("\t")
                        Depolarization.append(float(depol_data[3]))
                        Depolarization_error.append(float(depol_data[4]))
            elif Ellipsometer_type == 4:
                print(f"Ellipsometer type: RC2")
                pattern_data = re.compile(rf"^{num}\t{num}\t{num}\t{num}\t{num}\t{num}$")
                pattern_depol = re.compile(rf"^dPolE\t{num}\t{num}\t{num}\t{num}$")
                for line in lines:
                    if pattern_data.match(line):
                        data = line.strip().split("\t")
                        Energy.append(float(data[0]))
                        AOI.append(float(data[1]))
                        Psi.append(float(data[2]))
                        Delta.append(float(data[3]))
                        Psi_error.append(float(data[4]))
                        Delta_error.append(float(data[5]))
                    if pattern_depol.match(line):
                        depol_data = line.strip().split("\t")
                        Depolarization.append(float(depol_data[3]))
                        Depolarization_error.append(float(depol_data[4]))
            elif Ellipsometer_type == 5:
                print(f"Ellipsometer type: V-VASE")
                pattern_data = re.compile(rf"^E\t{num}\t{num}\t{num}\t{num}\t{num}\t{num}$")
                pattern_depol = re.compile(rf"^dpolE\t{num}\t{num}\t{num}\t{num}$")
                for line in lines:
                    if pattern_data.search(line):
                        data = line.strip().split("\t")
                        if data[0].isalpha():
                            index = 1
                        else:
                            index = 0
                        Energy.append(float(data[index]))
                        AOI.append(float(data[index+1]))
                        Psi.append(float(data[index+2]))
                        Delta.append(float(data[index+3]))
                        Psi_error.append(float(data[index+4]))
                        Delta_error.append(float(data[index+5]))
                    if pattern_depol.match(line):
                        depol_data = line.strip().split("\t")
                        Depolarization.append(float(depol_data[3]))
                        Depolarization_error.append(float(depol_data[4]))
            elif Ellipsometer_type == 9:
                print(f"Ellipsometer type: RC2")
                pattern_data = re.compile(rf"^E\t{num}\t{num}\t{num}\t{num}\t{num}\t{num}$")
                pattern_depol = re.compile(rf"^dPolE\t{num}\t{num}\t{num}\t{num}$")
                for line in lines:
                    if pattern_data.match(line):
                        data = line.strip().split("\t")
                        Energy.append(float(data[1]))
                        AOI.append(float(data[2]))
                        Psi.append(float(data[3]))
                        Delta.append(float(data[4]))
                        Psi_error.append(float(data[5]))
                        Delta_error.append(float(data[6]))
                    if pattern_depol.match(line):
                        depol_data = line.strip().split("\t")
                        Depolarization.append(float(depol_data[3]))
                        Depolarization_error.append(float(depol_data[4]))
            data_table = pd.DataFrame(columns = ["eV", "1/cm", "nm", "AOI", "Psi", "Delta", "Psi_error", "Delta_error", "Depolarization", "Depolarization_error", "pseudo_e1", "pseudo_e2"])
            if energy_unit == "eV":
                data_table["eV"] = Energy
                data_table["1/cm"] = [E_to_cm(E) for E in Energy]
                data_table["nm"] = [E_to_nm(E) for E in Energy]
            elif energy_unit == "1/cm":
                data_table["1/cm"] = Energy
                data_table["eV"] = [cm_to_E(cm) for cm in Energy]
                data_table["nm"] = [cm_to_nm(cm) for cm in Energy]
            elif energy_unit == "nm":
                data_table["nm"] = Energy
                data_table["eV"] = [nm_to_E(nm) for nm in Energy]
                data_table["1/cm"] = [nm_to_cm(nm) for nm in Energy]
            data_table["AOI"] = AOI
            data_table["Psi"] = Psi
            data_table["Delta"] = Delta
            data_table["Psi_error"] = Psi_error
            data_table["Delta_error"] = Delta_error
            if Depolarization:
                data_table["Depolarization"] = Depolarization
                data_table["Depolarization_error"] = Depolarization_error
            data_table["pseudo_e1"] = [pseudodielectric_function(psi, delta, aoi)[0] for psi, delta, aoi in zip(data_table["Psi"], data_table["Delta"], data_table["AOI"])]
            data_table["pseudo_e2"] = [pseudodielectric_function(psi, delta, aoi)[1] for psi, delta, aoi in zip(data_table["Psi"], data_table["Delta"], data_table["AOI"])]
            data_table = data_table.sort_values(by="eV", ignore_index=True)
            data_table
            measured_datas.append(data_table)
    # if len(measurement_temperatures) == 0:
    #     measurement_temperatures = [300] * len(measured_datas)
    return measured_datas, measurement_temperatures
