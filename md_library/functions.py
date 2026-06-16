import matplotlib as mpl
import numpy as np
from typing import Any, Optional, Union
import os

def ColorFade(c1: Any, c2: Any, mix: float, temperature: Union[int, float], T_c: Optional[Union[int, float]]) -> str:
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    if T_c is not None and temperature == int(T_c):
        return "black"
    else:
        return mpl.colors.to_hex((1-mix) * c1 + mix * c2)

def GetLabel(temperature: Union[int, float], T_c: Optional[Union[int, float]]) -> str:
    if T_c is not None and temperature == int(T_c):
        return r"$T_\mathrm{C} = $" + str(temperature)
    else:
        return str(temperature)

def load_file_names(sample_name: str, file_names_list: list, specification: Optional[str] = None, suffix: Optional[str] = None) -> list:
    file_names = []
    for file_name in file_names_list:
        if sample_name in file_name and (specification is None or specification in file_name) and (suffix is None or file_name.endswith(suffix)):
            file_names.append(file_name)
    print(f"\033[34m Found \033[0m {len(file_names)} \033[34m files for sample \033[0m '{sample_name}':")
    print(file_names)
    return file_names

def load_fit_report(file_name, file_location, sample_name):
    fit_report_file_name = f"{file_name}_fit_report.txt"
    fit_report_file_location = os.path.join(file_location, f"{file_name}_fit_report.txt")
    if os.path.exists(fit_report_file_location):
        fit_report_data = open(fit_report_file_location, 'r').read()
        print(f"\033[34m Found file \033[0m {fit_report_file_name} \033[34m with fit report for sample \033[0m '{sample_name}'")
        return fit_report_data, fit_report_file_name, fit_report_file_location
    else:
        print(f"\033[31m No file \033[0m {fit_report_file_name} \033[31m with fit report for sample \033[0m '{sample_name}' \033[31m found. Fit report will be created after fitting. \033[0m")
        return None, fit_report_file_name, fit_report_file_location




    

