# %%
import json, yaml, sys, os
from cideMOD import (
    ModelOptions,
    CSI
)
with open('rendered_wano.yml') as file:
    wano_file1 = yaml.full_load(file)

n_cycles     = wano_file1["cideMOD parameters"]['n_cycles']
crate        = wano_file1["cideMOD parameters"]['crate']
v_min        = wano_file1["cideMOD parameters"]['v_min']
v_max        = wano_file1["cideMOD parameters"]['v_max']
cycling_lump = wano_file1["Safari 2009 case files"]['Cycling lumped']
cycling_work = wano_file1["Safari 2009 case files"]['Cycling workflow']
OGV_G = wano_file1["Safari 2009 case files"]['OCV_G']
OGV_LCO = wano_file1["Safari 2009 case files"]['OCV_LCO']


# check if the ../kmc_data.yml exists
try:
    #open the file
    with open('../kmc_data.yml') as file:
        wano_file2 = yaml.full_load(file)
except IOError:
    print("I/O error: ../kmc_data.yml file does not exist")
    sys.exit()
    
  
sei_delta0        = wano_file2["thickness"]*1e-9
sei_porosity     = wano_file2["porosity"]
    

# %%
overwrite = False
case = "Safari_2009"
data_path = os.getcwd()

model_options = ModelOptions(mode='P2D', clean_on_exit=False, solve_SEI=True)

# %%
c_nom = 1.8

# %%
I_app = crate * c_nom
t_max = 1 / crate * 1.25
cycling = {
    "name": "Discharge Cycle",
            "type": "Cycle",
            "count": n_cycles,
            "steps": [
                {
                    "name": "Discharge",
                    "type": "Current",
                    "value": -I_app,
                    "unit": "A",
                    "t_max": {"value": t_max, "unit": "h"},
                    "store_delay": -1,
                    "min_step": 30,
                    "adaptive": True,
                    "events": [
                        {
                            "type": "Voltage",
                            "value": 2,
                            "unit": "V",
                            "atol": 1e-4,
                            "rtol": 1e-3,
                            "goto": "Next"
                        }
                    ]
                },
                {
                    "name": "Charge",
                    "type": "Current",
                    "value": I_app,
                    "unit": "A",
                    "t_max": {"value": t_max, "unit": "h"},
                    "store_delay": -1,
                    "min_step": 30,
                    "adaptive": True,
                    "events": [
                        {
                            "type": "Voltage",
                            "value": 4.2,
                            "unit": "V",
                            "atol": 1e-4,
                            "rtol": 1e-3,
                            "goto": "CV"
                        }
                    ]
                }
            ]}

cycling_test_plan = {
    'initial_state': {
        'SOC': 1,
        'exterior_temperature': 298
    },
    'steps': [cycling]
}


# %% Diffusion limited case
cell = json.load(open(cycling_lump))
sei_dict = cell['negativeElectrode']['SEI']
sei_dict['solventDiffusion'].update({"value": 6.8e-21})
sei_dict['rateConstant'].update({"value": 1.36e-7})

ref_porosity = sei_dict['solventPorosity']['value']
ref_density = sei_dict['density']['value']
ref_conductivity = sei_dict['conductivity']['value']

# THIS WOULD BE AN INPUT FROM KMC -----
sei_density = ref_density / (1 - ref_porosity) * (1 - sei_porosity)
# THIS WOULD BE AN INPUT FROM KMC -----
sei_conductivity = ref_conductivity / (1 - ref_porosity) * (1 - sei_porosity)

sei_dict['solventPorosity'].update({"value": sei_porosity})
sei_dict['density'].update({"value": sei_density})
sei_dict['delta0'].update({"value": sei_delta0})
sei_dict['conductivity'].update({"value": sei_conductivity})

bms_dl = CSI(cell, model_options, data_path, name=f'{case}_dl', overwrite=overwrite)
bms_dl.read_test_plan(cycling_test_plan)
bms_dl.run_test_plan()

problem_dl = bms_dl.problem
