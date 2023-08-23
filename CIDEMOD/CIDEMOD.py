# %%
import json, yaml, sys, os
import numpy as np
import matplotlib.pyplot as plt
from cideMOD import (
    ModelOptions,
    CSI
)
with open('rendered_wano.yml') as file:
    wano_file1 = yaml.full_load(file)

n_cycles     = wano_file1["Testplan parameters"]['n_cycles']
crate        = wano_file1["Testplan parameters"]['crate']
v_min        = wano_file1["Testplan parameters"]['v_min']
v_max        = wano_file1["Testplan parameters"]['v_max']
case         = wano_file1["cideMOD parameters"]['Case name']
params       = wano_file1["cideMOD parameters"]['Parameters json']
OCP_neg      = wano_file1["cideMOD parameters"]['Negative electrode OCP']
OCP_pos      = wano_file1["cideMOD parameters"]['Positive electrode OCP']
SEI_model    = wano_file1["cideMOD parameters"]['SEI model']


# check if the ../kmc_data.yml exists
try:
    #open the file
    with open('../kmc_data.yml') as file:
        wano_file2 = yaml.full_load(file)
except IOError:
    print("I/O error: ../kmc_data.yml file does not exist")
    sys.exit()
    
  
sei_delta0       = wano_file2["thickness"]*1e-9  #nm unit
sei_porosity     = wano_file2["porosity"]
    

# %%
overwrite = False
case = case.replace(" ", "_")
data_path = os.getcwd()

if SEI_model != 'Solvent diffusion-limited':
    print("Currently only 'Solvent diffusion-limited' SEI model is available in this version.")
    SEI_model = 'Solvent diffusion-limited'

model_options = ModelOptions(mode='P2D', clean_on_exit=False, solve_SEI=True)

# %% Modify SEI parameters and create case problem
cell = json.load(open(params))

# Note that the OCP filenames must match the ones in the json

# Define new SEI parameters for KMC
sei_dict = cell['negativeElectrode']['SEI']

ref_porosity = sei_dict['solventPorosity']['value']
ref_density = sei_dict['density']['value']
ref_conductivity = sei_dict['conductivity']['value']

sei_density = ref_density / (1 - ref_porosity) * (1 - sei_porosity)
sei_conductivity = ref_conductivity / (1 - ref_porosity) * (1 - sei_porosity)

sei_dict['solventPorosity'].update({"value": sei_porosity})
sei_dict['density'].update({"value": sei_density})
sei_dict['delta0'].update({"value": sei_delta0})
sei_dict['conductivity'].update({"value": sei_conductivity})

bms_dl = CSI(cell, model_options, data_path, name=case, overwrite=overwrite)

# %%
c_nom = bms_dl.problem.Q

# %% Define test plan
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


# %% Run case
bms_dl.read_test_plan(cycling_test_plan)
bms_dl.run_test_plan()

problem_dl = bms_dl.problem

# %% Plot figures

# Extract problem variables
vars_ = ['time', 'current', 'voltage', 'capacity', 'delta_sei_a0', 'Q_sei_a0']
results = {var_: np.array(problem_dl.WH.get_global_variable(var_)) for var_ in vars_}
save_path = problem_dl.WH.save_path

# SEI properties
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), dpi=150)

index = np.append(np.where(np.diff(results['current'])<(-crate*c_nom))[0], 
                  len(results['current'])-1)
# index2 = np.append(np.where(np.diff(results['current'])>(crate*c_nom*1.2))[0], 
#                   len(results['current'])-1)
nCycles = len(index)
cycles = np.arange(1,nCycles+1)

av_thickness = results['delta_sei_a0'][index]
ax1.plot(cycles, av_thickness, "-.", label='cideMOD')
ax1.set_title("SEI thickness")
ax1.set_xlabel("Cycle nr.")
ax1.set_ylabel("SEI thickness [$\mu\,\mathrm{m}$]")
# ax1.legend(loc="best")
ax1.set_xlim([0, nCycles])

av_capacity_loss = results['Q_sei_a0'][index]
ax2.plot(cycles, av_capacity_loss/c_nom*100, "-.", label='cideMOD')
ax2.set_title("Capacity loss")
ax2.set_xlabel("Cycle nr.")
ax2.set_ylabel("Capacity loss [$\%$]")
# ax1.legend(loc="best")
ax2.set_xlim([0, nCycles])

plt.tight_layout()
plt.savefig(os.path.join(save_path,"SEI_results.png"))

# First and last cycle discharge curve
fig, ax1 = plt.subplots(1, 1, figsize=(6, 4), dpi=150)
dch_voltage = results['voltage'][results['current'] < 0]
dch_capacity = results['capacity'][results['current'] < 0]
index2 = np.where(np.diff(dch_voltage)>0)[0]
ax1.plot(dch_capacity[0:index2[0]], dch_voltage[0:index2[0]], 
         label="Cycle 1")
ax1.plot(dch_capacity[index2[-2]+1:index2[-1]]-dch_capacity[index2[-2]+1], 
         dch_voltage[index2[-2]+1:index2[-1]], 
         label=f"Cycle {n_cycles}")

ax1.legend(loc="best")
ax1.set_xlabel("Discharged capacity [Ah]")
ax1.set_ylabel("Voltage [V]")

plt.savefig(os.path.join(save_path,"voltage_vs_capacity.png"))
