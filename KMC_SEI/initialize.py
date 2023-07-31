from support import *
from collection import *
###############################################################################
# working diroctory
dirc = os.getcwd()
###############################################################################
# reading the input file from xml file as yaml
with open('rendered_wano.yml') as file:
    wano_file = yaml.full_load(file)

xdim     = wano_file["Lattice size"]['Xdim']
ydim     = wano_file["Lattice size"]['Ydim']
T        = wano_file["Temperature (K)"]['T']
SaveStep = wano_file["Save step"]['SaveStep']
 #barriers E1 to E17
E1      = wano_file["Activation energy barriers (eV)"]['E1']
E2      = wano_file["Activation energy barriers (eV)"]['E2']
E3      = wano_file["Activation energy barriers (eV)"]['E3']
E4      = wano_file["Activation energy barriers (eV)"]['E4']
E5      = wano_file["Activation energy barriers (eV)"]['E5']
E6      = wano_file["Activation energy barriers (eV)"]['E6']
E7      = wano_file["Activation energy barriers (eV)"]['E7']
E8      = wano_file["Activation energy barriers (eV)"]['E8']
E9      = wano_file["Activation energy barriers (eV)"]['E9']
E10     = wano_file["Activation energy barriers (eV)"]['E10']
E11     = wano_file["Activation energy barriers (eV)"]['E11']
E12     = wano_file["Activation energy barriers (eV)"]['E12']
E13     = wano_file["Activation energy barriers (eV)"]['E13']
E13     = wano_file["Activation energy barriers (eV)"]['E13']
E14     = wano_file["Activation energy barriers (eV)"]['E14']
E15     = wano_file["Activation energy barriers (eV)"]['E15']   
E16     = wano_file["Activation energy barriers (eV)"]['E16']
E17     = wano_file["Activation energy barriers (eV)"]['E17']
E18     = wano_file["Activation energy barriers (eV)"]['E18']
E19     = wano_file["Activation energy barriers (eV)"]['E19']
"""
calling the lattice
"""
lattice = Lattice(xdim, ydim)
points = len(lattice)
###############################################################################
'''
lattice consideration: EC as the initial concentration of active EC.
 it should be taken into account as first distribution.
First lattice is like
    up    : A
    down  : B
    middle: S
'''
for i in lattice:
    if i.coordinate[1] == 0:
        i.species = "E"
    elif i.coordinate[1] == ydim+1:
        i.species = "A"
    else:
        i.species = "S"
        
'''
Firts possbile events are E+S-> E+P which means having the first electron reduction product close to electrodes.
To do this all E should be collected in a list.
'''
EE = [i for i in lattice if i.species == "E"]
top_list = EE
'''at beginging kmc needs to take all species could react
we need all S and clusters '''
# others = [i for i in lattice if i.species in ["S"]]
#
# top_list  = [j for j in lattice if j.status == 'sei' and "S" in [k.species for k in j.nbr]]
#creating an empty list for collecting events at each step.
pre = []
#set the number of steps (KMC steps)
num=0
t = 0
# event counter and residence time
'''to keep track of events and their consumed time, two lists will be created'''
# C = [0] * len(Events)
# res_time = [0] * len(Events)
'''a list for counting the number of outgoing species
0 -> P the first electron reduction product
1 -> O organic SEI ingredient molecule
2 -> O2 two molecules of organic SEI ingredient
3 -> C  cluster of organic SEI ingredients
'''
# Counter = [0] * 4
frac_save = SaveStep
# to save output
vis_ = 0
lists = pd.DataFrame()
# a counter to keep track of last saved step
keep = 0
tracking = {}
post_list = []
