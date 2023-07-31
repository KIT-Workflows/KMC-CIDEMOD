"""
Required packages
"""
import itertools, os, random, sys, time, pickle
import numpy as np
import pandas as pd
from random import uniform, random, choice, sample, randrange,shuffle
from math import sqrt
from itertools import accumulate
import yaml
###############################################################################
'''                             FUNCTIONS                                   '''
###############################################################################
# Boltzman anf Planck constants
K_B = 8.617333262145e-5  # EVK-1
h = 4.135667696e-15  # EV.S
'''This function calculates the rate for included equations based on an
Arrhenius formula. It takes the activation energy barrier and temperature 
to calculate the rate'''
def bar_rate(e,T):
    return ((K_B * T)/h)*np.exp(-e/(K_B*T))
###############################################################################
'''This function handle the float number once they are going to be dumpped as a 
yanl file. Because the yaml file KMC dumps contains different types
of information, it is gonna be problematic if they are written directly into the yanl file'''
def round_float_list(float_list, decimal_points):
    float_list = [round(float(item),decimal_points) for item in float_list]
    return float_list
###############################################################################
# ---------------------------------SITE CLASS-----------------------------------
'''every site on lattice is an object of this class with 5 attributes for:

 coordinate: [x,y] pairs on lattice
 species   : name of the species on the site
 bonds     : if it is bonded species, then the inde of the other site
 nbr       : all 4 neighbors 
 status    : Keeping information about the occupation, formed SEI, etc.

'''
class site:

    def __init__(self, coordinate, species, bonds, nbr, status):
        self.coordinate = coordinate
        self.species = species
        self.nbr = nbr
        self.bonds = bonds
        self.status = status

    def setN(self, coordinate):
        self.coordinate = coordinate

    def setS(self, species):
        self.species = species

    def setnb(self, nbr):
        self.nbr = nbr

    def setb(self, bonds):
        self.bonds = bonds

    def setss(self, status):
        self.status = status

    def getC(self, coordinate):
        return self.coordinate

    def getS(self, species):
        return self.species

    def getnb(self, nbr):
        return self.nbr

    def getb(self, bonds):
        return self.bonds

    def getss(self, status):
        return self.status

    @staticmethod
    ###############################################################################
    # to species the color in visualization
    def colors(s):

        if Full(s) == "A":
            return 'grey'

        elif Full(s) == "E":
            return 'black'

        elif Full(s) == "P":
            return 'green'

        elif Full(s) == "S":
            return 'w'

        elif Full(s) == "O":
            return 'orange'

        elif Full(s) == "F":
            return 'red'

        elif Full(s) == "OO":
            return 'blue'

        elif Full(s) == "C":
            return 'm'
###############################################################################
# --------------------------------EVENT CLASS-----------------------------------
# events/reaction are objects of this class with four attributes
'''This class helps with creating event objects based on the main list of reactions.
This way, objects can be manipulated easily adding more details like the corresponding sites.
This can be handled using a subclass of this class.
Each event object can cary information about:

event reactants
event product
event barrier
event rate'''

class event:
    # Generic attributes
    def __init__(self ,reactant, product,barrier, rate):
        self.reactant = reactant
        self.rate     = rate
        self.product  = product
        self.barrier  = barrier
    # -------------------------------------------------------------------------------
    # assiging parameters in class

    def set_reactant(self, reactant):
        self.reactant = reactant

    def set_rate(self, rate):
        self.rate = rate

    def set_product(self, product):
        self.product = product

    def set_statt(self, barrier):
        self.barrier = barrier
 # ------------------------------------------------------------------------------
    # Get parameters

    def get_reactant(self, reactant):
        return self.reactant

    def get_rate(self, rate):
        return self.rate

    def get_product(self, product):
        return self.product

    def get_stat(self, stat):
        return self.stat
# -----------------------------NEW EVENT SUBCLASS-------------------------------
"""
this subclass helps in collecting events
"""
class new_event(event):
    def __init__(self, stat,reactant, product, rate, new_coord):
        super().__init__(stat,reactant, product, rate)
        self.new_coord = new_coord

    def set_newC(self, new_coord):
        self.new_coord = new_coord

    def get_newC(self, new_coord):
        return self.new_coord
# ----------------------------------- lattice-----------------------------------
"""
this function provides a NxN lattice with all grid points as site object which 
were defined before.
"""
def Lattice(xdim, ydim):
    # lattice points as a list of 2-elements sublists as [x,y]
    l = [[i, j] for i, j in itertools.product(range(1, xdim + 1), range(0, ydim + 2))]
    # initialize lattice as a list of site objects// coordinate,species,bonds,nbr,status
    sites = [site(i, [], [], [[i[0] + 1, i[1]], [i[0] - 1, i[1]], [i[0], i[1] + 1], [i[0], i[1] - 1]], [[], [], '', [0,'']])
             for i in l]
    # predefining the neighbors
    for i in sites:
        temp = i.nbr
        temp2 = [i for i in sites if i.coordinate in temp]
        i.nbr = temp2
    return sites
#------------------------------selection function----------------------------------
'''This is the function that selects the evens based on a search among all collected events(rates)'''
def ratee(i, r, R_index):
    if R_index[i - 1] < r <= R_index[i]:
        picked = R_index[i]
        # this index indicates the ith element in index list
        # e.g. if this one is 2 one should look at the index list
        # to find the true index for EVENT which is 3
        picked_index = i - 1
        return [picked, picked_index]
###############################################################################
'''This function returns a site corresponding species'''
def Full(x):
    if x.species in ["C","F"] or len(x.bonds) == 0:
        return x.species
    else:
        if x.bonds:
            temp = x.species
            for i in x.bonds:
                if i!=x:
                    temp = temp + i.species
            return temp
###############################################################################
'''As formed SEI components might be deposited at different region, this
function can return a list of sites for the SEI cluster'''
# cluster collector function
def cluster_x(single):
    listx = [single]
    stack = [j for j in single.bonds if j.species == "C"]
    while stack != []:
        picked = choice(stack)
        stack = [j for j in stack if j != picked]
        listx.append(picked)
        temp = [k for k in picked.bonds if k.species == "C" and k not in listx]
        stack.extend(temp)
    # to remove repetitive elements
    listx = list(dict.fromkeys(listx))
    return listx

#------------- a function to calculate thickness and porosity--------------------
def thickness(xdim, lattice):
    #thickness
    C = [i.coordinate for i in lattice if i.species ==  "C" and i.status[3][1] == 'SEI']
    F = [i.coordinate for i in lattice if i.species ==  "F"]
    sei_specs = C + F
    take_max_y = []
    for x in range(1,xdim+1):
        X = [i[1] for i in sei_specs if i[0] == x]
        if X:
            ymax = max(X)
            take_max_y.append(ymax)
            
    if take_max_y:
        thick = sum(take_max_y)/xdim
    else:
        thick = np.nan
        
    # porosity
    #creating lattice in proper size
    height = (thick)/2
    #fing the empty sites inside the SEI
    empties = [i for i in lattice if i.coordinate[1] < height and i.species not in ["C","F","E"]]
    porosity = len(empties)/(height*xdim)
    
    
    return [thick, porosity]

    
