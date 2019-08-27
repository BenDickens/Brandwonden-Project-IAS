from mesa.visualization.modules import CanvasGrid
from mesa.visualization.ModularVisualization import ModularServer

from mesa.visualization.modules import ChartModule
from mesa.visualization.UserParam import UserSettableParameter

from ABM import WoundModel
from agents import *


Neutrophil_slider = UserSettableParameter('slider', "Number of Neutrophils", 24, 1, 50, 1)
Macrophage_slider = UserSettableParameter('slider', "Number of Macrophages", 19, 1, 50, 1)
wound_size_slider = UserSettableParameter('slider', 'Wound Radius',13,1,25,1)
coagulation_slider = UserSettableParameter('slider', 'Proportion of Coagulation', 0.7, 0, 1, 0.1)

def agent_portrayal(agent):
    portrayal = {"Shape": "circle",
                 "Filled": "true",
                 "r": 0.5}


    if type(agent) is Endothelial:
        if agent.TGFb > 0.7 and agent.wound:
            portrayal["Color"] = "purple"
            portrayal["Layer"] = 0
        else:
            if agent.oxy > 75:
                portrayal["Color"] = "tan"
                portrayal["Layer"] = 0
            elif agent.oxy <= 75 and agent.oxy >= 25:
                portrayal["Color"] = "orange"
                portrayal["Layer"] = 0
            elif agent.oxy < 25:
                portrayal["Color"] = "red"
                portrayal["Layer"] = 0



    elif type(agent) is Neutrophil:
        if agent.energy > 0:
            portrayal["Color"] = "green"
            portrayal["Layer"] = 1
        else:
            portrayal = {}

    elif type(agent) is Macrophage:
        if agent.energy > 0:
            if agent.phenotype == 0:
                portrayal["Color"] = "blue"
                portrayal["Layer"] = 1
            else:
                portrayal["Color"] = "cyan"
                portrayal["Layer"] = 1
        else:
            portrayal = {}

    return portrayal

grid = CanvasGrid(agent_portrayal, 25, 25, 500, 500)


server = ModularServer(WoundModel,
                       [grid],
                       "Burn Wound Healing Model",
                       {"Neutrophils": Neutrophil_slider, "Macrophages": Macrophage_slider, "width": 25, "height": 25, "wound_radius": wound_size_slider, "coagulation": coagulation_slider})