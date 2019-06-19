from mesa.visualization.modules import CanvasGrid
from mesa.visualization.ModularVisualization import ModularServer

from mesa.visualization.modules import ChartModule
from mesa.visualization.UserParam import UserSettableParameter

from ABM import WoundModel
from agents import *


Neutrophil_slider = UserSettableParameter('slider', "Number of Neutrophils", 1000, 1, 5000, 1)
Macrophage_slider = UserSettableParameter('slider', "Number of Macrophages", 1000, 1, 5000, 1)
IL6_slider = UserSettableParameter('slider', "Number of IL6 cytokines", 1000, 1, 5000, 1)
IL10_slider = UserSettableParameter('slider', "Number of IL10 cytokines", 1000, 1, 5000, 1)
wound_size_slider = UserSettableParameter('slider', 'Wound Radius',25,1,50,1)

def agent_portrayal(agent):
    portrayal = {"Shape": "circle",
                 "Filled": "true",
                 "r": 0.5}
    if type(agent) is Injury:
        portrayal["Color"] = "red"
        portrayal["Layer"] = 0


    elif type(agent) is Neutrophil:
        if agent.energy > 0:
            portrayal["Color"] = "green"
            portrayal["Layer"] = 0
        else:
            portrayal = {}

    elif type(agent) is Macrophage:
        if agent.energy > 0:
            portrayal["Color"] = "blue"
            portrayal["Layer"] = 0
        else:
            portrayal = {}

    return portrayal

grid = CanvasGrid(agent_portrayal, 100, 100, 500, 500)


server = ModularServer(WoundModel,
                       [grid],
                       "Burn Wound healing model",
                       {"Neutrophils": Neutrophil_slider, "Macrophages": Macrophage_slider,"IL6": IL6_slider, "IL10": IL10_slider, "width": 100, "height": 100, "wound_radius": wound_size_slider})