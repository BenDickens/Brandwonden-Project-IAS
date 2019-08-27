'''
Wound healing ABM Prediction Model
================================
'''

from mesa.space import MultiGrid
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.datacollection import DataCollector

from agents import Neutrophil, Macrophage,  Endothelial
from coordinates import *
from schedule import RandomActivationByAgent

class WoundModel(Model):
    """An ABM wound healing model simulating inflammation and contraction."""
    def __init__(self, Neutrophils, Macrophages, width, height, wound_radius, coagulation):

        self.running = True
        self.neutrophils = Neutrophils
        self.macrophages = Macrophages
        self.grid = MultiGrid(width, height, True)
        self.schedule = RandomActivationByAgent(self)
        self.current_id = 0
        self.centre = (width//2, height//2)
        self.coagulation_size = wound_radius * coagulation 

        #create wound and non-wound region
        self.all_coordinates = all_coordinates(self.grid.width, self.grid.height)
        self.wound_radius = wound_radius
        #self.wound_coord = square_wound_area(self.grid.width, self.grid.height, self.wound_radius) # square shaped wound
        self.wound_coord = circle_wound_area(self.grid.width, self.grid.height, self.wound_radius, self.all_coordinates) #ellipse shaped wound
        self.coagulation_coord = circle_wound_area(self.grid.width, self.grid.height, self.coagulation_size, self.all_coordinates)
        self.non_wound_coord = []
        for coordinates in self.all_coordinates:
            if coordinates not in self.wound_coord:
                self.non_wound_coord += [coordinates]

        # Create non-wound Endothelial cells
        for i in range(len(self.non_wound_coord)):
            # Add the agent in the non wound-region
            coord = self.non_wound_coord[i]
            a = Endothelial(self.next_id(),(coord[0],coord[1]),100,0, self)
            self.schedule.add(a)
            self.grid.place_agent(a, (coord[0], coord[1]))

        # Create wound Endothelial-cells
        for i in range(len(self.wound_coord)):
            # Add the agent in the wound-region
            coord = self.wound_coord[i]
            if(coord in self.coagulation_coord):
                a = Endothelial(self.next_id(),(coord[0],coord[1]),0,1, self)
            else:
            	a = Endothelial(self.next_id(),(coord[0],coord[1]),50,1, self)
            self.schedule.add(a)
            self.grid.place_agent(a, (coord[0], coord[1]))

        # Create Neutrophils
        for i in range(self.neutrophils):

            # Add the agent in a non-wound region
            coord = int(self.random.randrange(len(self.non_wound_coord)))
            coord = self.non_wound_coord[coord]
            x = coord[0]
            y = coord[1]
            a = Neutrophil(self.next_id(), (x,y),self.centre, self)
            self.schedule.add(a)
            self.grid.place_agent(a, (x, y))

        # Create Macrophages
        for i in range(self.macrophages):

            # Add the agent in a non-wound region
            coord = int(self.random.randrange(len(self.non_wound_coord)))
            coord = self.non_wound_coord[coord]
            x = coord[0]
            y = coord[1]
            a = Macrophage(self.next_id(),(x,y), self)
            self.schedule.add(a)
            self.grid.place_agent(a, (x, y))

        self.running = True
        #self.datacollector.collect(self)

    def step(self):
        self.schedule.step()
        print('step')
        #self.datacollector.collect(self)

    def run_model(self, step_count=200):

        for i in range(step_count):
            self.step()





