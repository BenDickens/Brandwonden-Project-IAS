'''
Wound healing ABM Prediction Model
================================
'''

from mesa.space import MultiGrid
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.datacollection import DataCollector

from agents import Neutrophil, Macrophage, Injury
from coordinates import *
from schedule import RandomActivationByAgent

class WoundModel(Model):
    """An ABM wound healing model simulating inflammation and contraction."""
    def __init__(self, Neutrophils, Macrophages,IL6, IL10, width, height, wound_radius):

        self.running = True
        self.neutrophils = Neutrophils
        self.macrophages = Macrophages
        self.IL6 = IL6
        self.IL10 = IL10
        self.grid = MultiGrid(width, height, True)
        self.schedule = RandomActivationByAgent(self)
        self.current_id = 0
        #create wound and non-wound region
        self.all_coordinates = all_coordinates(self.grid.width, self.grid.height)
        self.wound_radius = wound_radius
        #self.wound_coord = square_wound_area(self.grid.width, self.grid.height, self.wound_radius) # square shaped wound
        self.wound_coord = circle_wound_area(self.grid.width, self.grid.height, self.wound_radius, self.all_coordinates) #ellipse shaped wound
        self.non_wound_coord = []
        for coordinates in self.all_coordinates:
            if coordinates not in self.wound_coord:
                self.non_wound_coord += [coordinates]

        # Create Injury
        for i in range(len(self.wound_coord)):

            # Add the agent in the wound-region
            coord = self.wound_coord[i]
            x = coord[0]
            y = coord[1]
            a = Injury(self.next_id(),(x,y), self)
            self.schedule.add(a)
            self.grid.place_agent(a, (x, y))

        # Create Neutrophils
        for i in range(self.neutrophils):

            # Add the agent in a non-wound region
            coord = int(self.random.randrange(len(self.non_wound_coord)))
            coord = self.non_wound_coord[coord]
            x = coord[0]
            y = coord[1]
            a = Neutrophil(self.next_id(), (x,y), self)
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
        print('hi')
        #self.datacollector.collect(self)

    def run_model(self, step_count=200):

        for i in range(step_count):
            self.step()





