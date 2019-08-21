from mesa import Agent

class Endothelial(Agent):
    """" An Endotheial non-mobile agent imbedded in the grid-space of the model.

    variables:
    oxy -> oxygen level of the epithelial cell ranging from 0 to 100 (indicates the 'health' status)
    Oxy simulates the effect of downstream ischemia
    other celltypes can spawn on this part of the grid whenever oxy = 100

    """
    def __init__(self, unique_id,pos, oxy, model):
        super().__init__(unique_id, model)
        self.oxy = oxy
        self.pos = pos
        

        #if self.oxy <= 33:
         #   self.ec_attach = 100
          #  self.ec_rolling = 3
        #else:
         #   self.ec_attach = 0
          #  self.ec_rolling = 0

    #def activation(self):
     #   self.ec_

    #def step(self):
     #   self.activation()
        # if self.energy > 0:
        #   self.lose_energy()


class Neutrophil(Agent):
    """ An agent with fixed initial energy."""
    def __init__(self, unique_id, pos, model):
        super().__init__(unique_id, model)
        self.energy = 1
        self.pos = pos

    def migration(self):
        possible_steps = self.model.grid.get_neighborhood(self.pos, moore=True, include_center=False)
        neighbors = self.model.grid.get_neighbors(self.pos, 1, include_center=False)

        #only migration over the non-wounded areas.
        for agent in neighbors:
            if type(agent) is Endothelial:
                if agent.oxy < 49:
                    possible_steps.remove(agent.pos)

        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)
        for agent in self.pos:
        	if type(agent) is Endothelial:
        		if agent.oxy > 0:
        			agent.oxy = agent.oxy-10


    def burst(self):
        """
        This function simulates PMN respiratory burst and thus its cytotoxic/bactericidal effect.
        It does this primarily by updating a patch variable called "cytotox."
        In the presented model this represents overall free radical species (subsequent models will differentiate these species).
        It is updated at a value of 10 or "TNF," whichever is greater.  "Cytotox" has the following functions:
        """
        cellmates = self.model.grid.get_cell_list_contents([self.pos])
        if len(cellmates) > 1:
            other = self.random.choice(cellmates)
            other.energy += 1
            self.energy -= 1

    def update_cytokines(self):
        """" Updates Cytokine levels """

    def heal(self):
        """" heals EC """

    def apoptosis(self):
        "Neutrophil apoptosis"
        cellmates = self.model.grid.get_cell_list_contents([self.pos])
        if len(cellmates) > 1:
            other = self.random.choice(cellmates)
            other.energy += 1
            self.energy -= 1

    def step(self):
        self.migration()
        #if self.energy > 0:
         #   self.lose_energy()

class Macrophage(Agent):
    """ A Macrophage agent"""
    def __init__(self, unique_id,pos, model):
        super().__init__(unique_id, model)
        self.energy = 1
        self.pos = pos

    def move(self):
        possible_steps = self.model.grid.get_neighborhood(self.pos, moore=True, include_center=False)
        neighbors = self.model.grid.get_neighbors(self.pos, 1, include_center=False)

        #only migration over the non-wounded areas.
        for agent in neighbors:
            if type(agent) is Endothelial:
                if agent.oxy < 55:
                    possible_steps.remove(agent.pos)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)

    def step(self):
    	self.move()

class IL6(Agent):
    """ A IL6 (pro-inflamm) agent"""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.energy = 1

class IL10(Agent):
    """ A IL10 (anti-inflamm) agent"""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.energy = 1



