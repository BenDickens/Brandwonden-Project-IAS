from mesa import Agent

class Injury(Agent):
    """ An agent with fixed initial energy."""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)


class Neutrophil(Agent):
    """ An agent with fixed initial energy."""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.energy = 1

    def move(self):
        possible_steps = self.model.grid.get_neighborhood(self.pos,include_center=False)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)

    #activated or not?
    def lose_energy(self):
        cellmates = self.model.grid.get_cell_list_contents([self.pos])
        if len(cellmates) > 1:
            other = self.random.choice(cellmates)
            other.energy += 1
            self.energy -= 1

    def step(self):
        self.move()
        #if self.energy > 0:
         #   self.lose_energy()

class Macrophage(Agent):
    """ A Macrophage agent"""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.energy = 1

    def move(self):
        possible_steps = self.model.grid.get_neighborhood(
            self.pos,
            include_center=False)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)

class IL6(Agent):
    """ A IL6 (pro-inflamm) agent"""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.energy = 1

class IL10(Agent):
    """ A IL6 (anti-inflamm) agent"""
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.energy = 1


