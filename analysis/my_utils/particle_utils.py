import numpy as np
SPEED_OF_LIGHT = 299.792458 # mm / ns

class Particle:
    def __init__(self, name, legend, color, mass, momentum, track_length):
        self.name = name
        self.legend = legend
        self.color = color
        self.mass = mass
        self.momentum = momentum
        self.track_length = track_length
        self.tof = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + (mass/momentum)**2 )
        self.beta = track_length/(SPEED_OF_LIGHT * self.tof)
        self.gamma = 1/np.sqrt(1 - self.beta**2)
