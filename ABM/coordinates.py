# determining the x,y coordinates of the whole model and of the wound area

import numpy as np

def square_wound_area(width, height, woundsize):
    # retrieve squared wound area coordinates
    wound_coord = []
    for x in range(int(((width - woundsize) / 2)), int((width - (width - woundsize) / 2))):
        for y in range(int(((height - woundsize) / 2)), int((height - (height - woundsize) / 2))):
            wound_coord += [(x, y)]
    return wound_coord


def circle_wound_area(width, height, radius, all_coordinates):
    # retrieve circular wound area coordinates
    a = radius # in case of ellipse -> a is different from b
    b = radius
    x0 = width / 2  # x center, half width
    y0 = height / 2  # y center, half height
    x = np.linspace(0, width, width)  # x values of interest
    y = np.linspace(0, height, height)[:, None]  # y values of interest, as a "column" array
    ellipse = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 <= 1; ellipse = np.ndarray.tolist(ellipse); ellipse = [item for sublist in ellipse for item in sublist]  # create True/False list for wound coordinates

    wound_coord = []
    for i in range(len(all_coordinates)):
        if ellipse[i] is True:
            wound_coord += [all_coordinates[i]]
    return wound_coord


def all_coordinates(width, height):
    all_coordinates = []
    for x in range(width):
        for y in range(height):
            all_coordinates += [(x, y)]
    return all_coordinates