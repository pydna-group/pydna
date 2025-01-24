# -*- coding: utf-8 -*-
from PIL import Image, ImageDraw, ImageFont, ImageFilter

# Notice the 'from PIL' at the start of the line
import numpy as np
import math

IMG_SIZE = 210, 600
LANES = np.zeros((2, IMG_SIZE[1]), dtype=np.int)
LANE_WIDTH = 80 * 2.0 / 3.0


def add_band(lane, bases, mass):
    lane_width = 80 * 2.0 / 3.0
    max_intensity = 256
    log = math.log(bases, 10)
    # logartigmic-exponential fitting
    peak_centre = int(-240.52 * log**2 + 717.16 * log + 224.7) * 2.0 / 3.0
    height = mass / (80 * 30.0 * log)
    max_spread = 10
    # handle primers
    if bases < 50:
        peak_centre += 50
        max_spread *= 4
        max_intensity /= 10
    band_spread = max_spread / log
    for i in range(max_spread, 0, -1):
        y1 = peak_centre - i
        y2 = peak_centre + i
        intensity = gaussian(y1, height, peak_centre, band_spread) * max_intensity
        add_to_lane(lane, y1, y2, intensity)


# a is peak height, b is peak centre and c is peak spread:
def gaussian(x, a, b, c):
    frac = float(((x - b) ** 2)) / (2 * (c**2))
    return a * math.exp(-frac)


def add_to_lane(lane, y1, y2, intensity):
    for y in range(int(y1), int(y2)):
        LANES[lane][y] += intensity


def normalise_intensity():
    global LANES
    for i, lane in enumerate(LANES):
        max_intensity = np.amax(LANES[i])
        if max_intensity > 256:
            LANES[i] = np.multiply(LANES[i], 256)
            LANES[i] = np.divide(LANES[i], max_intensity)


def draw_lanes(draw):
    normalise_intensity()
    for i, lane in enumerate(LANES):
        x1 = (75 + i * 130) * 2.0 / 3.0
        x2 = x1 + LANE_WIDTH + 2 * 2.0 / 3.0
        for y, intensity in enumerate(lane):
            y1 = y
            y2 = y + 1
            draw.rectangle((x1, y1, x2, y2), fill=(intensity, intensity, intensity))


def draw_label(draw, bases):
    log = math.log(bases, 10)
    peak_centre = int(-240.52 * log**2 + 717.16 * log + 220) * 2.0 / 3.0
    label = str(bases) + " -"
    label = label if len(label) >= 6 else str(" " * (6 - len(label))) + label
    # font = ImageFont.truetype("arial.ttf", 16)
    draw.text((5, peak_centre), str(label), (255, 255, 255))  # , font=font)


def draw_mock_band(draw):
    add_band(1, 1200, 35)
    add_band(1, 1000, 100)
    add_band(1, 830, 27)
    add_band(1, 800, 24)
    add_band(1, 780, 21)
    add_band(1, 330, 27)
    add_band(1, 300, 24)
    add_band(1, 270, 21)
    add_band(1, 10, 6)
    # add_band(1, 20, 100)


def draw_ladder(draw):
    f = 80
    add_band(0, 1517, 45 * f)
    add_band(0, 1200, 35 * f)
    add_band(0, 1000, 95 * f)
    add_band(0, 900, 27 * f)
    add_band(0, 800, 24 * f)
    add_band(0, 700, 21 * f)
    add_band(0, 600, 18 * f)
    add_band(0, 500, 97 * f)
    add_band(0, 400, 38 * f)
    add_band(0, 300, 29 * f)
    add_band(0, 200, 25 * f)
    add_band(0, 100, 48 * f)


def draw_labels(draw):
    draw_label(draw, 1517)
    draw_label(draw, 1200)
    draw_label(draw, 1000)
    draw_label(draw, 900)
    draw_label(draw, 800)
    draw_label(draw, 700)
    draw_label(draw, 600)
    draw_label(draw, 500)
    draw_label(draw, 400)
    draw_label(draw, 300)
    draw_label(draw, 200)
    draw_label(draw, 100)


def generate_gel(dna):
    im = Image.new("RGB", (IMG_SIZE), "#ddd")
    draw = ImageDraw.Draw(im)
    # Draw the background
    draw.rectangle((0, 0, IMG_SIZE), fill=(0, 0, 0))
    # Draw all LANES
    draw_ladder(draw)
    # draw_mock_band(draw)
    for bases, mass in dna:
        add_band(1, bases, mass)
    draw_lanes(draw)
    im = im.filter(ImageFilter.GaussianBlur(3))
    # Draw labels
    draw = ImageDraw.Draw(im)
    draw_labels(draw)
    return im


if __name__ == "__main__":
    length = 49
    mass = 700
    dna = [(length, mass)]
    im = generate_gel(dna)
    print(im)
    im.show()
    from IPython.display import display

    # display(im)

    # filename = "image_" + str(mass) + ".png"
    # im.save(filename)
