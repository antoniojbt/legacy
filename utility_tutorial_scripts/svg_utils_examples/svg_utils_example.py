
# See: http://cairosvg.org/
# svg_utils: https://neuroscience.telenczuk.pl/?p=331
# https://github.com/btel/svg_utils

import svgutils.transform as sg
import sys
import cairosvg


# Create new SVG figure to place plots onto:
fig = sg.SVGFigure("16cm", "6.5cm")

# Load e.g. matpotlib-generated figures, can be any svg file though:
fig1 = sg.fromfile('sigmoid_fit.svg')
fig2 = sg.fromfile('anscombe.svg')

# Get the plot objects:
plot1 = fig1.getroot()
plot2 = fig2.getroot()
plot2.moveto(280, 0, scale=0.5)
#plot2.moveto(0, 280, scale = 1.0)


# Add text labels:
txt1 = sg.TextElement(25,20, "A", size=12, weight="bold")
txt2 = sg.TextElement(305,20, "B", size=12, weight="bold")

# Append plots and labels to figure:
fig.append([plot1, plot2])
fig.append([txt1, txt2])

# Save generated SVG files:
fig.save("fig_final.svg")

# Convert SVG images to other formats if needed:
cairosvg.svg2pdf(url='fig_final.svg', write_to='fig_final.pdf')

