#!/usr/bin/python

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import sys
path=sys.argv[1]
fileName=sys.argv[2]

# read medium-surf.vtk
mediumsurfvtk = LegacyVTKReader(FileNames=[path+'/'+fileName])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1100, 825]

# get display properties
mediumsurfvtkDisplay = GetDisplayProperties(mediumsurfvtk, view=renderView1)

# get color transfer function/color map for 'SkinSpacing'
skinSpacingLUT = GetColorTransferFunction('SkinSpacing')
skinSpacingLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 3.0578250885009766, 0.865003, 0.865003, 0.865003, 6.115650177001953, 0.705882, 0.0156863, 0.14902]
skinSpacingLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'SkinSpacing'
skinSpacingPWF = GetOpacityTransferFunction('SkinSpacing')
skinSpacingPWF.Points = [0.0, 0.0, 0.5, 0.0, 6.115650177001953, 1.0, 0.5, 0.0]
skinSpacingPWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
skinSpacingLUT.ApplyPreset('Blue to Red Rainbow', True)

# convert to log space
skinSpacingLUT.MapControlPointsToLogSpace()

# Properties modified on skinSpacingLUT
skinSpacingLUT.UseLogScale = 1

# Rescale transfer function
skinSpacingLUT.RescaleTransferFunction(0.05, 20.0)

# Properties modified on skinSpacingLUT
skinSpacingLUT.NumberOfTableValues = 9

# set active source
SetActiveSource(mediumsurfvtk)

# set scalar coloring
ColorBy(mediumsurfvtkDisplay, ('POINTS', 'SkinSpacing'))

# rescale color and/or opacity maps used to include current data range
mediumsurfvtkDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
mediumsurfvtkDisplay.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
skinSpacingLUT.RescaleTransferFunction(0.01, 5.0)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-347.5, 393.0, 41.0]
renderView1.CameraFocalPoint = [5377.5, -4260.0, -509.0]
renderView1.CameraViewUp = [-0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 1950


# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

scalarbar = GetScalarBar(skinSpacingLUT, renderView1)
scalarbar.Orientation = 'Horizontal'
scalarbar.Position = [0.05,0.85]
scalarbar.TitleColor = [1.0, 1.0, 1.0]
scalarbar.LabelColor = [1.0, 1.0, 1.0]
scalarbar.TitleFontSize = 9
scalarbar.LabelFontSize = 8


#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

SaveScreenshot(path+"/medium-nose.png")

from PyQt4 import QtCore
QtCore.QCoreApplication.quit()
