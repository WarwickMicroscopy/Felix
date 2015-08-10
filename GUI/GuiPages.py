#!/usr/bin/python
# -*- coding: utf-8 -*-

import wx
import wx.lib.agw.floatspin as FS
import math
import os
import shutil
import sys
import GuiPages
import Bin2Tiff
import FileCtrl


class FlagPanel(wx.Panel):

  def __init__(self, parent):

    wx.Panel.__init__(self, parent)

    title = wx.StaticText(self, wx.ID_ANY, 'Flags')

#=========================================================================

    flag1 = {'name': 'IWriteFLAG', 'choices': [
        'Silent', 'Crucial information', 'Basic information', 'All information'], 'object type': 'CHOICE', 'default': 'All information'}
    flag2 = {'name': 'IScatterFactorMethodFLAG', 'choices': [
        'Kirkland', 'Doyle-Turner', 'Peng', 'Lobato'], 'object type': 'CHOICE', 'default': 'Kirkland'}
    flag3 = {'name': 'IMaskFLAG', 'choices': [],
             'object type': 'CHECKBOX', 'default': 0}
    flag4 = {'name': 'IZolzFLAG', 'choices': [],
             'object type': 'CHECKBOX', 'default': 0}
    flag5 = {'name': 'IAbsorbFLAG', 'choices': [
        'None', 'Proportional'], 'object type': 'CHOICE', 'default': 'Proportional'}
    flag6 = {'name': 'IAnisoDebyeWallerFLAG', 'choices': [
        '0'], 'object type': 'CHOICE', 'default': '0'}
    flag7 = {'name': 'IPseudoCubicFLAG', 'choices': [
        '0'], 'object type': 'CHOICE', 'default': '0'}
    flag8 = {'name': 'IXDirectionFLAG', 'choices': [
        'Automatic', 'Manual'], 'object type': 'CHOICE', 'default': 'Automatic'}

    flags = []
    self.flags = flags

    flags.append(flag1)
    flags.append(flag2)
    flags.append(flag3)
    flags.append(flag4)
    flags.append(flag5)
    flags.append(flag6)
    flags.append(flag7)
    flags.append(flag8)

    # check
    print "Checking flags\n"
    for flag in flags:
      print 'Checking flag: {0}.\n'.format(flag['name'])
      if type(flag['name']) != str:
        sys.exit("Incorrect value for name in flag (please use a string)\n")
      if type(flag['choices']) != list:
        sys.exit(
            "Incorrect value for the choices in flag (please use a list of strings, or an empty list for a checkbox)\n")
      for choice in flag['choices']:
        if flag['choices'] == False or type(choice) != str:
          sys.exit("Incorrect value for choice in flag (please use a string)\n")
      if flag['object type'] != 'CHOICE' and flag['object type'] != 'CHECKBOX':
        sys.exit(
            "Incorrect value for object type in flag (use CHECKBOX or CHOICE\n")
      if flag['object type'] == 'CHOICE':
        if type(flag['default']) != str:
          sys.exit(
              "Incorrect value for default in flag (use a string for CHOICE types\n")
      if flag['object type'] == 'CHECKBOX':
        if flag['default'] != 0 and flag['default'] != 1:
          sys.exit(
              "Incorrect value for default in flag (use a 0 or 1 for checkbox off or on, respectively\n")


#=========================================================================

    # print flagnames

    # Making the code a bit more future proof by generating flag layout
    # from a list of flag names
    flagnumber = len(flags)
    print 'The number of flags: {0}.\n'.format(flagnumber)
    numberOfRows = int(math.ceil(flagnumber / 3.0))
    print 'The number of rows: {0}.\n'.format(numberOfRows)

    # Make some lists
    flagObjectsLabels = []
    flagObjectsChoices = []
    SizerObjects = []

    self.flagObjectsChoices = flagObjectsChoices

    # Finds the number of empty slots on the bottom row to add spacer later
    spacerNo = (3 - (flagnumber % 3)) % 3
    print 'The number of spacers: {0}.\n'.format(spacerNo)

    # Adds a sizer for each row to a sizer list
    for x in range(0, numberOfRows):
      SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))

    # Adds labels and choice objects to each respective lists
    for flag in flags:
      choices = flag['choices']
      print(choices)
      flagObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, flag['name']))
      flagIndex = flags.index(flag)
      if flag['object type'] == 'CHOICE':
        flagObjectsChoices.append(wx.Choice(self, wx.ID_ANY, size=(100, -1),
                                            choices=choices, name=flag['name']))
        flagObjectsChoices[flagIndex].SetStringSelection(flag['default'])

      elif flag['object type'] == 'CHECKBOX':
        flagObjectsChoices.append(wx.CheckBox(self, wx.ID_ANY, size=(100, -1),
                                              name=flag['name']))
        if flag['default'] == 1:
          flagObjectsChoices[flagIndex].SetValue(True)

    # Adds the objects from the lists to the respective rows in the sizer list
    for flagNo in range(0, flagnumber):
      row = int(math.floor(flagNo / 3))
      SizerObjects[row].Add(flagObjectsLabels[flagNo], 2, wx.ALL, 5)
      SizerObjects[row].Add(flagObjectsChoices[flagNo], 1, wx.ALL, 5)

    # Adds spacers if necessary
    if spacerNo != 0:
      for x in range(0, spacerNo):
        SizerObjects[numberOfRows - 1].AddStretchSpacer(3)

    # Set up overall and title sizers
    topflagSizer = wx.BoxSizer(wx.VERTICAL)
    flagTitleSizer = wx.BoxSizer(wx.HORIZONTAL)

    flagTitleSizer.Add(title, 0, wx.ALL, 5)

    topflagSizer.Add(flagTitleSizer, 0, wx.CENTER)

    # Add sizers from list to topflagSizer
    for sizerNo in range(0, numberOfRows):
      topflagSizer.Add(SizerObjects[sizerNo], 0, wx.CENTER)

    self.SetSizer(topflagSizer)
    topflagSizer.Fit(self)


class RadiusPanel(wx.Panel):

  def __init__(self, parent):

    wx.Panel.__init__(self, parent)

    title = wx.StaticText(self, wx.ID_ANY, 'Radius of Beam')
    RadiusLabel = wx.StaticText(self, wx.ID_ANY, 'Radius of Beam in Pixels')

    RadiusSizer = wx.BoxSizer(wx.HORIZONTAL)
    RadiusTitleSizer = wx.BoxSizer(wx.HORIZONTAL)
    topRadiusSizer = wx.BoxSizer(wx.VERTICAL)

    self.IPixelCount = FS.FloatSpin(self, size=(60, -1), value=64, min_val=0, max_val=512,
                                    increment=64, agwStyle=FS.FS_RIGHT)
    self.IPixelCount.SetDigits(0)
    RadiusSizer.Add(RadiusLabel, 3, wx.ALL, 5)
    RadiusSizer.Add(self.IPixelCount, 1, wx.ALL, 5)
    RadiusSizer.AddStretchSpacer(4)

    RadiusTitleSizer.Add(title, 0, wx.ALL, 5)

    topRadiusSizer.Add(RadiusTitleSizer, 0, wx.CENTER)
    topRadiusSizer.Add(RadiusSizer, 0)

    self.SetSizer(topRadiusSizer)
    topRadiusSizer.Fit(self)


class BeamPanel(wx.Panel):

  def __init__(self, parent):

    wx.Panel.__init__(self, parent)

    BeamControlList = []
    self.BeamControlList = BeamControlList

#=========================================================================

    # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
    # ['name', default value, increment, min, max, 'type']
    # with type referring to a 1 for spinctrl or 2 for a float spin!
    # NB: spin does not need increment, so just put 0!
    self.BeamControl1 = ['IMinReflectionPool', '15', 0, 0, 100000, 1]
    self.BeamControl2 = ['IMinStrongBeams', '7', 0, 0, 100000, 1]
    self.BeamControl3 = ['IMinWeakBeams', '5', 0, 0, 100000, 1]
    self.BeamControl4 = ['RBSBMax', '0.1', 0.1, 0, 100000, 2]
    self.BeamControl5 = ['RBSPMax', '0.1', 0.1, 0, 100000, 2]

    # Add them to list (of lists) - Need to find a better method for this
    BeamControlList.append(self.BeamControl1)
    BeamControlList.append(self.BeamControl2)
    BeamControlList.append(self.BeamControl3)
    BeamControlList.append(self.BeamControl4)
    BeamControlList.append(self.BeamControl5)

    # check
    print "Checking beam controls\n"
    for BeamCtrl in BeamControlList:
      print 'Checking BeamCtrl: {0}.\n'.format(BeamCtrl[0])
      if type(BeamCtrl[0]) != str:
        sys.exit("Incorrect value for name in BeamCtrl (please use a string)\n")
      if type(BeamCtrl[5]) != int:
        sys.exit("Incorrect value for type in BeamCtrl (please use a 1 or a 2)\n")
      if BeamCtrl[5] == 1:
        if type(BeamCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in BeamCtrl (please use a string)\n")
        if type(BeamCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in BeamCtrl (please use an int)\n")
        if type(BeamCtrl[3]) != int:
          sys.exit("Incorrect value for min in BeamCtrl (please use a int)\n")
        if type(BeamCtrl[4]) != int:
          sys.exit("Incorrect value for max in BeamCtrl (please use a int)\n")
      if BeamCtrl[5] == 2:
        if type(BeamCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in BeamCtrl (please use a string)\n")
        if type(BeamCtrl[2]) != float and type(BeamCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in BeamCtrl (please use an float or int)\n")
        if type(BeamCtrl[3]) != float and type(BeamCtrl[3]) != int:
          sys.exit(
              "Incorrect value for min in BeamCtrl (please use a float or int)\n")
        if type(BeamCtrl[4]) != float and type(BeamCtrl[4]) != int:
          sys.exit(
              "Incorrect value for max in BeamCtrl (please use a float or int)\n")

#=========================================================================

    title = wx.StaticText(self, wx.ID_ANY, 'Beam Selection')

    # Making the code a bit more future proof by generating beam layout
    # from a list of beam names
    beamnumber = len(BeamControlList)
    print 'The number of beam controls: {0}.\n'.format(beamnumber)
    numberOfRows = int(math.ceil(beamnumber / 3.0))
    print 'The number of rows: {0}.\n'.format(numberOfRows)

    # Make some lists
    SizerObjects = []
    beamObjectsLabels = []
    beamObjectsControls = []

    self.SizerObjects = SizerObjects
    self.beamObjectsLabels = beamObjectsLabels
    self.beamObjectsControls = beamObjectsControls

    # Finds the number of empty slots on the bottom row to add spacer later
    spacerNo = (3 - (beamnumber % 3)) % 3
    print 'The number of spacers: {0}.\n'.format(spacerNo)

    # Adds a sizer for each row to a sizer list
    for x in range(0, numberOfRows):
      SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))

    # Adds labels and choice objects to each respective lists
    for beam in BeamControlList:
      beamname = beam[0]
      beamvalue = beam[1]
      beamincrement = beam[2]
      beammin = beam[3]
      beammax = beam[4]
      beamtype = beam[5]

      currentIndex = BeamControlList.index(beam)

      beamObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, beamname))

      if beamtype == 1:
        beamObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY, size=(60, -1),
                                               value=beamvalue, min=beammin, max=beammax))
        print('Added a spin!')
      elif beamtype == 2:
        beamObjectsControls.append(FS.FloatSpin(self, size=(60, -1),
                                                value=beamvalue, increment=beamincrement, min_val=beammin,
                                                max_val=beammax, agwStyle=FS.FS_RIGHT))

        beamObjectsControls[currentIndex].SetFormat("%f")
        beamObjectsControls[currentIndex].SetDigits(1)
        print('Added a float spin!')

    # Adds the objects from the lists to the respective rows in the sizer list
    for beamNo in range(0, beamnumber):
      row = int(math.floor(beamNo / 3))
      SizerObjects[row].Add(beamObjectsLabels[beamNo], 3, wx.ALL, 5)
      SizerObjects[row].Add(beamObjectsControls[beamNo], 1, wx.ALL, 5)

    # Adds spacers if necessary
    if spacerNo != 0:
      for x in range(0, spacerNo):
        SizerObjects[numberOfRows - 1].AddStretchSpacer(4)

    # Set up overall and title sizers
    topbeamSizer = wx.BoxSizer(wx.VERTICAL)
    beamTitleSizer = wx.BoxSizer(wx.HORIZONTAL)

    beamTitleSizer.Add(title, 0, wx.ALL, 5)

    topbeamSizer.Add(beamTitleSizer, 0, wx.CENTER)

    # Add sizers from list to topbeamSizer
    for sizerNo in range(0, numberOfRows):
      topbeamSizer.Add(SizerObjects[sizerNo], 0)

    self.SetSizer(topbeamSizer)
    topbeamSizer.Fit(self)


class crystalPanel(wx.Panel):

  def __init__(self, parent):

    wx.Panel.__init__(self, parent)

    crystalControlList = []

#=========================================================================

    # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
    # ['name', default value, increment, min, mac, 'type', *NUMBER OF DIGITS*]
    # with type referring to a 1 for spinctrl or 2 for a float spin!
    # NB: spin does not need increment, so just put 0!
    crystalControl1 = ['RDebyeWallerConstant', '0.467', 0.001, 0, 100000, 2, 4]
    crystalControl2 = ['RAbsorptionPer', '2.9', 0.1, 0, 100000, 2, 1]

    # Add them to list (of lists) - Need to find a better method for this
    crystalControlList.append(crystalControl1)
    crystalControlList.append(crystalControl2)

    # check
    print "Checking crystal controls\n"
    for CrystalCtrl in crystalControlList:
      print 'Checking CrystalCtrl: {0}.\n'.format(CrystalCtrl[0])
      if type(CrystalCtrl[0]) != str:
        sys.exit("Incorrect value for name in CrystalCtrl (please use a string)\n")
      if type(CrystalCtrl[5]) != int:
        sys.exit("Incorrect value for type in CrystalCtrl (please use a 1 or a 2)\n")
      if CrystalCtrl[5] == 1:
        if type(CrystalCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in CrystalCtrl (please use a string)\n")
        if type(CrystalCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in CrystalCtrl (please use an int)\n")
        if type(CrystalCtrl[3]) != int:
          sys.exit("Incorrect value for min in CrystalCtrl (please use a int)\n")
        if type(CrystalCtrl[4]) != int:
          sys.exit("Incorrect value for max in CrystalCtrl (please use a int)\n")
      if CrystalCtrl[5] == 2:
        if type(CrystalCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in CrystalCtrl (please use a string)\n")
        if type(CrystalCtrl[2]) != float and type(CrystalCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in CrystalCtrl (please use an float or int)\n")
        if type(CrystalCtrl[3]) != float and type(CrystalCtrl[3]) != int:
          sys.exit(
              "Incorrect value for min in CrystalCtrl (please use a float or int)\n")
        if type(CrystalCtrl[4]) != float and type(CrystalCtrl[4]) != int:
          sys.exit(
              "Incorrect value for max in CrystalCtrl (please use a float or int)\n")
        if type(CrystalCtrl[6]) != int:
          sys.exit(
              "Incorrect value for no. of digits in CrystalCtrl (please use an int)\n")
#=========================================================================

    title = wx.StaticText(self, wx.ID_ANY, 'Crystal Settings')

    # Making the code a bit more future proof by generating crystal layout
    # from a list of crystal names
    crystalnumber = len(crystalControlList)
    print 'The number of crystal controls: {0}.\n'.format(crystalnumber)
    numberOfRows = int(math.ceil(crystalnumber / 3.0))
    print 'The number of rows: {0}.\n'.format(numberOfRows)

    # Make some lists
    SizerObjects = []
    self.crystalObjectsLabels = []
    self.crystalObjectsControls = []

    # Finds the number of empty slots on the bottom row to add spacer later
    spacerNo = (3 - (crystalnumber % 3)) % 3
    print 'The number of spacers: {0}.\n'.format(spacerNo)

    # Adds a sizer for each row to a sizer list
    for x in range(0, numberOfRows):
      SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))

    # Adds labels and choice objects to each respective lists
    for crystal in crystalControlList:
      crystalname = crystal[0]
      crystalvalue = crystal[1]
      crystalincrement = crystal[2]
      crystalmin = crystal[3]
      crystalmax = crystal[4]
      crystaltype = crystal[5]
      crystaldigits = crystal[6]

      currentIndex = crystalControlList.index(crystal)

      self.crystalObjectsLabels.append(
          wx.StaticText(self, wx.ID_ANY, crystalname))

      if crystaltype == 1:
        self.crystalObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY,
                                                       size=(60, -1), value=crystalvalue, min=crystalmin, max=crystalmax))
        print('Added a spin!')
      elif crystaltype == 2:
        self.crystalObjectsControls.append(FS.FloatSpin(self, size=(60, -1),
                                                        value=crystalvalue, increment=crystalincrement,
                                                        min_val=crystalmin, max_val=crystalmax, agwStyle=FS.FS_RIGHT))

        self.crystalObjectsControls[currentIndex].SetFormat("%f")
        self.crystalObjectsControls[currentIndex].SetDigits(crystaldigits)
        print('Added a float spin!')

    # Adds the objects from the lists to the respective rows in the sizer list
    for crystalNo in range(0, crystalnumber):
      row = int(math.floor(crystalNo / 3))
      SizerObjects[row].Add(self.crystalObjectsLabels[crystalNo], 3, wx.ALL, 5)
      SizerObjects[row].Add(self.crystalObjectsControls[
                            crystalNo], 1, wx.ALL, 5)

    # Adds spacers if necessary
    if spacerNo != 0:
      for x in range(0, spacerNo):
        SizerObjects[numberOfRows - 1].AddStretchSpacer(4)

    # Set up overall and title sizers
    topcrystalSizer = wx.BoxSizer(wx.VERTICAL)
    crystalTitleSizer = wx.BoxSizer(wx.HORIZONTAL)

    crystalTitleSizer.Add(title, 0, wx.ALL, 5)

    topcrystalSizer.Add(crystalTitleSizer, 0, wx.CENTER)

    # Add sizers from list to topcrystalSizer
    for sizerNo in range(0, numberOfRows):
      topcrystalSizer.Add(SizerObjects[sizerNo], 0)

    self.SetSizer(topcrystalSizer)
    topcrystalSizer.Fit(self)


class microscopePanel(wx.Panel):

  def __init__(self, parent):

    wx.Panel.__init__(self, parent)

    microscopeControlList = []

#=========================================================================

    # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
    # ['name', default value, increment, min, mac, 'type']
    # with type referring to a 1 for spinctrl or 2 for a float spin!
    # NB: spin does not need increment, so just put 0!
    microscopeControl1 = ['ROuterConvergenceAngle', '3.0', 0.1, 0, 50, 2]
    microscopeControl2 = ['RInnerConvergenceAngle', '0.0', 0.1, 0, 50, 2]
    microscopeControl3 = [
        'IIncidentBeamDirectionX', '1', 0, -100000, 100000, 1]
    microscopeControl4 = [
        'IIncidentBeamDirectionY', '1', 0, -100000, 100000, 1]
    microscopeControl5 = [
        'IIncidentBeamDirectionZ', '1', 0, -100000, 100000, 1]
    microscopeControl6 = ['IXDirectionX', '1', 0, -100000, 100000, 1]
    microscopeControl7 = ['IXDirectionY', '1', 0, -100000, 100000, 1]
    microscopeControl8 = ['IXDirectionZ', '1', 0, -100000, 100000, 1]
    microscopeControl9 = ['INormalDirectionX', '1', 0, -100000, 100000, 1]
    microscopeControl10 = ['INormalDirectionY', '1', 0, -100000, 100000, 1]
    microscopeControl11 = ['INormalDirectionZ', '1', 0, -100000, 100000, 1]
    microscopeControl12 = ['RAcceleratingVoltage', '200.0', 0.1, 0, 100000, 2]
    microscopeControl13 = ['RAcceptanceAngle', '0.0', 0.1, 0, 180, 2]

    # Add them to list (of lists) - Need to find a better method for this
    microscopeControlList.append(microscopeControl1)
    microscopeControlList.append(microscopeControl2)
    microscopeControlList.append(microscopeControl3)
    microscopeControlList.append(microscopeControl4)
    microscopeControlList.append(microscopeControl5)
    microscopeControlList.append(microscopeControl6)
    microscopeControlList.append(microscopeControl7)
    microscopeControlList.append(microscopeControl8)
    microscopeControlList.append(microscopeControl9)
    microscopeControlList.append(microscopeControl10)
    microscopeControlList.append(microscopeControl11)
    microscopeControlList.append(microscopeControl12)
    microscopeControlList.append(microscopeControl13)

    # check
    print "Checking microscope controls\n"
    for MicroscopeCtrl in microscopeControlList:
      print 'Checking MicroscopeCtrl: {0}.\n'.format(MicroscopeCtrl[0])
      if type(MicroscopeCtrl[0]) != str:
        sys.exit(
            "Incorrect value for name in MicroscopeCtrl (please use a string)\n")
      if type(MicroscopeCtrl[5]) != int:
        sys.exit(
            "Incorrect value for type in MicroscopeCtrl (please use a 1 or a 2)\n")
      if MicroscopeCtrl[5] == 1:
        if type(MicroscopeCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in MicroscopeCtrl (please use a string)\n")
        if type(MicroscopeCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in MicroscopeCtrl (please use an int)\n")
        if type(MicroscopeCtrl[3]) != int:
          sys.exit("Incorrect value for min in MicroscopeCtrl (please use a int)\n")
        if type(MicroscopeCtrl[4]) != int:
          sys.exit("Incorrect value for max in MicroscopeCtrl (please use a int)\n")
      if MicroscopeCtrl[5] == 2:
        if type(MicroscopeCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in MicroscopeCtrl (please use a string)\n")
        if type(MicroscopeCtrl[2]) != float and type(MicroscopeCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in MicroscopeCtrl (please use an float or int)\n")
        if type(MicroscopeCtrl[3]) != float and type(MicroscopeCtrl[3]) != int:
          sys.exit(
              "Incorrect value for min in MicroscopeCtrl (please use a float or int)\n")
        if type(MicroscopeCtrl[4]) != float and type(MicroscopeCtrl[4]) != int:
          sys.exit(
              "Incorrect value for max in MicroscopeCtrl (please use a float or int)\n")

#=========================================================================

    title = wx.StaticText(self, wx.ID_ANY, 'Microscope Selection')

    # Making the code a bit more future proof by generating microscope layout
    # from a list of microscope names
    microscopenumber = len(microscopeControlList)
    print 'The number of microscope controls: {0}.\n'.format(microscopenumber)
    numberOfRows = int(math.ceil(microscopenumber / 3.0))
    print 'The number of rows: {0}.\n'.format(numberOfRows)

    # Make some lists
    SizerObjects = []
    self.microscopeObjectsLabels = []
    self.microscopeObjectsControls = []

    # Finds the number of empty slots on the bottom row to add spacer later
    spacerNo = (3 - (microscopenumber % 3)) % 3
    print 'The number of spacers: {0}.\n'.format(spacerNo)

    # Adds a sizer for each row to a sizer list
    for x in range(0, numberOfRows):
      SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))

    # Adds labels and choice objects to each respective lists
    for microscope in microscopeControlList:
      microscopename = microscope[0]
      microscopevalue = microscope[1]
      microscopeincrement = microscope[2]
      microscopemin = microscope[3]
      microscopemax = microscope[4]
      microscopetype = microscope[5]

      currentIndex = microscopeControlList.index(microscope)

      self.microscopeObjectsLabels.append(
          wx.StaticText(self, wx.ID_ANY, microscopename))

      if microscopetype == 1:
        self.microscopeObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY,
                                                          size=(60, -1), value=microscopevalue,
                                                          min=microscopemin, max=microscopemax))
        print('Added a spin!')
      elif microscopetype == 2:
        self.microscopeObjectsControls.append(FS.FloatSpin(self, size=(60, -1),
                                                           value=microscopevalue, increment=microscopeincrement,
                                                           min_val=microscopemin, max_val=microscopemax, agwStyle=FS.FS_RIGHT))

        self.microscopeObjectsControls[currentIndex].SetFormat("%f")
        self.microscopeObjectsControls[currentIndex].SetDigits(1)
        print('Added a float spin!')

    # Adds the objects from the lists to the respective rows in the sizer list
    for microscopeNo in range(0, microscopenumber):
      row = int(math.floor(microscopeNo / 3))
      SizerObjects[row].Add(self.microscopeObjectsLabels[
                            microscopeNo], 3, wx.ALL, 5)
      SizerObjects[row].Add(self.microscopeObjectsControls[
                            microscopeNo], 1, wx.ALL, 5)

    # Adds spacers if necessary
    if spacerNo != 0:
      for x in range(0, spacerNo):
        SizerObjects[numberOfRows - 1].AddStretchSpacer(4)

    # Set up overall and title sizers
    topmicroscopeSizer = wx.BoxSizer(wx.VERTICAL)
    microscopeTitleSizer = wx.BoxSizer(wx.HORIZONTAL)

    microscopeTitleSizer.Add(title, 0, wx.ALL, 5)

    topmicroscopeSizer.Add(microscopeTitleSizer, 0, wx.CENTER)

    # Add sizers from list to topmicroscopeSizer
    for sizerNo in range(0, numberOfRows):
      topmicroscopeSizer.Add(SizerObjects[sizerNo], 0)

    self.SetSizer(topmicroscopeSizer)
    topmicroscopeSizer.Fit(self)


class imagePanel(wx.Panel):

  def __init__(self, parent):

    wx.Panel.__init__(self, parent)

    imageControlList = []

#=========================================================================

    # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
    # ['name', default value, increment, min, mac, 'type', NUMBER OF DIGITS]
    # with type referring to a 1 for spinctrl or 2 for a float spin, 3 for checkbox!
    # NB: spin does not need increment, so just put 0!
    imageControl1 = ['RInitialThickness', '1000.0', 1, 0, 100000, 2, 1]
    imageControl2 = ['RFinalThickness', '1000.0', 1, 0, 100000, 2, 1]
    imageControl3 = ['RDeltaThickness', '10.0', 1, 0, 100000, 2, 1]
    imageControl4 = ['IReflectOut', '7', 0, 0, 100000, 1, 1]
    IImageFLAG1 = ['Montage', True, 0, 0, 0, 3, 0]
    IImageFLAG2 = ['Stack Reflections', False, 0, 0, 0, 3, 0]
    IImageFLAG3 = ['Amplitude and Phase', False, 0, 0, 0, 3, 0]

    # Add them to list (of lists) - Need to find a better method for this
    imageControlList.append(imageControl1)
    imageControlList.append(imageControl2)
    imageControlList.append(imageControl3)
    imageControlList.append(imageControl4)
    imageControlList.append(IImageFLAG1)
    imageControlList.append(IImageFLAG2)
    imageControlList.append(IImageFLAG3)

    # check
    print "Checking image controls\n"
    for ImageCtrl in imageControlList:
      print 'Checking ImageCtrl: {0}.\n'.format(ImageCtrl[0])
      if type(ImageCtrl[0]) != str:
        sys.exit("Incorrect value for name in ImageCtrl (please use a string)\n")
      if type(ImageCtrl[5]) != int:
        sys.exit("Incorrect value for type in ImageCtrl (please use a 1 or a 2)\n")
      if ImageCtrl[5] == 1:
        if type(ImageCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in ImageCtrl (please use a string)\n")
        if type(ImageCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in ImageCtrl (please use an int)\n")
        if type(ImageCtrl[3]) != int:
          sys.exit("Incorrect value for min in ImageCtrl (please use a int)\n")
        if type(ImageCtrl[4]) != int:
          sys.exit("Incorrect value for max in ImageCtrl (please use a int)\n")
      if ImageCtrl[5] == 2:
        if type(ImageCtrl[1]) != str:
          sys.exit(
              "Incorrect value for default value in ImageCtrl (please use a string)\n")
        if type(ImageCtrl[2]) != float and type(ImageCtrl[2]) != int:
          sys.exit(
              "Incorrect value for increment in ImageCtrl (please use an float or int)\n")
        if type(ImageCtrl[3]) != float and type(ImageCtrl[3]) != int:
          sys.exit(
              "Incorrect value for min in ImageCtrl (please use a float or int)\n")
        if type(ImageCtrl[4]) != float and type(ImageCtrl[4]) != int:
          sys.exit(
              "Incorrect value for max in ImageCtrl (please use a float or int)\n")
        if type(ImageCtrl[6]) != int:
          sys.exit(
              "Incorrect value for no. of digits in ImageCtrl (please use an int)\n")

#=========================================================================

    title = wx.StaticText(self, wx.ID_ANY, 'image Settings')

    # Making the code a bit more future proof by generating image layout
    # from a list of image names
    imagenumber = len(imageControlList)
    print 'The number of image controls: {0}.\n'.format(imagenumber)
    numberOfRows = int(math.ceil(imagenumber / 3.0))
    print 'The number of rows: {0}.\n'.format(numberOfRows)

    # Make some lists
    SizerObjects = []
    self.imageObjectsLabels = []
    self.imageObjectsControls = []

    # Finds the number of empty slots on the bottom row to add spacer later
    spacerNo = (3 - (imagenumber % 3)) % 3
    print 'The number of spacers: {0}.\n'.format(spacerNo)

    # Adds a sizer for each row to a sizer list
    for x in range(0, numberOfRows):
      SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))

    # Adds labels and choice objects to each respective lists
    for image in imageControlList:
      imagename = image[0]
      imagevalue = image[1]
      imageincrement = image[2]
      imagemin = image[3]
      imagemax = image[4]
      imagetype = image[5]
      imagedigits = image[6]

      currentIndex = imageControlList.index(image)

      self.imageObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, imagename))

      if imagetype == 1:
        self.imageObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY,
                                                     size=(60, -1), value=imagevalue, min=imagemin, max=imagemax))
        print('Added a spin!')
      elif imagetype == 2:
        self.imageObjectsControls.append(FS.FloatSpin(self, size=(60, -1),
                                                      value=imagevalue, increment=imageincrement, min_val=imagemin,
                                                      max_val=imagemax, agwStyle=FS.FS_RIGHT))

        self.imageObjectsControls[currentIndex].SetFormat("%f")
        self.imageObjectsControls[currentIndex].SetDigits(imagedigits)
        print('Added a float spin!')
      elif imagetype == 3:
        self.imageObjectsControls.append(wx.CheckBox(self, wx.ID_ANY, size=(60, -1),
                                                     name=imagename))
        self.imageObjectsControls[currentIndex].SetValue(imagevalue)

    # Adds the objects from the lists to the respective rows in the sizer list
    for imageNo in range(0, imagenumber):
      row = int(math.floor(imageNo / 3))
      SizerObjects[row].Add(self.imageObjectsLabels[imageNo], 3, wx.ALL, 5)
      SizerObjects[row].Add(self.imageObjectsControls[imageNo], 1, wx.ALL, 5)

    # Adds spacers if necessary
    if spacerNo != 0:
      for x in range(0, spacerNo):
        SizerObjects[numberOfRows - 1].AddStretchSpacer(4)

    # Set up overall and title sizers
    topimageSizer = wx.BoxSizer(wx.VERTICAL)
    imageTitleSizer = wx.BoxSizer(wx.HORIZONTAL)

    imageTitleSizer.Add(title, 0, wx.ALL, 5)

    topimageSizer.Add(imageTitleSizer, 0, wx.CENTER)

    # Add sizers from list to topimageSizer
    for sizerNo in range(0, numberOfRows):
      topimageSizer.Add(SizerObjects[sizerNo], 0)

    self.SetSizer(topimageSizer)
    topimageSizer.Fit(self)


class optionPanel(wx.Panel):

  def __init__(self, parent, main):
    wx.Panel.__init__(self, parent)

    # Probably should change this at some point
    self.main = main

    title = wx.StaticText(self, wx.ID_ANY, 'Options')

    # optionSizer            = wx.BoxSizer(wx.HORIZONTAL)
    optionTitleSizer = wx.BoxSizer(wx.HORIZONTAL)
    topOptionSizer = wx.BoxSizer(wx.VERTICAL)

    optionTitleSizer.Add(title, 0, wx.ALL, 5)

    # Text box
    self.CIFtextbox = wx.TextCtrl(
        self, -1, os.getcwd(), size=(400, -1), style=wx.TE_RIGHT)
    textSizer = wx.BoxSizer(wx.HORIZONTAL)
    textSizer.Add(self.CIFtextbox, 4, wx.ALL, 5)
    # textSizer.AddStretchSpacer(4)

    # Buttons
    Run = wx.Button(self, label='Run')
    Cancel = wx.Button(self, label='Cancel')
    CIFFile = wx.Button(self, label='Browse')
    InputFile = wx.Button(self, label='Save Input File')
    InputLoad = wx.Button(self, label='Load Input File')

    buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
    buttonSizer.Add(Run, 0, wx.ALL, 5)
    buttonSizer.Add(Cancel, 0, wx.ALL, 5)
    buttonSizer.Add(CIFFile, 0, wx.ALL, 5)
    buttonSizer.Add(InputFile, 0, wx.ALL, 5)
    buttonSizer.Add(InputLoad, 0, wx.ALL, 5)

    # CSC Checkbox
    CSC = wx.CheckBox(self, label='CSC')
    CSC.SetValue(True)
    checkSizer = wx.BoxSizer(wx.HORIZONTAL)
    checkSizer.Add(CSC, 0, wx.ALL, 5)

    # Number of cores
    coreLabel = wx.StaticText(self, label='MpiCores')
    self.MPICores = wx.SpinCtrl(self, size=(60, -1), value='1', min=0, max=100)
    coreSizer = wx.BoxSizer(wx.HORIZONTAL)
    coreSizer.Add(coreLabel, 3, wx.ALL, 5)
    coreSizer.Add(self.MPICores, 1, wx.ALL, 5)
    # coreSizer.AddStretchSpacer(4)

    topOptionSizer.Add(optionTitleSizer, 0, wx.CENTER)
    topOptionSizer.Add(textSizer, 0, wx.CENTER)
    topOptionSizer.Add(buttonSizer, 0, wx.CENTER)
    topOptionSizer.Add(checkSizer, 0, wx.CENTER)
    topOptionSizer.Add(coreSizer, 0, wx.CENTER)

    # topOptionSizer.Add(optionSizer, 0)
    self.SetSizer(topOptionSizer)
    topOptionSizer.Fit(self)

    #InputFileInst = FileCtrl.InputFileCreateObj(self)
    # the various functions of the buttons - run felix, write input file,
    # cancel, and browse file
    Run.Bind(wx.EVT_BUTTON, self.CIFCreateWrapper)
    Cancel.Bind(wx.EVT_BUTTON, self.OnCloseWrapper)
    CIFFile.Bind(wx.EVT_BUTTON, self.OnCifWrapper)
    InputFile.Bind(wx.EVT_BUTTON, self.InpCreateWrapper)
    InputLoad.Bind(wx.EVT_BUTTON, self.LoadInputFileWrapper)

  def CIFCreateWrapper(self, event):
    FileCtrl.CIFCreate(self, event)

  def OnCloseWrapper(self, event):
    FileCtrl.OnClose(self, event)

  def OnCifWrapper(self, event):
    FileCtrl.OnCif(self, event)

  def InpCreateWrapper(self, event):
    FileCtrl.InpCreate(self, event)

  def LoadInputFileWrapper(self, event):
    FileCtrl.LoadInputFile(self, event)

class ViewerPanel(wx.Panel):

  def __init__(self, parent):
    wx.Panel.__init__(self, parent)
    self.PhotoMaxSize = 240
    self.createWidgets()

  def createWidgets(self):
    instructions = 'Browse for an image'
    img = wx.EmptyImage()
    self.imageCtrl = wx.StaticBitmap(self, wx.ID_ANY,
                                     wx.BitmapFromImage(img))

    instructLbl = wx.StaticText(self, label=instructions)
    self.photoTxt = wx.TextCtrl(self, size=(200, -1))
    browseBtn = wx.Button(self, label='Browse')
    browseBtn.Bind(wx.EVT_BUTTON, self.onBrowse)

    self.mainSizer = wx.BoxSizer(wx.VERTICAL)
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.mainSizer.Add(wx.StaticLine(self, wx.ID_ANY),
                       0, wx.ALL | wx.EXPAND, 5)
    self.mainSizer.Add(instructLbl, 0, wx.ALL, 5)
    self.sizer.Add(self.photoTxt, 0, wx.ALL, 5)
    self.sizer.Add(browseBtn, 0, wx.ALL, 5)
    self.mainSizer.Add(self.sizer, 0, wx.ALL, 5)
    self.mainSizer.Add(self.imageCtrl, 0, wx.ALL, 5)

    self.SetSizer(self.mainSizer)

    self.Layout()

  def onBrowse(self, event):
    """
    Browse for file
    """
    wildcard = "TIFF files (*.tif)|*.tif"
    dialog = wx.FileDialog(None, "Choose a file",
                           wildcard=wildcard,
                           style=wx.OPEN)
    if dialog.ShowModal() == wx.ID_OK:
      self.photoTxt.SetValue(dialog.GetPath())
    dialog.Destroy()
    self.onView()


  def onView(self):
    filepath = self.photoTxt.GetValue()
    img = wx.Image(filepath, wx.BITMAP_TYPE_ANY)
    # scale the image, preserving the aspect ratio
    W = img.GetWidth()
    H = img.GetHeight()
    if W > H:
        NewW = self.PhotoMaxSize
        NewH = self.PhotoMaxSize * H / W
    else:
        NewH = self.PhotoMaxSize
        NewW = self.PhotoMaxSize * W / H
    img = img.Scale(NewW,NewH)

    self.imageCtrl.SetBitmap(wx.BitmapFromImage(img))
    self.Refresh()
