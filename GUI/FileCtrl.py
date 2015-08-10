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


# Class defining gui which selects where the user wants to save the input file
# When Felix is run, this file will be copied as felix.inp


class WriteInputDialog(wx.Dialog):

  def __init__(self, parent, id, title):
    wx.Dialog.__init__(self, parent, id, title, size=(400, 200))
    self.UserInterface2()

  def UserInterface2(self):

    wx.StaticText(self, label="Choose location for Input file",
                  pos=(5, 20), style=wx.ALIGN_CENTRE_HORIZONTAL)
    self.InputFileTextBox = wx.TextCtrl(
        self, -1, os.getcwd(), pos=(5, 100), size=(250, -1), style=wx.TE_RIGHT)
    InputFileBrowse = wx.Button(
        self, label='Save', pos=(270, 100), size=(60, -1))
    InputFileBrowse.Bind(wx.EVT_BUTTON, self.OnINPBrowse)

    ReturnButton = wx.Button(self, label='OK', pos=(200, 130), size=(60, -1))
    CancelButton = wx.Button(
        self, label='Cancel', pos=(300, 130), size=(60, -1))

    ReturnButton.Bind(wx.EVT_BUTTON, self.ReturnSuccess)
    CancelButton.Bind(wx.EVT_BUTTON, self.ReturnCancel)

  def OnINPBrowse(self, event):
    InputFileDialog = wx.FileDialog(self, "Save Input file", "", "",
                                    "Inp Files (*.inp)|*.inp", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

    if InputFileDialog.ShowModal() == wx.ID_CANCEL:
      # IErr=-1
      return()

    UserInputPath = InputFileDialog.GetPath()

    self.InputFileTextBox.SetValue(UserInputPath)

  def ReturnSuccess(self, event):

    self.SInpPath, self.SInpFilename = os.path.split(
        self.InputFileTextBox.GetValue())

    if self.SInpFilename.find(".inp") == -1:
      wx.MessageBox('.inp file not found, please save a .inp file',
                    'Error', wx.OK | wx.ICON_ERROR)
      return()

    self.Close()

  def ReturnCancel(self, event):

    self.Close()




# subroutine which re-names the selected cif file and opens a directory
# for the input files - then runs Felix
def CIFCreate(event):

  # opens new (re-existing) directory and puts all input files into it
  #***ensure file and path check needs to be put here****
  cpath, cfilename = os.path.split(self.CIFtextbox.GetValue())

  if cfilename.find(".cif") == -1:
    wx.MessageBox('.cif file not found, please select a .cif file', 'Error',
                  wx.OK | wx.ICON_ERROR)

    # creates working directory
  cfilename = cfilename.replace(" ", "")
  print cfilename
  dir = cfilename.rstrip('.cif') + "_" + str(self.main.notebook.page3.beamObjectsControls[0].GetValue())\
      + "_" + \
      str(self.main.notebook.page3.beamObjectsControls[
          1].GetValue()) + "_input_directory"
  if os.path.exists(dir):
    shutil.rmtree(dir)
  os.makedirs(dir)

  # copy cif file from user specified location to new directory
  shutil.copy2(self.CIFtextbox.GetValue(), dir + "/felix.cif")

  # copy sca file from samples folder - TO BE CHANGED, only for Richard's Use
  shutil.copy2("../samples/Si/felix.sca", dir + "/felix.sca")

  InputFileSwitch = 1

  # Write input file
  self.WriteInputFile(dir, InputFileSwitch)
  wx.MessageBox('Files created successfully, click okay to run Felix', 'Info',
                wx.OK | wx.CANCEL | wx.ICON_INFORMATION)

  # Change to working directory
  os.chdir(dir)

  # Get Value of mpicores
  NumberofCores = self.MPICores.GetValue()

  # Run in parallel or single core
  if NumberofCores == 1:
    os.system("../felixsim")  # single core
  else:
    os.system("mpirun -n " + str(NumberofCores) + " ../felixsim")  # parallel

  self.Close(True)



#TEST CLASSSSSSSS_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DOESN'T~~~WORK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class InputFileCreateObj(GuiPages.optionPanel):

  def __init__(self, parent, main):
    wx.Panel.__init__(self, parent)

  def InpCreate(self, event):

    # if os.path.exists(dir):
    #  shutil.rmtree(dir)
    # os.makedirs(dir)

    InputFileCreate = WriteInputDialog(self, -1, 'Save File As')
    InputFileBlock = InputFileCreate.ShowModal()

    dir = InputFileCreate.SInpPath

    # InputFileCreate.Destroy()

    InputFileSwitch = 2

    self.WriteInputFile(dir, InputFileSwitch)

    x.MessageBox('Input File successfully written',
                'Info', wx.OK | wx.ICON_INFORMATION)
    # wx.MessageBox('Input File successfully written',
    #              'Info', wx.OK | wx.ICON_INFORMATION)

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

def WriteInputFile(dir, InputFileSwitch):

  FelixInpFilename = dir + "/felix.inp"

  print FelixInpFilename

  inpfile = open(FelixInpFilename, "wb")

  # opening text in input file
  inpfile.write(
      "# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:" + "\n")
  inpfile.write("# ------------------------------------\n")
  inpfile.write("\n")
  inpfile.write("# ------------------------------------\n")
  inpfile.write("# felixsim input\n")
  inpfile.write("\n")
  inpfile.write("# control flags\n")

  # FLAGS=================================================================

  # IWriteFLAG
  IWriteFLAGVal = self.main.notebook.page1.flagObjectsChoices[
      0].GetCurrentSelection()
  inpfile.write("IWriteFLAG                = ")
  inpfile.write(str(IWriteFLAGVal))
  inpfile.write("\n")

  # IImageFLAG
  inpfile.write("IImageFLAG                = ")

  IImageFLAGVal0 = self.main.notebook.page6.imageObjectsControls[
      4].GetValue()
  if IImageFLAGVal0 == True:
    inpfile.write("0")

  IImageFLAGVal1 = self.main.notebook.page6.imageObjectsControls[
      5].GetValue()
  if IImageFLAGVal1 == True:
    inpfile.write("1")

  IImageFLAGVal2 = self.main.notebook.page6.imageObjectsControls[
        6].GetValue()
  if IImageFLAGVal2 == True:
    inpfile.write("2")

  inpfile.write("\n")

  # IScatterFactorMethodFLAG
  IScatterFactorMethodFLAGVal = self.main.notebook.page1.flagObjectsChoices[
      1].GetCurrentSelection()
  inpfile.write("IScatterFactorMethodFLAG  = ")
  inpfile.write(str(IScatterFactorMethodFLAGVal))
  inpfile.write("\n")

  # IMaskFLAG
  IMaskFLAGVal = self.main.notebook.page1.flagObjectsChoices[2].GetValue()
  inpfile.write("IMaskFLAG                 = ")
  if IMaskFLAGVal == True:
    IMaskFLAGBool = 1
  elif IMaskFLAGVal == False:
    IMaskFLAGBool = 0
  inpfile.write(str(IMaskFLAGBool))
  inpfile.write("\n")

  # IZolzFLAG
  IZolzFLAGVal = self.main.notebook.page1.flagObjectsChoices[3].GetValue()
  inpfile.write("IZolzFLAG                 = ")
  if IZolzFLAGVal == True:
    IZolzFLAGBool = 1
  elif IZolzFLAGVal == False:
    IZolzFLAGBool = 0
  inpfile.write(str(IZolzFLAGBool))
  inpfile.write("\n")

  # IAbsorbFLAG
  IAbsorbFLAGVal = self.main.notebook.page1.flagObjectsChoices[
      4].GetCurrentSelection()
  inpfile.write("IAbsorbFLAG               = ")
  inpfile.write(str(IAbsorbFLAGVal))
  inpfile.write("\n")

  # IAnisoDebyeWallerFLAG
  IAnisoDebyeWallerFLAGVal = self.main.notebook.page1.flagObjectsChoices[
      5].GetCurrentSelection()
  inpfile.write("IAnisoDebyeWallerFLAG     = ")
  inpfile.write(str(IAnisoDebyeWallerFLAGVal))
  inpfile.write("\n")

  # IPseudoCubicFLAG
  IPseudoCubicFLAGVal = self.main.notebook.page1.flagObjectsChoices[
      6].GetCurrentSelection()
  inpfile.write("IPseudoCubicFLAG          = ")
  inpfile.write(str(IPseudoCubicFLAGVal))
  inpfile.write("\n")

  # IXDirectionFLAG
  IXDirectionFLAGVal = self.main.notebook.page1.flagObjectsChoices[
      7].GetCurrentSelection()
  inpfile.write("IXDirectionFLAG           = ")
  inpfile.write(str(IXDirectionFLAGVal))
  inpfile.write("\n")

  # RADIUS OF BEAM=========================================================
  inpfile.write("\n# radius of the beam in pixels\n")

  # IPixelCount
  IPixelCountVal = int(self.main.notebook.page2.IPixelCount.GetValue())
  inpfile.write("IPixelCount               = ")
  inpfile.write(str(IPixelCountVal))
  inpfile.write("\n")

  # BEAM SELECTION=========================================================
  inpfile.write("\n# beam selection criteria\n")

  # IMinReflectionPool
  IMinReflectionPoolVal = self.main.notebook.page3.beamObjectsControls[
      0].GetValue()
  inpfile.write("IMinReflectionPool        = ")
  inpfile.write(str(IMinReflectionPoolVal))
  inpfile.write("\n")

  # IMinStrongBeams
  IMinStrongBeamsVal = self.main.notebook.page3.beamObjectsControls[
      1].GetValue()
  inpfile.write("IMinStrongBeams           = ")
  inpfile.write(str(IMinStrongBeamsVal))
  inpfile.write("\n")

  # IMinWeakBeams
  IMinWeakBeamsVal = self.main.notebook.page3.beamObjectsControls[
      2].GetValue()
  inpfile.write("IMinWeakBeams             = ")
  inpfile.write(str(IMinWeakBeamsVal))
  inpfile.write("\n")

  # RBSBMax
  RBSBMaxVal = self.main.notebook.page3.beamObjectsControls[3].GetValue()
  inpfile.write("RBSBMax                   = ")
  inpfile.write(str(RBSBMaxVal))
  inpfile.write("\n")

  # RBSPMax
  RBSPMaxVal = self.main.notebook.page3.beamObjectsControls[4].GetValue()
  inpfile.write("RBSPMax                   = ")
  inpfile.write(str(RBSPMaxVal))
  inpfile.write("\n")

  # CRYSTAL SETTINGS======================================================
  inpfile.write("\n# crystal settings\n")

  # RDebyeWallerConstant
  RDebyeWallerConstantVal = self.main.notebook.page4.crystalObjectsControls[
      0].GetValue()
  inpfile.write("RDebyeWallerConstant      = ")
  inpfile.write(str(RDebyeWallerConstantVal))
  inpfile.write("\n")

  # RAbsorptionPer
  RAbsorptionPerVal = self.main.notebook.page4.crystalObjectsControls[
      1].GetValue()
  inpfile.write("RAbsorptionPer            = ")
  inpfile.write(str(RAbsorptionPerVal))
  inpfile.write("\n")

  # MICROSCOPE SETTINGS====================================================
  inpfile.write("\n# microscope settings\n")

  # ROuterConvergenceAngle
  ROuterConvergenceAngleVal = self.main.notebook.page5.microscopeObjectsControls[
      0].GetValue()
  inpfile.write("ROuterConvergenceAngle    = ")
  inpfile.write(str(ROuterConvergenceAngleVal))
  inpfile.write("\n")

  # ROuterConvergenceAngle
  RInnerConvergenceAngleVal = self.main.notebook.page5.microscopeObjectsControls[
      1].GetValue()
  inpfile.write("RInnerConvergenceAngle    = ")
  inpfile.write(str(RInnerConvergenceAngleVal))
  inpfile.write("\n")

  # IIncidentBeamDirectionX
  IIncidentBeamDirectionXVal = self.main.notebook.page5.microscopeObjectsControls[
      2].GetValue()
  inpfile.write("IIncidentBeamDirection    = [")
  inpfile.write(str(IIncidentBeamDirectionXVal))

  # IIncidentBeamDirectionY
  IIncidentBeamDirectionYVal = self.main.notebook.page5.microscopeObjectsControls[
      3].GetValue()
  inpfile.write(",")
  inpfile.write(str(IIncidentBeamDirectionYVal))

  # IIncidentBeamDirectionZ
  IIncidentBeamDirectionZVal = self.main.notebook.page5.microscopeObjectsControls[
      4].GetValue()
  inpfile.write(",")
  inpfile.write(str(IIncidentBeamDirectionZVal))
  inpfile.write("]\n")

  # IXDirectionX
  IXDirectionXVal = self.main.notebook.page5.microscopeObjectsControls[
      5].GetValue()
  inpfile.write("IXDirection               = [")
  inpfile.write(str(IXDirectionXVal))

  # IXDirectionY
  IXDirectionYVal = self.main.notebook.page5.microscopeObjectsControls[
      6].GetValue()
  inpfile.write(",")
  inpfile.write(str(IXDirectionYVal))

  # IXDirectionZ
  IXDirectionZVal = self.main.notebook.page5.microscopeObjectsControls[
      7].GetValue()
  inpfile.write(",")
  inpfile.write(str(IXDirectionZVal))
  inpfile.write("]\n")

  # INormalDirectionX
  INormalDirectionXVal = self.main.notebook.page5.microscopeObjectsControls[
      8].GetValue()
  inpfile.write("INormalDirectionX         = [")
  inpfile.write(str(INormalDirectionXVal))

  # INormalDirectionY
  INormalDirectionYVal = self.main.notebook.page5.microscopeObjectsControls[
      9].GetValue()
  inpfile.write(",")
  inpfile.write(str(INormalDirectionYVal))

  # INormalDirectionZ
  INormalDirectionZVal = self.main.notebook.page5.microscopeObjectsControls[
      10].GetValue()
  inpfile.write(",")
  inpfile.write(str(INormalDirectionZVal))
  inpfile.write("]\n")

  # RAcceleratingVoltage (kV)
  RAcceleratingVoltageVal = self.main.notebook.page5.microscopeObjectsControls[
      11].GetValue()
  inpfile.write("RAcceleratingVoltage (kV) = ")
  inpfile.write(str(RAcceleratingVoltageVal))
  inpfile.write("\n")

  # RAcceptanceAngle (deg)
  RAcceptanceAngleVal = self.main.notebook.page5.microscopeObjectsControls[
      12].GetValue()
  inpfile.write("RAcceptanceAngle (deg)    = ")
  inpfile.write(str(RAcceptanceAngleVal))
  inpfile.write("\n")

  # IMAGE OUTPUT OPTIONS====================================================
  inpfile.write("\n# Image Output Options\n \n")

  # RInitialThickness
  RInitialThicknessVal = self.main.notebook.page6.imageObjectsControls[
      0].GetValue()
  inpfile.write("RInitialThickness         = ")
  inpfile.write(str(RInitialThicknessVal))
  inpfile.write("\n")

  # RFinalThickness
  RFinalThicknessVal = self.main.notebook.page6.imageObjectsControls[
      1].GetValue()
  inpfile.write("RFinalThickness           = ")
  inpfile.write(str(RFinalThicknessVal))
  inpfile.write("\n")

  # RDeltaThickness
  RDeltaThicknessVal = self.main.notebook.page6.imageObjectsControls[
      2].GetValue()
  inpfile.write("RDeltaThickness           = ")
  inpfile.write(str(RDeltaThicknessVal))
  inpfile.write("\n")

  # IReflectOut
  IReflectOutVal = self.main.notebook.page6.imageObjectsControls[
      3].GetValue()
  inpfile.write("IReflectOut               = ")
  inpfile.write(str(IReflectOutVal))
  inpfile.write("\n")

  inpfile.close()

def LoadInputFile(event):
  InputFileDialog = wx.FileDialog(self, "Load input file", "", "",
                                    "INP files (*.inp)|*.inp", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

  # If the User selects cancel
  if InputFileDialog.ShowModal() == wx.ID_CANCEL:
    return

  # get the filename and path from user
  self.User_INPfilename = InputFileDialog.GetFilename()
  self.User_INPpath = InputFileDialog.GetPath()

  self.SetGUIFromFile(self.User_INPpath)

  # shutil.copy2(self.User_INPpath,dir+"/felix.cif")

  # self.main.option.CIFtextbox.SetValue(self.User_CIFpath)

def SetGUIFromFile(dir):

  FelixInpLoad = dir
  inpLoadFile = open(FelixInpLoad, "rb")

  FileLineList = []

  for line in inpLoadFile:
    FileLineList.append(line[28:])
    print(line[28:])

  # IWriteFLAG
  IWriteFLAGSet = int(FileLineList[7])
  self.main.notebook.page1.flagObjectsChoices[0].SetSelection(IWriteFLAGSet)

  # IImageFLAG
  IImageFlagSet = FileLineList[8]
  print('\n')
  print(IImageFlagSet)
  IImageFlagSize = len(IImageFlagSet) - 1
  print('\n')
  print(IImageFlagSize)

  IImageFlagValues = [False, False, False]

  for x in range(0, IImageFlagSize):
    index = int(IImageFlagSet[x])
    IImageFlagValues[index] = True

  print(IImageFlagValues)

  for y in range(0, 3):
    self.main.notebook.page6.imageObjectsControls[
        (y + 4)].SetValue(IImageFlagValues[y])

  # IScatterFactorMethodFLAG
  IScatterFactorMethodFLAGSet = int(FileLineList[9])
  self.main.notebook.page1.flagObjectsChoices[
      1].SetSelection(IScatterFactorMethodFLAGSet)

  # IMaskFLAG
  IMaskFlagSet = FileLineList[10]
  if IMaskFlagSet[0] == '0':
    self.main.notebook.page1.flagObjectsChoices[2].SetValue(False)
  elif IMaskFlagSet[0] == '1':
    self.main.notebook.page1.flagObjectsChoices[2].SetValue(True)

  # IZolzFLAG
  IZolzFLAGSet = FileLineList[11]
  if IZolzFLAGSet[0] == '0':
    self.main.notebook.page1.flagObjectsChoices[3].SetValue(False)
  elif IZolzFLAGSet[0] == '1':
    self.main.notebook.page1.flagObjectsChoices[3].SetValue(True)

  # IAbsorbFLAG
  IAbsorbFLAGSet = int(FileLineList[12])
  self.main.notebook.page1.flagObjectsChoices[4].SetSelection(IAbsorbFLAGSet)

  # IAnisoDebyeWallerFLAG
  IAnisoDebyeWallerFLAGSet = int(FileLineList[13])
  self.main.notebook.page1.flagObjectsChoices[
      5].SetSelection(IAnisoDebyeWallerFLAGSet)

  # IPseudoCubicFLAG
  IPseudoCubicFLAGSet = int(FileLineList[14])
  self.main.notebook.page1.flagObjectsChoices[
      6].SetSelection(IPseudoCubicFLAGSet)

  # IXDirectionFLAG
  IXDirectionFLAGSet = int(FileLineList[15])
  self.main.notebook.page1.flagObjectsChoices[
      7].SetSelection(IXDirectionFLAGSet)

  # RADIUS OF BEAM=========================================================
  # IPixelCount
  IPixelCountSet = float(FileLineList[18])
  self.main.notebook.page2.IPixelCount.SetValue(IPixelCountSet)

  # BEAM SELECTION=========================================================
  # IMinReflectionPool
  IMinReflectionPoolSet = int(FileLineList[21])
  self.main.notebook.page3.beamObjectsControls[
      0].SetValue(IMinReflectionPoolSet)

  # IMinStrongBeams
  IMinStrongBeamsSet = int(FileLineList[22])
  self.main.notebook.page3.beamObjectsControls[
      1].SetValue(IMinStrongBeamsSet)

  # IMinWeakBeams
  IMinWeakBeamsSet = int(FileLineList[23])
  self.main.notebook.page3.beamObjectsControls[2].SetValue(IMinWeakBeamsSet)

  # RBSBMax
  RBSBMaxSet = float(FileLineList[24])
  self.main.notebook.page3.beamObjectsControls[3].SetValue(RBSBMaxSet)

  # RBSPMax
  RBSPMaxSet = float(FileLineList[25])
  self.main.notebook.page3.beamObjectsControls[4].SetValue(RBSPMaxSet)

  # CRYSTAL SETTINGS======================================================

  # RDebyeWallerConstant
  RDebyeWallerConstantSet = float(FileLineList[28])
  self.main.notebook.page4.crystalObjectsControls[
      0].SetValue(RDebyeWallerConstantSet)

  # RAbsorptionPer
  RAbsorptionPerSet = float(FileLineList[29])
  self.main.notebook.page4.crystalObjectsControls[
      1].SetValue(RAbsorptionPerSet)

  # MICROSCOPE SETTINGS====================================================

  # ROuterConvergenceAngle
  ROuterConvergenceAngleSet = float(FileLineList[32])
  self.main.notebook.page5.microscopeObjectsControls[
      0].SetValue(ROuterConvergenceAngleSet)

  # ROuterConvergenceAngle
  RInnerConvergenceAngleSet = float(FileLineList[33])
  self.main.notebook.page5.microscopeObjectsControls[
      1].SetValue(RInnerConvergenceAngleSet)

  # Parsing the beam direction input
  IIncidentBeamDirection = FileLineList[34]
  IIncidentBeamDirection = IIncidentBeamDirection.lstrip('[')
  IIncidentBeamDirection = IIncidentBeamDirection.rstrip('\n')
  IIncidentBeamDirection = IIncidentBeamDirection.rstrip('')
  IIncidentBeamDirection = IIncidentBeamDirection.rstrip(']')
  IIncidentBeamDirection = IIncidentBeamDirection.split(',')

  IXDirection = FileLineList[35]
  IXDirection = IXDirection.lstrip('[')
  IXDirection = IXDirection.rstrip('\n')
  IXDirection = IXDirection.rstrip('')
  IXDirection = IXDirection.rstrip(']')
  IXDirection = IXDirection.split(',')

  INormalDirection = FileLineList[36]
  INormalDirection = INormalDirection.lstrip('[')
  INormalDirection = INormalDirection.rstrip('\n')
  INormalDirection = INormalDirection.rstrip('')
  INormalDirection = INormalDirection.rstrip(']')
  INormalDirection = INormalDirection.split(',')

  # IIncidentBeamDirectionX
  IIncidentBeamDirectionXSet = int(IIncidentBeamDirection[0])
  self.main.notebook.page5.microscopeObjectsControls[
      2].SetValue(IIncidentBeamDirectionXSet)

  # IIncidentBeamDirectionY
  IIncidentBeamDirectionYSet = int(IIncidentBeamDirection[1])
  self.main.notebook.page5.microscopeObjectsControls[
      3].SetValue(IIncidentBeamDirectionYSet)

  # IIncidentBeamDirectionZ
  IIncidentBeamDirectionZSet = int(IIncidentBeamDirection[2])
  self.main.notebook.page5.microscopeObjectsControls[
      4].SetValue(IIncidentBeamDirectionZSet)

  # IXDirectionX
  IXDirectionXSet = int(IXDirection[0])
  self.main.notebook.page5.microscopeObjectsControls[
      5].SetValue(IXDirectionXSet)

  # IXDirectionY
  IXDirectionYSet = int(IXDirection[1])
  self.main.notebook.page5.microscopeObjectsControls[
      6].SetValue(IXDirectionYSet)

  # IXDirectionZ
  IXDirectionZSet = int(IXDirection[2])
  self.main.notebook.page5.microscopeObjectsControls[
      7].SetValue(IXDirectionZSet)

  # INormalDirectionX
  INormalDirectionXSet = int(INormalDirection[0])
  self.main.notebook.page5.microscopeObjectsControls[
      8].SetValue(INormalDirectionXSet)

  # INormalDirectionY
  INormalDirectionYSet = int(INormalDirection[1])
  self.main.notebook.page5.microscopeObjectsControls[
      9].SetValue(INormalDirectionYSet)

  # INormalDirectionZ
  INormalDirectionZSet = int(INormalDirection[2])
  self.main.notebook.page5.microscopeObjectsControls[
      10].SetValue(INormalDirectionZSet)

  # RAcceleratingVoltage (kV)
  RAcceleratingVoltageSet = float(FileLineList[37])
  self.main.notebook.page5.microscopeObjectsControls[
      11].SetValue(RAcceleratingVoltageSet)

  # RAcceptanceAngle (deg)
  RAcceptanceAngleSet = float(FileLineList[38])
  self.main.notebook.page5.microscopeObjectsControls[
      12].SetValue(RAcceptanceAngleSet)

  # IMAGE OUTPUT OPTIONS====================================================

  # RInitialThickness
  RInitialThicknessSet = float(FileLineList[42])
  self.main.notebook.page6.imageObjectsControls[
      0].SetValue(RInitialThicknessSet)

  # RFinalThickness
  RFinalThicknessSet = float(FileLineList[43])
  self.main.notebook.page6.imageObjectsControls[
      1].SetValue(RFinalThicknessSet)

  # RDeltaThickness
  RDeltaThicknessSet = float(FileLineList[44])
  self.main.notebook.page6.imageObjectsControls[
      2].SetValue(RDeltaThicknessSet)

  # IReflectOut
  IReflectOutSet = int(FileLineList[45])
  self.main.notebook.page6.imageObjectsControls[3].SetValue(IReflectOutSet)

  inpLoadFile.close()

def OnClose(event):

  self.Close(True)

def OnCif(event):

  CIFFileDialog = wx.FileDialog(self, "Load CIF file", "", "",
                                  "CIF files (*.cif)|*.cif", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

  # If the User selects cancel
  if CIFFileDialog.ShowModal() == wx.ID_CANCEL:
    return

  # get the filename and path from user
  self.User_CIFfilename = CIFFileDialog.GetFilename()
  self.User_CIFpath = CIFFileDialog.GetPath()

  self.main.option.CIFtextbox.SetValue(self.User_CIFpath)
