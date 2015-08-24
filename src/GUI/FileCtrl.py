#!/usr/bin/python
# -*- coding: utf-8 -*-

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# felixsim
#
# Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
#
# (C) 2013/14, all right reserved
#
# Version: :VERSION:
# Date:    :DATE:
# Time:    :TIME:
# Status:  :RLSTATUS:
# Build:   :BUILD:
# Author:  :AUTHOR:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  This file is part of felixsim.
#
#  felixsim is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  felixsim is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import wx
import wx.lib.agw.floatspin as FS
import math
import os
import shutil
import sys
import GuiPages
# import Bin2Tiff
import Felix_gui
import datetime


# Class defining gui which selects where the user wants to save the input file
# When Felix is run, this file will be copied as felix.inp

def RunFelix(parent, CIFPath, OutputPath):

  # opens new (re-existing) directory and puts all input files into it
  #***ensure file and path check needs to be put here****
  filePath = os.path.realpath(__file__)
  fileDirectory, fileName = os.path.split(filePath)

  if parent.main.notebook.page6.imageObjectsControls[4].GetValue() == False and parent.main.notebook.page6.imageObjectsControls[5].GetValue() == False and parent.main.notebook.page6.imageObjectsControls[6].GetValue() == False:
    wx.MessageBox('Please select montage, stack reflections, or amplitude and phase from the image panel (at least one must be selected).', 'Error',
                  wx.OK | wx.ICON_ERROR)
    return

  parent.main.wiki.felixStatus = 1

  if OutputPath == None:

    wx.MessageBox('Please select an output directory', 'Error',
                  wx.OK | wx.ICON_ERROR)

    OutputPath = OutputDirSelect(parent)

  if CIFPath == None:
    wx.MessageBox('Please load a .cif file', 'Error',
                  wx.OK | wx.ICON_ERROR)
    CIFPath = OnCif(parent)

  if CIFPath !=None and os.path.exists(CIFPath) == False:
    wx.MessageBox('.cif file not found, please load a .cif file', 'Error',
                  wx.OK | wx.ICON_ERROR)
    CIFPath = OnCif(parent)

  # creates working directory
  cpath, cfilename = os.path.split(CIFPath)
  cname = cfilename[:-4]
  OutputDirectory = OutputPath + "/" + cname + "_output"
  workingDir = OutputDirectory + "/temp/Working_Directory/"

  count = 0
  OriginalDir = OutputDirectory
  while os.path.exists(OutputDirectory):
    count = count + 1
    OutputDirectory = OriginalDir + "(" + str(count) + ")"

  if os.path.exists(workingDir):
    shutil.rmtree(workingDir)
  os.makedirs(workingDir)

  # Change to working directory
  os.chdir(workingDir)

  # copy cif file from user specified location to new directory
  shutil.copy2(CIFPath, workingDir + "/felix.cif")

  # Write input file
  InpName = "felix.inp"
  WriteInputFile(workingDir, parent, InpName)
  wx.MessageBox('Files created successfully, click okay to run Felix', 'Info',
                wx.OK | wx.CANCEL | wx.ICON_INFORMATION)

  # Get Value of mpicores
  NumberofCores = parent.MPICores.GetValue()

  # Run in parallel or single core
  now = datetime.datetime.now().strftime("%I_%M%p_%B_%d_%Y")
  logname = "log_" + now
  if os.path.exists(OutputDirectory + "/logs") == False:
    os.makedirs(OutputDirectory + "/logs")

  print "Files before run: " + str(os.listdir(workingDir))

  if NumberofCores == 1:
    os.system("../../../../../felixsim > " + OutputDirectory + "/logs/" + logname)  # single core
  else:
    os.system("mpirun -n " + str(NumberofCores) +
              " ../../../../../felixsim  > " + OutputDirectory + "/logs" + logname)  # parallel

  print "Files after run: " + str(os.listdir(workingDir))

  shutil.copytree(workingDir, OutputDirectory)
  os.rename(OutputDirectory + "/felix.cif", OutputDirectory + "/" + cfilename)

  if os.path.exists(workingDir):
    shutil.rmtree(workingDir)

  os.chdir(fileDirectory)

  wx.MessageBox('Felixsim successfully run!',
                'Info', wx.OK | wx.ICON_INFORMATION)

  parent.main.wiki.felixStatus = 0

  return

  # os.chdir("../")
  # Bin2Tiff.convert(dir, '8', 'tif', '1')

  # parent.main.viewer.onView(dir)
  # parent.Close(True)


def InpCreate(parent):

  if parent.main.notebook.page6.imageObjectsControls[4].GetValue() == False and parent.main.notebook.page6.imageObjectsControls[5].GetValue() == False and parent.main.notebook.page6.imageObjectsControls[6].GetValue() == False:
    wx.MessageBox('Please select montage, stack reflections, or amplitude and phase from the image panel (at least one must be selected).', 'Error',
                  wx.OK | wx.ICON_ERROR)
    return

  InputSaveDialog = wx.FileDialog(parent, "Save Input file", "", "felix.inp",
                                  "Inp Files (*.inp)|*.inp", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

  InputSaveDialog.ShowModal()

  if InputSaveDialog == wx.ID_CANCEL:
    return

  User_InpSavePath, User_InpSaveName = os.path.split(InputSaveDialog.GetPath())

  SaveStatus = WriteInputFile(User_InpSavePath, parent, User_InpSaveName)

  if SaveStatus == True:
    wx.MessageBox('Input File successfully written',
                  'Info', wx.OK | wx.ICON_INFORMATION)


def WriteInputFile(dir, parent, name):

  FelixInpFilename = dir + "/" + name

  if os.path.exists(FelixInpFilename):
    os.remove(FelixInpFilename)

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
  IWriteFLAGVal = parent.main.notebook.page1.flagObjectsChoices[
      0].GetCurrentSelection()
  inpfile.write("IWriteFLAG                = ")
  inpfile.write(str(IWriteFLAGVal))
  inpfile.write("\n")

  # IImageFLAG
  inpfile.write("IImageFLAG                = ")

  IImageFLAGVal0 = parent.main.notebook.page6.imageObjectsControls[
      4].GetValue()
  if IImageFLAGVal0 == True:
    inpfile.write("0")

  IImageFLAGVal1 = parent.main.notebook.page6.imageObjectsControls[
      5].GetValue()
  if IImageFLAGVal1 == True:
    inpfile.write("1")

  IImageFLAGVal2 = parent.main.notebook.page6.imageObjectsControls[
      6].GetValue()
  if IImageFLAGVal2 == True:
    inpfile.write("2")

  inpfile.write("\n")

  # IScatterFactorMethodFLAG
  IScatterFactorMethodFLAGVal = parent.main.notebook.page1.flagObjectsChoices[
      1].GetCurrentSelection()
  inpfile.write("IScatterFactorMethodFLAG  = ")
  inpfile.write(str(IScatterFactorMethodFLAGVal))
  inpfile.write("\n")

  # IMaskFLAG
  IMaskFLAGVal = parent.main.notebook.page1.flagObjectsChoices[2].GetValue()
  inpfile.write("IMaskFLAG                 = ")
  if IMaskFLAGVal == True:
    IMaskFLAGBool = 1
  elif IMaskFLAGVal == False:
    IMaskFLAGBool = 0
  inpfile.write(str(IMaskFLAGBool))
  inpfile.write("\n")

  # IZolzFLAG
  IZolzFLAGVal = parent.main.notebook.page1.flagObjectsChoices[3].GetValue()
  inpfile.write("IZolzFLAG                 = ")
  if IZolzFLAGVal == True:
    IZolzFLAGBool = 1
  elif IZolzFLAGVal == False:
    IZolzFLAGBool = 0
  inpfile.write(str(IZolzFLAGBool))
  inpfile.write("\n")

  # IAbsorbFLAG
  IAbsorbFLAGVal = parent.main.notebook.page1.flagObjectsChoices[
      4].GetCurrentSelection()
  inpfile.write("IAbsorbFLAG               = ")
  inpfile.write(str(IAbsorbFLAGVal))
  inpfile.write("\n")

  # IAnisoDebyeWallerFLAG
  IAnisoDebyeWallerFLAGVal = parent.main.notebook.page1.flagObjectsChoices[
      5].GetCurrentSelection()
  inpfile.write("IAnisoDebyeWallerFLAG     = ")
  inpfile.write(str(IAnisoDebyeWallerFLAGVal))
  inpfile.write("\n")

  # IPseudoCubicFLAG
  IPseudoCubicFLAGVal = parent.main.notebook.page1.flagObjectsChoices[
      6].GetCurrentSelection()
  inpfile.write("IPseudoCubicFLAG          = ")
  inpfile.write(str(IPseudoCubicFLAGVal))
  inpfile.write("\n")

  # IXDirectionFLAG
  IXDirectionFLAGVal = parent.main.notebook.page1.flagObjectsChoices[
      7].GetCurrentSelection()
  inpfile.write("IXDirectionFLAG           = ")
  inpfile.write(str(IXDirectionFLAGVal))
  inpfile.write("\n")

  # RADIUS OF BEAM=========================================================
  inpfile.write("\n# radius of the beam in pixels\n")

  # IPixelCount
  IPixelCountVal = int(parent.main.notebook.page2.IPixelCount.GetValue())
  inpfile.write("IPixelCount               = ")
  inpfile.write(str(IPixelCountVal))
  inpfile.write("\n")

  # BEAM SELECTION=========================================================
  inpfile.write("\n# beam selection criteria\n")

  # IMinReflectionPool
  IMinReflectionPoolVal = parent.main.notebook.page3.beamObjectsControls[
      0].GetValue()
  inpfile.write("IMinReflectionPool        = ")
  inpfile.write(str(IMinReflectionPoolVal))
  inpfile.write("\n")

  # IMinStrongBeams
  IMinStrongBeamsVal = parent.main.notebook.page3.beamObjectsControls[
      1].GetValue()
  inpfile.write("IMinStrongBeams           = ")
  inpfile.write(str(IMinStrongBeamsVal))
  inpfile.write("\n")

  # IMinWeakBeams
  IMinWeakBeamsVal = parent.main.notebook.page3.beamObjectsControls[
      2].GetValue()
  inpfile.write("IMinWeakBeams             = ")
  inpfile.write(str(IMinWeakBeamsVal))
  inpfile.write("\n")

  # RBSBMax
  RBSBMaxVal = parent.main.notebook.page3.beamObjectsControls[3].GetValue()
  inpfile.write("RBSBMax                   = ")
  inpfile.write(str(RBSBMaxVal))
  inpfile.write("\n")

  # RBSPMax
  RBSPMaxVal = parent.main.notebook.page3.beamObjectsControls[4].GetValue()
  inpfile.write("RBSPMax                   = ")
  inpfile.write(str(RBSPMaxVal))
  inpfile.write("\n")

  # CRYSTAL SETTINGS======================================================
  inpfile.write("\n# crystal settings\n")

  # RDebyeWallerConstant
  RDebyeWallerConstantVal = parent.main.notebook.page4.crystalObjectsControls[
      0].GetValue()
  inpfile.write("RDebyeWallerConstant      = ")
  inpfile.write(str(RDebyeWallerConstantVal))
  inpfile.write("\n")

  # RAbsorptionPer
  RAbsorptionPerVal = parent.main.notebook.page4.crystalObjectsControls[
      1].GetValue()
  inpfile.write("RAbsorptionPer            = ")
  inpfile.write(str(RAbsorptionPerVal))
  inpfile.write("\n")

  # MICROSCOPE SETTINGS====================================================
  inpfile.write("\n# microscope settings\n")

  # ROuterConvergenceAngle
  ROuterConvergenceAngleVal = parent.main.notebook.page5.microscopeObjectsControls[
      0].GetValue()
  inpfile.write("ROuterConvergenceAngle    = ")
  inpfile.write(str(ROuterConvergenceAngleVal))
  inpfile.write("\n")

  # ROuterConvergenceAngle
  RInnerConvergenceAngleVal = parent.main.notebook.page5.microscopeObjectsControls[
      1].GetValue()
  inpfile.write("RInnerConvergenceAngle    = ")
  inpfile.write(str(RInnerConvergenceAngleVal))
  inpfile.write("\n")

  # IIncidentBeamDirectionX
  IIncidentBeamDirectionXVal = parent.main.notebook.page5.microscopeObjectsControls[
      2].GetValue()
  inpfile.write("IIncidentBeamDirection    = [")
  inpfile.write(str(IIncidentBeamDirectionXVal))

  # IIncidentBeamDirectionY
  IIncidentBeamDirectionYVal = parent.main.notebook.page5.microscopeObjectsControls[
      3].GetValue()
  inpfile.write(",")
  inpfile.write(str(IIncidentBeamDirectionYVal))

  # IIncidentBeamDirectionZ
  IIncidentBeamDirectionZVal = parent.main.notebook.page5.microscopeObjectsControls[
      4].GetValue()
  inpfile.write(",")
  inpfile.write(str(IIncidentBeamDirectionZVal))
  inpfile.write("]\n")

  # IXDirectionX
  IXDirectionXVal = parent.main.notebook.page5.microscopeObjectsControls[
      5].GetValue()
  inpfile.write("IXDirection               = [")
  inpfile.write(str(IXDirectionXVal))

  # IXDirectionY
  IXDirectionYVal = parent.main.notebook.page5.microscopeObjectsControls[
      6].GetValue()
  inpfile.write(",")
  inpfile.write(str(IXDirectionYVal))

  # IXDirectionZ
  IXDirectionZVal = parent.main.notebook.page5.microscopeObjectsControls[
      7].GetValue()
  inpfile.write(",")
  inpfile.write(str(IXDirectionZVal))
  inpfile.write("]\n")

  # INormalDirectionX
  INormalDirectionXVal = parent.main.notebook.page5.microscopeObjectsControls[
      8].GetValue()
  inpfile.write("INormalDirectionX         = [")
  inpfile.write(str(INormalDirectionXVal))

  # INormalDirectionY
  INormalDirectionYVal = parent.main.notebook.page5.microscopeObjectsControls[
      9].GetValue()
  inpfile.write(",")
  inpfile.write(str(INormalDirectionYVal))

  # INormalDirectionZ
  INormalDirectionZVal = parent.main.notebook.page5.microscopeObjectsControls[
      10].GetValue()
  inpfile.write(",")
  inpfile.write(str(INormalDirectionZVal))
  inpfile.write("]\n")

  # RAcceleratingVoltage (kV)
  RAcceleratingVoltageVal = parent.main.notebook.page5.microscopeObjectsControls[
      11].GetValue()
  inpfile.write("RAcceleratingVoltage (kV) = ")
  inpfile.write(str(RAcceleratingVoltageVal))
  inpfile.write("\n")

  # RAcceptanceAngle (deg)
  RAcceptanceAngleVal = parent.main.notebook.page5.microscopeObjectsControls[
      12].GetValue()
  inpfile.write("RAcceptanceAngle (deg)    = ")
  inpfile.write(str(RAcceptanceAngleVal))
  inpfile.write("\n")

  # IMAGE OUTPUT OPTIONS====================================================
  inpfile.write("\n# Image Output Options\n \n")

  # RInitialThickness
  RInitialThicknessVal = parent.main.notebook.page6.imageObjectsControls[
      0].GetValue()
  inpfile.write("RInitialThickness         = ")
  inpfile.write(str(RInitialThicknessVal))
  inpfile.write("\n")

  # RFinalThickness
  RFinalThicknessVal = parent.main.notebook.page6.imageObjectsControls[
      1].GetValue()
  inpfile.write("RFinalThickness           = ")
  inpfile.write(str(RFinalThicknessVal))
  inpfile.write("\n")

  # RDeltaThickness
  RDeltaThicknessVal = parent.main.notebook.page6.imageObjectsControls[
      2].GetValue()
  inpfile.write("RDeltaThickness           = ")
  inpfile.write(str(RDeltaThicknessVal))
  inpfile.write("\n")

  # IReflectOut
  IReflectOutVal = parent.main.notebook.page6.imageObjectsControls[
      3].GetValue()
  inpfile.write("IReflectOut               = ")
  inpfile.write(str(IReflectOutVal))
  inpfile.write("\n")

  inpfile.close()

  return True


def LoadInputFile(parent):
  InputFileDialog = wx.FileDialog(parent, "Load input file", "", "",
                                  "INP files (*.inp)|*.inp", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

  # If the User selects cancel
  if InputFileDialog.ShowModal() == wx.ID_CANCEL:
    return

  # get the filename and path from user
  User_INPfilename = InputFileDialog.GetFilename()
  User_INPpath = InputFileDialog.GetPath()

  SetGUIFromFile(User_INPpath, parent)

  # shutil.copy2(self.User_INPpath,dir+"/felix.cif")

  # self.main.option.CIFtextbox.SetValue(self.User_CIFpath)


def SetGUIFromFile(dir, parent):

  FelixInpLoad = dir
  inpLoadFile = open(FelixInpLoad, "rb")

  FileLineList = []

  for line in inpLoadFile:
    if line.find('=') != -1:
      FileLineList.append(line[line.find('=') + 2:])
    else:
      FileLineList.append(line)

  # IWriteFLAG
  IWriteFLAGSet = int(FileLineList[7])
  parent.main.notebook.page1.flagObjectsChoices[0].SetSelection(IWriteFLAGSet)

  # IImageFLAG
  IImageFlagSet = FileLineList[8]
  IImageFlagSize = len(IImageFlagSet) - 1

  IImageFlagValues = [False, False, False]

  for x in range(0, IImageFlagSize):
    index = int(IImageFlagSet[x])
    IImageFlagValues[index] = True

  for y in range(0, 3):
    parent.main.notebook.page6.imageObjectsControls[
        (y + 4)].SetValue(IImageFlagValues[y])

  # IScatterFactorMethodFLAG
  IScatterFactorMethodFLAGSet = int(FileLineList[9])
  parent.main.notebook.page1.flagObjectsChoices[
      1].SetSelection(IScatterFactorMethodFLAGSet)

  # IMaskFLAG
  IMaskFlagSet = FileLineList[10]
  if IMaskFlagSet[0] == '0':
    parent.main.notebook.page1.flagObjectsChoices[2].SetValue(False)
  elif IMaskFlagSet[0] == '1':
    parent.main.notebook.page1.flagObjectsChoices[2].SetValue(True)

  # IZolzFLAG
  IZolzFLAGSet = FileLineList[11]
  if IZolzFLAGSet[0] == '0':
    parent.main.notebook.page1.flagObjectsChoices[3].SetValue(False)
  elif IZolzFLAGSet[0] == '1':
    parent.main.notebook.page1.flagObjectsChoices[3].SetValue(True)

  # IAbsorbFLAG
  IAbsorbFLAGSet = int(FileLineList[12])
  parent.main.notebook.page1.flagObjectsChoices[4].SetSelection(IAbsorbFLAGSet)

  # IAnisoDebyeWallerFLAG
  IAnisoDebyeWallerFLAGSet = int(FileLineList[13])
  parent.main.notebook.page1.flagObjectsChoices[
      5].SetSelection(IAnisoDebyeWallerFLAGSet)

  # IPseudoCubicFLAG
  IPseudoCubicFLAGSet = int(FileLineList[14])
  parent.main.notebook.page1.flagObjectsChoices[
      6].SetSelection(IPseudoCubicFLAGSet)

  # IXDirectionFLAG
  IXDirectionFLAGSet = int(FileLineList[15])
  parent.main.notebook.page1.flagObjectsChoices[
      7].SetSelection(IXDirectionFLAGSet)

  # RADIUS OF BEAM=========================================================
  # IPixelCount
  IPixelCountSet = float(FileLineList[18])
  parent.main.notebook.page2.IPixelCount.SetValue(IPixelCountSet)

  # BEAM SELECTION=========================================================
  # IMinReflectionPool
  IMinReflectionPoolSet = int(FileLineList[21])
  parent.main.notebook.page3.beamObjectsControls[
      0].SetValue(IMinReflectionPoolSet)

  # IMinStrongBeams
  IMinStrongBeamsSet = int(FileLineList[22])
  parent.main.notebook.page3.beamObjectsControls[
      1].SetValue(IMinStrongBeamsSet)

  # IMinWeakBeams
  IMinWeakBeamsSet = int(FileLineList[23])
  parent.main.notebook.page3.beamObjectsControls[2].SetValue(IMinWeakBeamsSet)

  # RBSBMax
  RBSBMaxSet = float(FileLineList[24])
  parent.main.notebook.page3.beamObjectsControls[3].SetValue(RBSBMaxSet)

  # RBSPMax
  RBSPMaxSet = float(FileLineList[25])
  parent.main.notebook.page3.beamObjectsControls[4].SetValue(RBSPMaxSet)

  # CRYSTAL SETTINGS======================================================

  # RDebyeWallerConstant
  RDebyeWallerConstantSet = float(FileLineList[28])
  parent.main.notebook.page4.crystalObjectsControls[
      0].SetValue(RDebyeWallerConstantSet)

  # RAbsorptionPer
  RAbsorptionPerSet = float(FileLineList[29])
  parent.main.notebook.page4.crystalObjectsControls[
      1].SetValue(RAbsorptionPerSet)

  # MICROSCOPE SETTINGS====================================================

  # ROuterConvergenceAngle
  ROuterConvergenceAngleSet = float(FileLineList[32])
  parent.main.notebook.page5.microscopeObjectsControls[
      0].SetValue(ROuterConvergenceAngleSet)

  # ROuterConvergenceAngle
  RInnerConvergenceAngleSet = float(FileLineList[33])
  parent.main.notebook.page5.microscopeObjectsControls[
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
  parent.main.notebook.page5.microscopeObjectsControls[
      2].SetValue(IIncidentBeamDirectionXSet)

  # IIncidentBeamDirectionY
  IIncidentBeamDirectionYSet = int(IIncidentBeamDirection[1])
  parent.main.notebook.page5.microscopeObjectsControls[
      3].SetValue(IIncidentBeamDirectionYSet)

  # IIncidentBeamDirectionZ
  IIncidentBeamDirectionZSet = int(IIncidentBeamDirection[2])
  parent.main.notebook.page5.microscopeObjectsControls[
      4].SetValue(IIncidentBeamDirectionZSet)

  # IXDirectionX
  IXDirectionXSet = int(IXDirection[0])
  parent.main.notebook.page5.microscopeObjectsControls[
      5].SetValue(IXDirectionXSet)

  # IXDirectionY
  IXDirectionYSet = int(IXDirection[1])
  parent.main.notebook.page5.microscopeObjectsControls[
      6].SetValue(IXDirectionYSet)

  # IXDirectionZ
  IXDirectionZSet = int(IXDirection[2])
  parent.main.notebook.page5.microscopeObjectsControls[
      7].SetValue(IXDirectionZSet)

  # INormalDirectionX
  INormalDirectionXSet = int(INormalDirection[0])
  parent.main.notebook.page5.microscopeObjectsControls[
      8].SetValue(INormalDirectionXSet)

  # INormalDirectionY
  INormalDirectionYSet = int(INormalDirection[1])
  parent.main.notebook.page5.microscopeObjectsControls[
      9].SetValue(INormalDirectionYSet)

  # INormalDirectionZ
  INormalDirectionZSet = int(INormalDirection[2])
  parent.main.notebook.page5.microscopeObjectsControls[
      10].SetValue(INormalDirectionZSet)

  # RAcceleratingVoltage (kV)
  RAcceleratingVoltageSet = float(FileLineList[37])
  parent.main.notebook.page5.microscopeObjectsControls[
      11].SetValue(RAcceleratingVoltageSet)

  # RAcceptanceAngle (deg)
  RAcceptanceAngleSet = float(FileLineList[38])
  parent.main.notebook.page5.microscopeObjectsControls[
      12].SetValue(RAcceptanceAngleSet)

  # IMAGE OUTPUT OPTIONS====================================================

  # RInitialThickness
  RInitialThicknessSet = float(FileLineList[42])
  parent.main.notebook.page6.imageObjectsControls[
      0].SetValue(RInitialThicknessSet)

  # RFinalThickness
  RFinalThicknessSet = float(FileLineList[43])
  parent.main.notebook.page6.imageObjectsControls[
      1].SetValue(RFinalThicknessSet)

  # RDeltaThickness
  RDeltaThicknessSet = float(FileLineList[44])
  parent.main.notebook.page6.imageObjectsControls[
      2].SetValue(RDeltaThicknessSet)

  # IReflectOut
  IReflectOutSet = int(FileLineList[45])
  parent.main.notebook.page6.imageObjectsControls[3].SetValue(IReflectOutSet)

  inpLoadFile.close()


def OnCif(parent):

  filePath = os.path.realpath(__file__)
  fileDirectory, fileName = os.path.split(filePath)

  CIFFileDialog = wx.FileDialog(parent, "Load CIF file", "", "",
                                "CIF files (*.cif)|*.cif", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

  CIFFileDialog.ShowModal()

  # If the User selects cancel
  if CIFFileDialog == wx.ID_CANCEL:
    return

  # get the filename and path from user
  User_CIFfilename = CIFFileDialog.GetFilename()
  User_CIFpath = CIFFileDialog.GetPath()

  # parent.main.option.CIFtextbox.SetValue(User_CIFpath)
  return User_CIFpath


def OutputDirSelect(parent):

  filePath = os.path.realpath(__file__)
  fileDirectory, fileName = os.path.split(filePath)
  OutputDirDialog = wx.DirDialog(
      parent, "Select output directory")
  OutputDirDialog.ShowModal()

  # If the User selects cancel
  if OutputDirDialog == wx.ID_CANCEL:
    return

  # get the filename and path from user
  User_OutputPath = OutputDirDialog.GetPath()

  return User_OutputPath
