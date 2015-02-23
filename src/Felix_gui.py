# -*- coding: utf-8 -*-
"""

@author: Alex
"""

#!/usr/bin/python
# -*- coding: utf-8 -*-

#module loads - wx is gui interpreter
#             - FS is floatspin attribute in wx
#             - os provides the means to create files
#             - shutil adds functions to os 
import wx
import wx.lib.agw.floatspin as FS
import os
import shutil


#Class defining gui which selects where the user wants to save the input file
#When Felix is run, this file will be copied as felix.inp
class WriteInputDialog(wx.Dialog):

    def __init__(self, parent,id, title):
        wx.Dialog.__init__(self,parent,id,title,size=(400,200))

    
        self.UserInterface2()
        #self.SetSize(400,200)
        #self.SetTitle("Save Input File As")

    def UserInterface2(self):


        wx.StaticText(self, label="Choose location for Input file", pos=(5,20),style=wx.ALIGN_CENTRE_HORIZONTAL)
        self.InputFileTextBox=wx.TextCtrl(self, -1,os.getcwd(), pos=(5,100), size=(250,-1),style=wx.TE_RIGHT)
        InputFileBrowse=wx.Button(self, label='Save', pos=(270,100),size=(60,-1))
        InputFileBrowse.Bind(wx.EVT_BUTTON, self.OnINPBrowse)

        ReturnButton=wx.Button(self, label='OK', pos=(200,130),size=(60,-1))
        CancelButton=wx.Button(self, label='Cancel', pos=(300,130),size=(60,-1))

        ReturnButton.Bind(wx.EVT_BUTTON,self.ReturnSuccess)
        CancelButton.Bind(wx.EVT_BUTTON,self.ReturnCancel)
        

    def OnINPBrowse(self,event):
        InputFileDialog=wx.FileDialog(self, "Save Input file", "", "",
                                       "Inp Files (*.inp)|*.inp", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if InputFileDialog.ShowModal() == wx.ID_CANCEL:
            #IErr=-1
            return()#IErr)

        UserInputPath=InputFileDialog.GetPath()

        self.InputFileTextBox.SetValue(UserInputPath)



    def ReturnSuccess(self,event):

        global SInpPath
        global SInpFilename
        SInpPath, SInpFilename = os.path.split(self.InputFileTextBox.GetValue())
       
        if SInpFilename.find(".inp")==-1:
            wx.MessageBox('.inp file not found, please save a .inp file','Error', wx.OK | wx.ICON_ERROR)
            return()
       
        self.Close()
     

    def ReturnCancel(self,event):

        self.Close()
            
        

    

class Felix(wx.Frame):
           
    def __init__(self, *args, **kw):
        super(Felix, self).__init__(*args, **kw) 
        
        self.InitUI()
    #User Interface set up    
    def InitUI(self):   

        pnl = wx.Panel(self)
        
        #Borders for gui set up
        wx.StaticBox(pnl, label = 'FLAGs', pos=(5, 5), size=(650, 150))
        wx.StaticBox(pnl, label = 'Radius of the Beam in Pixels', pos=(5,160), size=(200,70))
        wx.StaticBox(pnl, label = 'Beam Selection', pos=(5,250), size=(440,130))
        wx.StaticBox(pnl, label = 'Crystal Settings', pos=(450,250), size=(205,130))
        wx.StaticBox(pnl, label = 'Microscope Settings', pos=(5,390), size=(650,130))
        wx.StaticBox(pnl, label = 'Image Output Options', pos=(5,530), size=(400,100))
        wx.StaticBox(pnl, pos = (450,510), size =(205,120))
        
        #list for flag options
        flags= ['0','1','2','3']
        writeflags = ['0','1','2','3','4','10','100','101','102','103','104','110']
        
        #Flag Gui section
        wx.StaticText(pnl, label='IWriteFLAG ', pos=(15, 30))
        self.IWriteFLAG = wx.ComboBox(pnl, value='0', pos=(140, 25),choices=writeflags, size=(60, -1),
                    style=wx.CB_DROPDOWN)        
        
        
        self.IImageFLAG = wx.ComboBox(pnl, pos=(140, 55),value="0", choices=flags[0:3],
            style=wx.CB_DROPDOWN,size=(60, -1))
        #self.st = wx.StaticText(pnl,label='',pos=(15,40))
        #IImageFLAG.Bind(wx.EVT_COMBOBOX, self.OnSelect)
        
        wx.StaticText(pnl, label='IImageFLAG ', pos=(15, 60))
        
        wx.StaticText(pnl, label='IOutputFLAG', pos=(15, 90))
        self.IOutputFLAG = wx.ComboBox(pnl, value='0', pos=(140, 85),choices=flags[0], size=(60, -1),
                    style=wx.CB_DROPDOWN)
       
        wx.StaticText(pnl, label='IBinorTextFLAG', pos=(15, 120))
        self.IBinorTextFLAG = wx.ComboBox(pnl, value='0', pos=(140, 115),choices=flags[0], size=(60, -1),
                    style=wx.CB_DROPDOWN)
                    
        wx.StaticText(pnl, label='IScatterFactorMethodFLAG', pos=(220, 30))
        self.IScatterFactorMethodFLAG = wx.ComboBox(pnl, value='0', pos=(380, 25),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)
                             
        wx.StaticText(pnl, label='IZolzFLAG', pos=(220, 60))
        self.IZolzFLAG = wx.ComboBox(pnl, value='1', pos=(380, 55),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)
                    
        wx.StaticText(pnl, label='ICentralBeamFLAG', pos=(220, 90))
        self.ICentralBeamFLAG = wx.ComboBox(pnl, value='1', pos=(380, 85),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)
        
        wx.StaticText(pnl, label='IMaskFLAG', pos=(220, 120))
        self.IMaskFLAG = wx.ComboBox(pnl, value='0', pos=(380, 115),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)
                    
        wx.StaticText(pnl, label='IAbsorbFLAG', pos=(460, 30))
        self.IAbsorbFLAG = wx.ComboBox(pnl, value='1', pos=(600, 25),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)           
                    
        wx.StaticText(pnl, label='IAnisoDebyeWallerFLAG', pos=(460, 60))
        self.IAnisoDebyeWallerFLAG = wx.ComboBox(pnl, value='0', pos=(600, 55),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)
                    
        wx.StaticText(pnl, label='IPseudoCubicFLAG', pos=(460, 90))
        self.IPseudoCubicFLAG = wx.ComboBox(pnl, value='0', pos=(600, 85),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN)
                    
        wx.StaticText(pnl, label='IXDirectionFLAG', pos=(460, 120))
        self.IXDirectionFLAG = wx.ComboBox(pnl, value='0', pos=(600, 115),choices=flags[0:2], size=(60, -1),
                    style=wx.CB_DROPDOWN) 
                    
        #Radius of beam
        wx.StaticText(pnl, label='IPixelCount', pos=(15, 185))
        self.IPixelCount = FS.FloatSpin(pnl, pos=(100,185), size=(60,-1), value=64, min_val=0, max_val=512,
                                   increment=64,agwStyle=FS.FS_RIGHT)
        self.IPixelCount.SetFormat("%f")
        self.IPixelCount.SetDigits(0)
        
        #Beam Selection
        wx.StaticText(pnl, label='IMinReflectionPool',pos=(15, 270))
        self.IMinReflectionPool = wx.SpinCtrl(pnl, pos=(135,265), size=(60,-1), value='15', min=0, max=100000)

        wx.StaticText(pnl, label='IMinStrongBeams',pos=(15, 300))
        self.IMinStrongBeams = wx.SpinCtrl(pnl, pos=(135,295), size=(60,-1), value='7', min=0, max=100000)
        
        wx.StaticText(pnl, label='IMinWeakBeams',pos=(15, 330))
        self.IMinWeakBeams = wx.SpinCtrl(pnl, pos=(135,325), size=(60,-1), value='5', min=0, max=100000)
        
        wx.StaticText(pnl, label='RBSBMax',pos=(220, 270))
        self.RBSBMax = FS.FloatSpin(pnl, pos=(380,265), size=(60,-1), value=0.1, increment=0.1, min_val=0,
                               max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RBSBMax.SetFormat("%f")
        self.RBSBMax.SetDigits(1)                                                    
        
        wx.StaticText(pnl, label='RBSPMax',pos=(220, 300))
        self.RBSPMax = FS.FloatSpin(pnl, pos=(380,295), size=(60,-1), value=0.1, increment=0.1, min_val=0,
                               max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RBSPMax.SetFormat("%f")
        self.RBSPMax.SetDigits(1)
                
        wx.StaticText(pnl, label='RConvergenceTolerance',pos=(220, 330))
        self.RConvergenceTolerance = FS.FloatSpin(pnl, pos=(380,325), size=(60,-1), value=0.1, increment=0.1, min_val=0,
                               max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RConvergenceTolerance.SetFormat("%f")
        self.RConvergenceTolerance.SetDigits(1)
        
        
        #Crystal Settings
        wx.StaticText(pnl, label='RDebyeWallerConstant', pos=(460, 270))
        self.RDebyeWallerConstant = FS.FloatSpin(pnl, value=0.4668, pos=(590, 265),size=(60,-1),
                                            increment=0.001, min_val=0,max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RDebyeWallerConstant.SetFormat("%f")
        self.RDebyeWallerConstant.SetDigits(4) 

           
        wx.StaticText(pnl, label='RAbsorptionPer', pos=(460, 300))
        self.RAbsorptionPer = FS.FloatSpin(pnl, value=2.9, pos=(590, 295),size=(60,-1),
                                            increment=0.1, min_val=0,max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RAbsorptionPer.SetFormat("%f")
        self.RAbsorptionPer.SetDigits(1) 
        
        
        #Microscope Settings
        
        wx.StaticText(pnl, label='ROuterConvergenceAngle', pos=(15, 415))
        self.ROuterConvergenceAngle = FS.FloatSpin(pnl, value=3.0, pos=(155, 410),size=(60,-1),
                                            increment=0.1, min_val=0,max_val=50, agwStyle=FS.FS_RIGHT)
        self.ROuterConvergenceAngle.SetFormat("%f")
        self.ROuterConvergenceAngle.SetDigits(1)     
        
        wx.StaticText(pnl, label='RInnerConvergenceAngle', pos=(15, 445))
        self.RInnerConvergenceAngle = FS.FloatSpin(pnl, value=0.0, pos=(155, 440),size=(60,-1),
                                            increment=0.1, min_val=0,max_val=50, agwStyle=FS.FS_RIGHT)
        self.RInnerConvergenceAngle.SetFormat("%f")
        self.RInnerConvergenceAngle.SetDigits(1)
        
        wx.StaticText(pnl, label='RAcceleratingVoltage', pos=(15, 475))
        self.RAcceleratingVoltage = FS.FloatSpin(pnl, value=200.0, pos=(155, 470),size=(60,-1),
                                            increment=0.1, min_val=0,max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RAcceleratingVoltage.SetFormat("%f")
        self.RAcceleratingVoltage.SetDigits(1)
        
        wx.StaticText(pnl, label='IIncidentBeamDirectionX',pos=(220, 415))
        self.IIncidentBeamDirectionX = wx.SpinCtrl(pnl, pos=(380,410), size=(60,-1), value='1', min=-100000, max=100000)
        
        wx.StaticText(pnl, label='IIncidentBeamDirectionY',pos=(220, 445))
        self.IIncidentBeamDirectionY = wx.SpinCtrl(pnl, pos=(380,440), size=(60,-1), value='1', min=-100000, max=100000)
        
        wx.StaticText(pnl, label='IIncidentBeamDirectionZ',pos=(220, 475))
        self.IIncidentBeamDirectionZ = wx.SpinCtrl(pnl, pos=(380,470), size=(60,-1), value='1', min=-100000, max=100000)
        
        
        wx.StaticText(pnl, label='IXDirectionX',pos=(460, 415))
        self.IXDirectionX = wx.SpinCtrl(pnl, pos=(590,410), size=(60,-1), value='1', min=-100000, max=100000)
        
        wx.StaticText(pnl, label='IXDirectionY',pos=(460, 445))
        self.IXDirectionY = wx.SpinCtrl(pnl, pos=(590,440), size=(60,-1), value='-1', min=-100000, max=100000)
           
        wx.StaticText(pnl, label='IXDirectionZ',pos=(460, 475))
        self.IXDirectionZ = wx.SpinCtrl(pnl, pos=(590,470), size=(60,-1), value='1', min=-100000, max=100000)
        
        wx.StaticText(pnl, label='INormalDirectionX',pos=(460, 535))
        self.INormalDirectionX = wx.SpinCtrl(pnl, pos=(590,530), size=(60,-1), value='1', min=-100000, max=100000)
        
        wx.StaticText(pnl, label='INormalDirectionY',pos=(460, 565))
        self.INormalDirectionY = wx.SpinCtrl(pnl, pos=(590,560), size=(60,-1), value='1', min=-100000, max=100000)
        
        wx.StaticText(pnl, label='INormalDirectionZ',pos=(460, 595))
        self.INormalDirectionZ = wx.SpinCtrl(pnl, pos=(590,590), size=(60,-1), value='1', min=-100000, max=100000)
        
        
        #Image Output Options
        wx.StaticText(pnl, label='RInitialThickness', pos=(15, 555))
        self.RInitialThickness = FS.FloatSpin(pnl, value=1000.0, pos=(155, 550),size=(60,-1),
                                            increment=1, min_val=0, max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RInitialThickness.SetFormat("%f")
        self.RInitialThickness.SetDigits(1)
        
        wx.StaticText(pnl, label='RFinalThickness', pos=(15, 585))
        self.RFinalThickness = FS.FloatSpin(pnl, value=1000.0, pos=(155, 580),size=(60,-1),
                                            increment=1, min_val=0, max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RFinalThickness.SetFormat("%f")
        self.RFinalThickness.SetDigits(1)
        
        wx.StaticText(pnl, label='RDeltaThickness', pos=(220, 555))
        self.RDeltaThickness = FS.FloatSpin(pnl, value=10.0, pos=(340, 550),size=(60,-1),
                                            increment=1, min_val=0, max_val=100000, agwStyle=FS.FS_RIGHT)
        self.RDeltaThickness.SetFormat("%f")
        self.RDeltaThickness.SetDigits(1)
        
        wx.StaticText(pnl, label='IReflectOut',pos=(220, 585))
        self.IReflectOut = wx.SpinCtrl(pnl, pos=(340,580), size=(60,-1), value='7', min=0, max=100000)
        
        
        #Beam selection gui
        
        #textbox
        self.CIFtextbox=wx.TextCtrl(pnl, -1,os.getcwd(), pos=(200, 650), size=(250,-1),style=wx.TE_RIGHT)

        #Buttons
        Run = wx.Button(pnl, label='Run', pos=(5, 650), size=(60, -1))
        Cancel = wx.Button(pnl, label='Cancel', pos=(100,650), size=(60,-1))
        CIFFile=wx.Button(pnl, label='Browse', pos=(500, 650), size=(60,-1))
        InputFile = wx.Button(pnl, label='Write Input file', pos=(5,620), size=(150,-1))
   
                                  
        #CSC Checkbox
        self.CSC=wx.CheckBox(pnl,label='CSC',pos=(5,700))
        self.CSC.SetValue(True)
        
        #Number of cores
        wx.StaticText(pnl, label='MpiCores',pos=(150, 700))
        self.MPICores = wx.SpinCtrl(pnl, pos=(250,690), size=(60,-1), value='1', min=0, max=100)

        #the various functions of the buttons - run felix, write input file, cancel, and browse file
        Run.Bind(wx.EVT_BUTTON, self.CIFCreate)
        Cancel.Bind(wx.EVT_BUTTON, self.OnClose)
        CIFFile.Bind(wx.EVT_BUTTON, self.OnCif)
        InputFile.Bind(wx.EVT_BUTTON, self.InpCreate)

    

        self.SetSize((700, 730))
        self.SetTitle('Felix')
        self.Centre()
        self.Show(True)

   
    #subroutine which re-names the selected cif file and opens a directory for the input files - then runs Felix
    def CIFCreate(self,event):
        
        #opens new (re-existing directory and puts all input files into it-***ensure file and path check needs to be put here****
        cpath, cfilename = os.path.split(self.CIFtextbox.GetValue())

        if cfilename.find(".cif")==-1:
            wx.MessageBox('.cif file not found, please select a .cif file','Error', wx.OK | wx.ICON_ERROR)

        #creates working directory
        dir = cfilename.rstrip('.cif')+"_"+str(self.IMinReflectionPool.GetValue())+"_"+str(self.IMinStrongBeams.GetValue())+"_input_directory"
        if os.path.exists(dir):
            shutil.rmtree(dir)
        os.makedirs(dir)
        
        #copy cif file from user specified location to new directory
        shutil.copy2(self.CIFtextbox.GetValue(),dir+"/felix.cif")

        #copy sca file from samples folder - TO BE CHANGED, only for Richard's Use
        shutil.copy2("../samples/Si/felix.sca",dir+"/felix.sca")

        InputFileSwitch=1
        
        #Write input file
        self.WriteInputFile(dir,InputFileSwitch)
        wx.MessageBox('Files created successfully, click okay to run Felix', 'Info', wx.OK | wx.CANCEL | wx.ICON_INFORMATION)

        #Change to working directory
        os.chdir(dir)
        
        #Get Value of mpicores
        NumberofCores=self.MPICores.GetValue()
        
        print "directory name:",dir
        
        #Run in parallel or single core
        if NumberofCores == 1:
            os.system("../felixsim") #single core
        else:
            os.system("mpirun -n "+str(NumberofCores)+" ../felixsim") #parallel

        self.Close(True)

        
    def InpCreate(self,event):
        
        dir="Felix_input_file"

        if os.path.exists(dir):
            shutil.rmtree(dir)
        os.makedirs(dir) 

        InputFileCreate=WriteInputDialog(self,-1,'Save File As')
        InputFileBlock=InputFileCreate.ShowModal()

                
        InputFileCreate.Destroy()

        InputFileSwitch=2

        self.WriteInputFile(dir,InputFileSwitch)

        wx.MessageBox('Input File successfully written', 'Info', wx.OK | wx.ICON_INFORMATION)
        
        
    def WriteInputFile(self,dir,InputFileSwitch):
        
        FelixInpFilename= dir+"/felix.inp"
    
        inpfile = open(FelixInpFilename,"wb")

        
        #opening text in input file
        inpfile.write("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:"+"\n")
        inpfile.write("# ------------------------------------\n")
        inpfile.write("\n")
        inpfile.write("# ------------------------------------\n")
        inpfile.write("# felixsim input\n")
        inpfile.write("\n")
        inpfile.write("# control flags\n")

    

        #Write into the input file the input variables and their respective values
        ##########################################################################
        #control flags
        valIWriteFLAG = self.IWriteFLAG.GetValue()
        inpfile.write("IWriteFLAG                = ")
        inpfile.write(valIWriteFLAG)
        inpfile.write("\n")
        
        
        valIImageFLAG = self.IImageFLAG.GetValue()
        inpfile.write("IImageFLAG                = ")
        inpfile.write(valIImageFLAG)
        inpfile.write("\n")
        
        valIOutputFLAG = self.IOutputFLAG.GetValue()
        inpfile.write("IOutputFLAG               = ") 
        inpfile.write(valIOutputFLAG)
        inpfile.write("\n")        
        
        valIBinorTextFLAG = self.IBinorTextFLAG.GetValue()
        inpfile.write ("IBinorTextFLAG            = ")
        inpfile.write (valIBinorTextFLAG)
        inpfile.write ("\n")
        
        valIScatterFactorMethodFLAG = self.IScatterFactorMethodFLAG.GetValue()
        inpfile.write ("IScatterFactorMethodFLAG  = ") 
        inpfile.write (valIScatterFactorMethodFLAG)
        inpfile.write ("\n")
        
        valICentralBeamFLAG = self.ICentralBeamFLAG.GetValue()
        inpfile.write ("ICentralBeamFLAG          = ")
        inpfile.write (valICentralBeamFLAG)
        inpfile.write ("\n")
        
        valIMaskFLAG = self.IMaskFLAG.GetValue()
        inpfile.write ("IMaskFLAG                 = ") 
        inpfile.write (valIMaskFLAG)
        inpfile.write ("\n")

        valIZolzFLAG = self.IZolzFLAG.GetValue()
        inpfile.write ("IZolzFLAG                 = "+valIZolzFLAG+"\n")
        
        valIAbsorbFLAG = self.IAbsorbFLAG.GetValue()
        inpfile.write ("IAbsorbFLAG               = ") 
        inpfile.write (valIAbsorbFLAG)
        inpfile.write ("\n")
        
        valIAnisoDebyeWallerFLAG = self.IAnisoDebyeWallerFLAG.GetValue()
        inpfile.write ("IAnisoDebyeWallerFLAG     = ") 
        inpfile.write (valIAnisoDebyeWallerFLAG)
        inpfile.write ("\n")

        inpfile.write("IBeamConvergenceFLAG      = 1\n")
        
        valIPseudoCubicFLAG = self.IPseudoCubicFLAG.GetValue()
        inpfile.write ("IPseudoCubicFLAG          = ") 
        inpfile.write (valIPseudoCubicFLAG)
        inpfile.write ("\n")
        
        valIXDirectionFLAG = self.IXDirectionFLAG.GetValue()
        inpfile.write ("IXDirectionFLAG           = ") 
        inpfile.write (valIXDirectionFLAG)
        inpfile.write ("\n")

        #section break text
        inpfile.write("\n# radius of the beam in pixels\n")
        
        #radius of the beam in pixels
        #for the float_spin and spin_ctrl we need to convert from float to string
        valIPixelCount = self.IPixelCount.GetValue()
        inpfile.write ("IPixelCount               = ") 
        svalIPixelCount=str(int(valIPixelCount))
        inpfile.write (svalIPixelCount)
        inpfile.write ("\n")
        
        #section break
        inpfile.write("\n# beam selection criteria\n")
        
        valIMinReflectionPool = self.IMinReflectionPool.GetValue()
        inpfile.write ("IMinReflectionPool        = ")
        svalIMinReflectionPool=str(valIMinReflectionPool)
        inpfile.write (svalIMinReflectionPool)
        inpfile.write ("\n")
        
        valIMinStrongBeams = self.IMinStrongBeams.GetValue()
        inpfile.write ("IMinStrongBeams           = ")
        svalIMinStrongBeams=str(valIMinStrongBeams)
        inpfile.write (svalIMinStrongBeams)
        inpfile.write ("\n")
        
        valIMinWeakBeams = self.IMinWeakBeams.GetValue()
        inpfile.write ("IMinWeakBeams             = ")
        svalIMinWeakBeams=str(valIMinWeakBeams)
        inpfile.write (svalIMinWeakBeams)
        inpfile.write ("\n")
        
        valRBSBMax = self.RBSBMax.GetValue()
        inpfile.write ("RBSBMax                   = ")
        svalRBSBMax = str(valRBSBMax)
        inpfile.write (svalRBSBMax)
        inpfile.write("\n") 
        
        valRBSPMax = self.RBSPMax.GetValue()
        inpfile.write ("RBSPMax                   = ")
        svalRBSPMax=str(valRBSPMax)
        inpfile.write (svalRBSPMax)
        inpfile.write ("\n")
        
        valRConvergenceTolerance = self.RConvergenceTolerance.GetValue()
        inpfile.write ("RConvergenceTolerance (%) = ")
        svalRConvergenceTolerance=str(valRConvergenceTolerance)
        inpfile.write (svalRConvergenceTolerance)
        inpfile.write ("\n")

        #section break
        inpfile.write("\n# crystal settings\n")

        
        valRDebyeWallerConstant = self.RDebyeWallerConstant.GetValue()
        inpfile.write ("RDebyeWallerConstant      = ")
        svalRDebyeWallerConstant=str(valRDebyeWallerConstant)
        inpfile.write (svalRDebyeWallerConstant)
        inpfile.write ("\n")
        
        valRAbsorptionPer = self.RAbsorptionPer.GetValue()
        inpfile.write ("RAbsorptionPer            = " )
        svalRAbsorptionPer=str(valRAbsorptionPer)
        inpfile.write (svalRAbsorptionPer)
        inpfile.write ("\n")
        
        #section break
        inpfile.write("\n# microscope settings\n")

        
        valROuterConvergenceAngle = self.ROuterConvergenceAngle.GetValue()
        inpfile.write ("ROuterConvergenceAngle    = " )
        svalROuterConvergenceAngle=str(valROuterConvergenceAngle)
        inpfile.write (svalROuterConvergenceAngle)
        inpfile.write ("\n")
        
        valRInnerConvergenceAngle = self.RInnerConvergenceAngle.GetValue()
        inpfile.write ("RInnerConvergenceAngle    = ")
        svalRInnerConvergenceAngle=str(valRInnerConvergenceAngle)
        inpfile.write (svalRInnerConvergenceAngle)
        inpfile.write ("\n")
        
        valIIncidentBeamDirectionX = self.IIncidentBeamDirectionX.GetValue()
        inpfile.write ("IIncidentBeamDirectionX   = ")
        svalIIncidentBeamDirectionX=str(valIIncidentBeamDirectionX)
        inpfile.write (svalIIncidentBeamDirectionX)
        inpfile.write ("\n")
        
        valIIncidentBeamDirectionY = self.IIncidentBeamDirectionY.GetValue()
        inpfile.write ("IIncidentBeamDirectionY   = " )
        svalIIncidentBeamDirectionY=str(valIIncidentBeamDirectionY)
        inpfile.write (svalIIncidentBeamDirectionY)
        inpfile.write ("\n")
        
        valIIncidentBeamDirectionZ = self.IIncidentBeamDirectionZ.GetValue()
        inpfile.write ("IIncidentBeamDirectionZ   = " )
        svalIIncidentBeamDirectionZ=str(valIIncidentBeamDirectionZ)
        inpfile.write (svalIIncidentBeamDirectionZ)
        inpfile.write ("\n")
        
        valIXDirectionX = self.IXDirectionX.GetValue()
        inpfile.write ("IXDirectionX              = " )
        svalIXDirectionX=str(valIXDirectionX)
        inpfile.write (svalIXDirectionX)
        inpfile.write ("\n")
        
        valIXDirectionY = self.IXDirectionY.GetValue()
        inpfile.write ("IXDirectionY              = " )
        svalIXDirectionY=str(valIXDirectionY)
        inpfile.write (svalIXDirectionY)
        inpfile.write ("\n")
        
        valIXDirectionZ = self.IXDirectionZ.GetValue()
        inpfile.write ("IXDirectionZ              = " )
        svalIXDirectionZ=str(valIXDirectionZ)
        inpfile.write (svalIXDirectionZ)
        inpfile.write ("\n")
        
        valINormalDirectionX = self.INormalDirectionX.GetValue()
        inpfile.write ("INormalDirectionX         = " )
        svalINormalDirectionX=str(valINormalDirectionX)
        inpfile.write (svalINormalDirectionX)
        inpfile.write ("\n")
        
        valINormalDirectionY = self.INormalDirectionY.GetValue()
        inpfile.write ("INormalDirectionY         = " )
        svalINormalDirectionY=str(valINormalDirectionY)
        inpfile.write (svalINormalDirectionY)
        inpfile.write ("\n")
        
        valINormalDirectionZ = self.INormalDirectionZ.GetValue()
        inpfile.write ("INormalDirectionZ         = " )
        svalINormalDirectionZ=str(valINormalDirectionZ)
        inpfile.write (svalINormalDirectionZ)
        inpfile.write ("\n")
        
        valRAcceleratingVoltage = self.RAcceleratingVoltage.GetValue()
        inpfile.write ("RAcceleratingVoltage (kV) = ") 
        svalRAcceleratingVoltage=str(valRAcceleratingVoltage)
        inpfile.write (svalRAcceleratingVoltage)
        inpfile.write ("\n")
        
        #section break
        inpfile.write("\n# Image Output Options\n\n")

        
        valRInitialThickness = self.RInitialThickness.GetValue()
        inpfile.write ("RInitialThickness        = ")
        svalRInitialThickness=str(valRInitialThickness)
        inpfile.write (svalRInitialThickness)
        inpfile.write ("\n")
        
        valRFinalThickness = self.RFinalThickness.GetValue()
        inpfile.write ("RFinalThickness          = "  )
        svalRFinalThickness=str(valRFinalThickness)
        inpfile.write (svalRFinalThickness)
        inpfile.write ("\n")
        
        valRDeltaThickness = self.RDeltaThickness.GetValue()
        inpfile.write ("RDeltaThickness          = ")
        svalRDeltaThickness=str(valRDeltaThickness)
        inpfile.write (svalRDeltaThickness)
        inpfile.write ("\n")
        
        valIReflectOut = self.IReflectOut.GetValue()
        inpfile.write ("IReflectOut              = ") 
        svalIReflectOut=str(valIReflectOut)
        inpfile.write (svalIReflectOut)
        inpfile.write ("\n")
        ##########################################################################
        
        inpfile.close()
        
        if InputFileSwitch == 2:
            SInpFullDirectory= SInpPath+"/"+SInpFilename
            shutil.copy2(FelixInpFilename,SInpFullDirectory)
            print "here"
        

        #i = e.GetString()
        #self.IXDirection.SetLabel(i)       
        

    def OnClose(self, e):
        
        self.Close(True)    

    def OnCif(self,e):

        CIFFileDialog = wx.FileDialog(self, "Load CIF file", "", "",
                                       "CIF files (*.cif)|*.cif", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        #If the User selects cancel
        if CIFFileDialog.ShowModal() == wx.ID_CANCEL:
            return     

        # get the filename and path from user
        self.User_CIFfilename=CIFFileDialog.GetFilename()
        self.User_CIFpath=CIFFileDialog.GetPath()

        self.CIFtextbox.SetValue(self.User_CIFpath)

        
        
        #if not input_stream.IsOk():

        #   wx.LogError("Cannot open file '%s'."%openFileDialog.GetPath())
        #  return


        # def OnFileSelect(self,e):
        
        
def main():
    
    ex = wx.App()
    Felix(None)
    ex.MainLoop()    


if __name__ == '__main__':
    main()   
