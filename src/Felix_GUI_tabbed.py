import wx
import wx.lib.agw.floatspin as FS
import math
import os
import shutil

#Class defining gui which selects where the user wants to save the input file
#When Felix is run, this file will be copied as felix.inp
class WriteInputDialog(wx.Dialog):

    def __init__(self, parent,id, title):
        wx.Dialog.__init__(self,parent,id,title,size=(400,200))

    
        self.UserInterface2()
        #self.SetSize(400,200)
        #self.SetTitle("Save Input File As")[

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
            return()

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
            

class FlagPanel(wx.Panel):

    def __init__(self, parent):
        
        wx.Panel.__init__(self, parent)
        
        title = wx.StaticText(self, wx.ID_ANY, 'Flags') 
 
#=======================================================================================       
        # Lists for easily changing flags and choice options for wx.Choice
        flags = ['0','1','2','3']
        writeflags = ['0','1','2','3','4','10','100','101','102','103','104','110']
        flagnames = ['IWriteflag', 'IScatterFactorMethodflag', 'IAbsorbflag',
                        'IImageflag', 'IZolzflag', 'IAnisoDebyeWallerflag',
                        'IOutputflag', 'ICentralflagflag', 'IPseudoCubicflag',
                        'IBinorTextflag', 'Ben control']
        #How far into writeflags the choices need to be e.g. 2 = 0, 1; 4 = 0, 1, 2, 3
        flagChoices = [12, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2]
#=======================================================================================        
        
        self.flagnames = flagnames
        
        print flagnames
        
        # Making the code a bit more future proof by generating flag layout
        # from a list of flag names
        flagnumber = flagnames.__len__()
        print 'The number of flags: {}.\n'.format(flagnumber)
        numberOfRows = int(math.ceil(flagnumber / 3.0))
        print 'The number of rows: {}.\n'.format(numberOfRows)
        
        # Make some lists
        flagObjectsLabels   = []
        flagObjectsChoices  = []
        SizerObjects        = []
        
        self.flagObjectsChoices = flagObjectsChoices
        
        # Finds the number of empty slots on the bottom row to add spacer later
        spacerNo = (3 - (flagnumber % 3)) % 3
        print 'The number of spacers: {}.\n'.format(spacerNo)
        
        #Adds a sizer for each row to a sizer list
        for x in range(0, numberOfRows):
            SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))
        
        #Adds labels and choice objects to each respective lists
        for flag in flagnames:
            choices = flagChoices[flagnames.index(flag)]
            print(choices)
            flagObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, flag))
            flagObjectsChoices.append(wx.Choice(self, wx.ID_ANY, size=(60, -1), 
                    choices=writeflags[0:choices], name=flag)) 
        
        #Adds the objects from the lists to the respective rows in the sizer list
        for flagNo in range(0, flagnumber):
            row = int(math.floor(flagNo / 3))
            SizerObjects[row].Add(flagObjectsLabels[flagNo], 3, wx.ALL, 5)
            SizerObjects[row].Add(flagObjectsChoices[flagNo], 1, wx.ALL, 5)
        
        #Adds spacers if necessary    
        if spacerNo != 0:
            for x in range(0, spacerNo):
                SizerObjects[numberOfRows - 1].AddStretchSpacer(4)
        
        #Set up overall and title sizers    
        topflagSizer            = wx.BoxSizer(wx.VERTICAL)
        flagTitleSizer          = wx.BoxSizer(wx.HORIZONTAL)
        
        flagTitleSizer.Add(title, 0, wx.ALL, 5)
        
        topflagSizer.Add(flagTitleSizer, 0, wx.CENTER)
        
        #Add sizers from list to topflagSizer
        for sizerNo in range(0, numberOfRows):
            topflagSizer.Add(SizerObjects[sizerNo], 0, wx.CENTER)
        
        self.SetSizer(topflagSizer)
        topflagSizer.Fit(self)
        
class RadiusPanel(wx.Panel):

    def __init__(self, parent):
        
        wx.Panel.__init__(self, parent)
        
        title = wx.StaticText(self, wx.ID_ANY, 'Radius of Beam')
        RadiusLabel = wx.StaticText(self, wx.ID_ANY, 'Radius of Beam in Pixels')
        
        RadiusSizer            = wx.BoxSizer(wx.HORIZONTAL)
        RadiusTitleSizer       = wx.BoxSizer(wx.HORIZONTAL)
        topRadiusSizer         = wx.BoxSizer(wx.VERTICAL)
        
        self.IPixelCount = FS.FloatSpin(self, size=(60, -1), value=64, min_val=0, max_val=512,
                                   increment=64, agwStyle=FS.FS_RIGHT)
        self.IPixelCount.SetFormat("%f")
        self.IPixelCount.SetDigits(0)
        RadiusSizer.Add(RadiusLabel, 3, wx.ALL, 5)
        RadiusSizer.Add(IPixelCount, 1, wx.ALL, 5)
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
        
#===============================================================================
        
        # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
        # ['name', default value, increment, min, mac, 'type']
        # with type referring to a 1 for spinctrl or 2 for a float spin!
        # NB: spin does not need increment, so just put 0!
        self.BeamControl1 = ['IMinReflectionPool', '15', 0, 0, 100000, 1]
        self.BeamControl2 = ['IMinStrongBeams', '7', 0, 0, 100000, 1]
        self.BeamControl3 = ['IMinWeakBeams', '5', 0, 0, 100000, 1]
        self.BeamControl4 = ['RBSBMax', '0.1', 0.1, 0, 100000, 2]
        self.BeamControl5 = ['RBSPMax', '0.1', 0.1, 0, 100000, 2]
        self.BeamControl6 = ['RConvergenceTolerance', '0.1', 0.1, 0, 100000, 2]
        
        #Add them to list (of lists) - Need to find a better method for this
        BeamControlList.append(self.BeamControl1)
        BeamControlList.append(self.BeamControl2)
        BeamControlList.append(self.BeamControl3)
        BeamControlList.append(self.BeamControl4)
        BeamControlList.append(self.BeamControl5)
        BeamControlList.append(self.BeamControl6)
        
#===============================================================================
        
        
        title = wx.StaticText(self, wx.ID_ANY, 'Beam Selection')
        
        # Making the code a bit more future proof by generating beam layout
        # from a list of beam names
        beamnumber = BeamControlList.__len__()
        print 'The number of beam controls: {}.\n'.format(beamnumber)
        numberOfRows = int(math.ceil(beamnumber / 3.0))
        print 'The number of rows: {}.\n'.format(numberOfRows)
        
        # Make some lists
        SizerObjects        = []
        beamObjectsLabels    = []
        beamObjectsControls  = []
        
        self.SizerObjects = SizerObjects
        self.beamObjectsLabels = beamObjectsLabels
        self.beamObjectsControls = beamObjectsControls
        
        # Finds the number of empty slots on the bottom row to add spacer later
        spacerNo = (3 - (beamnumber % 3)) % 3
        print 'The number of spacers: {}.\n'.format(spacerNo)
        
        #Adds a sizer for each row to a sizer list
        for x in range(0, numberOfRows):
            SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))
        
        #Adds labels and choice objects to each respective lists
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
            
        #Adds the objects from the lists to the respective rows in the sizer list
        for beamNo in range(0, beamnumber):
            row = int(math.floor(beamNo / 3))
            SizerObjects[row].Add(beamObjectsLabels[beamNo], 3, wx.ALL, 5)
            SizerObjects[row].Add(beamObjectsControls[beamNo], 1, wx.ALL, 5)
        
        #Adds spacers if necessary    
        if spacerNo != 0:
            for x in range(0, spacerNo):
                SizerObjects[numberOfRows - 1].AddStretchSpacer(4)
        
        #Set up overall and title sizers    
        topbeamSizer            = wx.BoxSizer(wx.VERTICAL)
        beamTitleSizer          = wx.BoxSizer(wx.HORIZONTAL)
        
        beamTitleSizer.Add(title, 0, wx.ALL, 5)
        
        topbeamSizer.Add(beamTitleSizer, 0, wx.CENTER)
        
        #Add sizers from list to topbeamSizer
        for sizerNo in range(0, numberOfRows):
            topbeamSizer.Add(SizerObjects[sizerNo], 0)
        
        self.SetSizer(topbeamSizer)
        topbeamSizer.Fit(self)

class crystalPanel(wx.Panel):


    def __init__(self, parent):
        
        wx.Panel.__init__(self, parent)
        
        crystalControlList = []
        
#===============================================================================
        
        # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
        # ['name', default value, increment, min, mac, 'type', NUMBER OF DIGITS]
        # with type referring to a 1 for spinctrl or 2 for a float spin!
        # NB: spin does not need increment, so just put 0!
        crystalControl1 = ['RDebyeWallerConstant', '0.467', 0.001, 0, 100000, 2, 3]
        crystalControl2 = ['RAbsorptionPer', '2.9', 0.1, 0, 100000, 2, 1]
        
        
        #Add them to list (of lists) - Need to find a better method for this
        crystalControlList.append(crystalControl1)
        crystalControlList.append(crystalControl2)
#===============================================================================
        
        
        title = wx.StaticText(self, wx.ID_ANY, 'Crystal Settings')
        
        # Making the code a bit more future proof by generating crystal layout
        # from a list of crystal names
        crystalnumber = crystalControlList.__len__()
        print 'The number of crystal controls: {}.\n'.format(crystalnumber)
        numberOfRows = int(math.ceil(crystalnumber / 3.0))
        print 'The number of rows: {}.\n'.format(numberOfRows)
        
        # Make some lists
        SizerObjects        = []
        crystalObjectsLabels    = []
        crystalObjectsControls  = []
        
        # Finds the number of empty slots on the bottom row to add spacer later
        spacerNo = (3 - (crystalnumber % 3)) % 3
        print 'The number of spacers: {}.\n'.format(spacerNo)
        
        #Adds a sizer for each row to a sizer list
        for x in range(0, numberOfRows):
            SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))
        
        #Adds labels and choice objects to each respective lists
        for crystal in crystalControlList:
            crystalname = crystal[0]
            crystalvalue = crystal[1]
            crystalincrement = crystal[2]
            crystalmin = crystal[3]
            crystalmax = crystal[4]
            crystaltype = crystal[5]
            crystaldigits = crystal[6]
            
            currentIndex = crystalControlList.index(crystal)
            
            crystalObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, crystalname))
            
            if crystaltype == 1:
                crystalObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY, 
                    size=(60, -1), value=crystalvalue, min=crystalmin, max=crystalmax))
                print('Added a spin!')
            elif crystaltype == 2:
                crystalObjectsControls.append(FS.FloatSpin(self, size=(60, -1), 
                    value=crystalvalue, increment=crystalincrement, 
                    min_val=crystalmin, max_val=crystalmax, agwStyle=FS.FS_RIGHT))
                
                crystalObjectsControls[currentIndex].SetFormat("%f")
                crystalObjectsControls[currentIndex].SetDigits(crystaldigits)
                print('Added a float spin!')
            
        #Adds the objects from the lists to the respective rows in the sizer list
        for crystalNo in range(0, crystalnumber):
            row = int(math.floor(crystalNo / 3))
            SizerObjects[row].Add(crystalObjectsLabels[crystalNo], 3, wx.ALL, 5)
            SizerObjects[row].Add(crystalObjectsControls[crystalNo], 1, wx.ALL, 5)
        
        #Adds spacers if necessary    
        if spacerNo != 0:
            for x in range(0, spacerNo):
                SizerObjects[numberOfRows - 1].AddStretchSpacer(4)
        
        #Set up overall and title sizers    
        topcrystalSizer            = wx.BoxSizer(wx.VERTICAL)
        crystalTitleSizer          = wx.BoxSizer(wx.HORIZONTAL)
        
        crystalTitleSizer.Add(title, 0, wx.ALL, 5)
        
        topcrystalSizer.Add(crystalTitleSizer, 0, wx.CENTER)
        
        #Add sizers from list to topcrystalSizer
        for sizerNo in range(0, numberOfRows):
            topcrystalSizer.Add(SizerObjects[sizerNo], 0)
        
        self.SetSizer(topcrystalSizer)
        topcrystalSizer.Fit(self)

class microscopePanel(wx.Panel):


    def __init__(self, parent):
        
        wx.Panel.__init__(self, parent)
        
        microscopeControlList = []
        
#===============================================================================
        
        # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
        # ['name', default value, increment, min, mac, 'type']
        # with type referring to a 1 for spinctrl or 2 for a float spin!
        # NB: spin does not need increment, so just put 0!
        microscopeControl1 = ['ROuterConvergenceAngle', '3.0', 0.1, 0, 50, 2]
        microscopeControl2 = ['RInnerConvergenceAngle', '0.0', 0.1, 0, 50, 2]
        microscopeControl3 = ['RAcceleratingVoltage', '200.0', 0.1, 0, 100000, 2]
        microscopeControl4 = ['IIncidentBeamDirectionX', '1', 0, -100000, 100000, 1]
        microscopeControl5 = ['IIncidentBeamDirectionY', '1', 0, -100000, 100000, 1]
        microscopeControl6 = ['IIncidentBeamDirectionZ', '1', 0, -100000, 100000, 1]
        microscopeControl7 = ['IXDirectionX', '1', 0, -100000, 100000, 1]
        microscopeControl8 = ['IXDirectionY', '1', 0, -100000, 100000, 1]
        microscopeControl9 = ['IXDirectionZ', '1', 0, -100000, 100000, 1]
        microscopeControl10 = ['INormalDirectionX', '1', 0, -100000, 100000, 1]
        microscopeControl11 = ['INormalDirectionY', '1', 0, -100000, 100000, 1]
        microscopeControl12 = ['INormalDirectionZ', '1', 0, -100000, 100000, 1]
        
        
        
        #Add them to list (of lists) - Need to find a better method for this
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
        
#===============================================================================
        
        title = wx.StaticText(self, wx.ID_ANY, 'Microscope Selection')
        
        # Making the code a bit more future proof by generating microscope layout
        # from a list of microscope names
        microscopenumber = microscopeControlList.__len__()
        print 'The number of microscope controls: {}.\n'.format(microscopenumber)
        numberOfRows = int(math.ceil(microscopenumber / 3.0))
        print 'The number of rows: {}.\n'.format(numberOfRows)
        
        # Make some lists
        SizerObjects        = []
        microscopeObjectsLabels    = []
        microscopeObjectsControls  = []
        
        # Finds the number of empty slots on the bottom row to add spacer later
        spacerNo = (3 - (microscopenumber % 3)) % 3
        print 'The number of spacers: {}.\n'.format(spacerNo)
        
        #Adds a sizer for each row to a sizer list
        for x in range(0, numberOfRows):
            SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))
        
        #Adds labels and choice objects to each respective lists
        for microscope in microscopeControlList:
            microscopename = microscope[0]
            microscopevalue = microscope[1]
            microscopeincrement = microscope[2]
            microscopemin = microscope[3]
            microscopemax = microscope[4]
            microscopetype = microscope[5]
            
            currentIndex = microscopeControlList.index(microscope)
            
            microscopeObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, microscopename))
            
            if microscopetype == 1:
                microscopeObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY, 
                        size=(60, -1), value=microscopevalue, 
                        min=microscopemin, max=microscopemax))
                print('Added a spin!')
            elif microscopetype == 2:
                microscopeObjectsControls.append(FS.FloatSpin(self, size=(60, -1), 
                    value=microscopevalue, increment=microscopeincrement, 
                    min_val=microscopemin, max_val=microscopemax, agwStyle=FS.FS_RIGHT))
                
                microscopeObjectsControls[currentIndex].SetFormat("%f")
                microscopeObjectsControls[currentIndex].SetDigits(1)
                print('Added a float spin!')
            
        #Adds the objects from the lists to the respective rows in the sizer list
        for microscopeNo in range(0, microscopenumber):
            row = int(math.floor(microscopeNo / 3))
            SizerObjects[row].Add(microscopeObjectsLabels[microscopeNo], 3, wx.ALL, 5)
            SizerObjects[row].Add(microscopeObjectsControls[microscopeNo], 1, wx.ALL, 5)
        
        #Adds spacers if necessary    
        if spacerNo != 0:
            for x in range(0, spacerNo):
                SizerObjects[numberOfRows - 1].AddStretchSpacer(4)
        
        #Set up overall and title sizers    
        topmicroscopeSizer            = wx.BoxSizer(wx.VERTICAL)
        microscopeTitleSizer          = wx.BoxSizer(wx.HORIZONTAL)
        
        microscopeTitleSizer.Add(title, 0, wx.ALL, 5)
        
        topmicroscopeSizer.Add(microscopeTitleSizer, 0, wx.CENTER)
        
        #Add sizers from list to topmicroscopeSizer
        for sizerNo in range(0, numberOfRows):
            topmicroscopeSizer.Add(SizerObjects[sizerNo], 0)
        
        self.SetSizer(topmicroscopeSizer)
        topmicroscopeSizer.Fit(self)

class imagePanel(wx.Panel):


    def __init__(self, parent):
        
        wx.Panel.__init__(self, parent)
        
        imageControlList = []
        
#===============================================================================
        
        # SET UP ALL THE CONTROLS! THIS IS THE FORMAT:
        # ['name', default value, increment, min, mac, 'type', NUMBER OF DIGITS]
        # with type referring to a 1 for spinctrl or 2 for a float spin!
        # NB: spin does not need increment, so just put 0!
        imageControl1 = ['RInitialThickness', '1000.0', 1, 0, 100000, 2, 1]
        imageControl2 = ['RFinalThickness', '1000.0', 1, 0, 100000, 2, 1]
        imageControl3 = ['RDeltaThickness', '10.0', 1, 0, 100000, 2, 1]
        imageControl4 = ['IReflectOut', '7', 0, 0, 100000, 1, 1]
        
        
        
        #Add them to list (of lists) - Need to find a better method for this
        imageControlList.append(imageControl1)
        imageControlList.append(imageControl2)
        imageControlList.append(imageControl3)
        imageControlList.append(imageControl4)
#===============================================================================
        
        
        title = wx.StaticText(self, wx.ID_ANY, 'image Settings')
        
        # Making the code a bit more future proof by generating image layout
        # from a list of image names
        imagenumber = imageControlList.__len__()
        print 'The number of image controls: {}.\n'.format(imagenumber)
        numberOfRows = int(math.ceil(imagenumber / 3.0))
        print 'The number of rows: {}.\n'.format(numberOfRows)
        
        # Make some lists
        SizerObjects        = []
        imageObjectsLabels    = []
        imageObjectsControls  = []
        
        # Finds the number of empty slots on the bottom row to add spacer later
        spacerNo = (3 - (imagenumber % 3)) % 3
        print 'The number of spacers: {}.\n'.format(spacerNo)
        
        #Adds a sizer for each row to a sizer list
        for x in range(0, numberOfRows):
            SizerObjects.append(wx.BoxSizer(wx.HORIZONTAL))
        
        #Adds labels and choice objects to each respective lists
        for image in imageControlList:
            imagename = image[0]
            imagevalue = image[1]
            imageincrement = image[2]
            imagemin = image[3]
            imagemax = image[4]
            imagetype = image[5]
            imagedigits = image[6]
            
            currentIndex = imageControlList.index(image)
            
            imageObjectsLabels.append(wx.StaticText(self, wx.ID_ANY, imagename))
            
            if imagetype == 1:
                imageObjectsControls.append(wx.SpinCtrl(self, id=wx.ID_ANY, 
                        size=(60, -1), value=imagevalue, min=imagemin, max=imagemax))
                print('Added a spin!')
            elif imagetype == 2:
                imageObjectsControls.append(FS.FloatSpin(self, size=(60, -1), 
                        value=imagevalue, increment=imageincrement, min_val=imagemin, 
                        max_val=imagemax, agwStyle=FS.FS_RIGHT))
                
                imageObjectsControls[currentIndex].SetFormat("%f")
                imageObjectsControls[currentIndex].SetDigits(imagedigits)
                print('Added a float spin!')
            
        #Adds the objects from the lists to the respective rows in the sizer list
        for imageNo in range(0, imagenumber):
            row = int(math.floor(imageNo / 3))
            SizerObjects[row].Add(imageObjectsLabels[imageNo], 3, wx.ALL, 5)
            SizerObjects[row].Add(imageObjectsControls[imageNo], 1, wx.ALL, 5)
        
        #Adds spacers if necessary    
        if spacerNo != 0:
            for x in range(0, spacerNo):
                SizerObjects[numberOfRows - 1].AddStretchSpacer(4)
        
        #Set up overall and title sizers    
        topimageSizer            = wx.BoxSizer(wx.VERTICAL)
        imageTitleSizer          = wx.BoxSizer(wx.HORIZONTAL)
        
        imageTitleSizer.Add(title, 0, wx.ALL, 5)
        
        topimageSizer.Add(imageTitleSizer, 0, wx.CENTER)
        
        #Add sizers from list to topimageSizer
        for sizerNo in range(0, numberOfRows):
            topimageSizer.Add(SizerObjects[sizerNo], 0)
        
        self.SetSizer(topimageSizer)
        topimageSizer.Fit(self)

class optionPanel(wx.Panel):
        
    def __init__(self, parent, main):
        wx.Panel.__init__(self, parent)
        
        
        #Probably should change this at some point
        self.main = main
        
        title = wx.StaticText(self, wx.ID_ANY, 'Options')
    
        #optionSizer            = wx.BoxSizer(wx.HORIZONTAL)
        optionTitleSizer       = wx.BoxSizer(wx.HORIZONTAL)
        topOptionSizer         = wx.BoxSizer(wx.VERTICAL)
        
    
    
    
        optionTitleSizer.Add(title, 0, wx.ALL, 5)
    
        #Text box
        CIFtextbox = wx.TextCtrl(self, -1, os.getcwd(), size=(400,-1), style=wx.TE_RIGHT)
        textSizer  = wx.BoxSizer(wx.HORIZONTAL)
        textSizer.Add(CIFtextbox, 4, wx.ALL, 5)
        #textSizer.AddStretchSpacer(4)

        #Buttons
        Run = wx.Button(self, label='Run')
        Cancel = wx.Button(self, label='Cancel')
        CIFFile=wx.Button(self, label='Browse')
        InputFile = wx.Button(self, label='Write Input file')
        
        buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonSizer.Add(Run, 0, wx.ALL, 5)
        buttonSizer.Add(Cancel, 0, wx.ALL, 5)
        buttonSizer.Add(CIFFile, 0, wx.ALL, 5)
        buttonSizer.Add(InputFile, 0, wx.ALL, 5)
                              
        #CSC Checkbox
        CSC=wx.CheckBox(self,label='CSC')
        CSC.SetValue(True)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)
        checkSizer.Add(CSC, 0, wx.ALL, 5)
    
        #Number of cores
        coreLabel = wx.StaticText(self, label='MpiCores')
        MPICores = wx.SpinCtrl(self, size=(60,-1), value='1', min=0, max=100)
        coreSizer = wx.BoxSizer(wx.HORIZONTAL)
        coreSizer.Add(coreLabel, 3, wx.ALL, 5)
        coreSizer.Add(MPICores, 1, wx.ALL, 5)
        #coreSizer.AddStretchSpacer(4)
        
        topOptionSizer.Add(optionTitleSizer, 0, wx.CENTER)
        topOptionSizer.Add(textSizer, 0, wx.CENTER)
        topOptionSizer.Add(buttonSizer, 0, wx.CENTER)
        topOptionSizer.Add(checkSizer, 0, wx.CENTER)
        topOptionSizer.Add(coreSizer, 0, wx.CENTER)
        
        #topOptionSizer.Add(optionSizer, 0)
        self.SetSizer(topOptionSizer)
        topOptionSizer.Fit(self)
        
        #the various functions of the buttons - run felix, write input file, cancel, and browse file
        Run.Bind(wx.EVT_BUTTON, self.CIFCreate)
        Cancel.Bind(wx.EVT_BUTTON, self.OnClose)
        CIFFile.Bind(wx.EVT_BUTTON, self.OnCif)
        InputFile.Bind(wx.EVT_BUTTON, self.InpCreate)

        
    #subroutine which re-names the selected cif file and opens a directory 
    #for the input files - then runs Felix
    def CIFCreate(self, event):
        
        #opens new (re-existing) directory and puts all input files into it
        #***ensure file and path check needs to be put here****
        cpath, cfilename = os.path.split(self.CIFtextbox.GetValue())

        if cfilename.find(".cif")==-1:
            wx.MessageBox('.cif file not found, please select a .cif file','Error', 
                        wx.OK | wx.ICON_ERROR)

        #creates working directory
        cfilename=cfilename.replace(" ","")
        print cfilename
        dir = cfilename.rstrip('.cif')+"_"+str(beamPanel.IMinReflectionPool.GetValue())\
                    +"_"+str(beamPanel.IMinStrongBeams.GetValue())+"_input_directory"
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
        wx.MessageBox('Files created successfully, click okay to run Felix', 'Info', \
                    wx.OK | wx.CANCEL | wx.ICON_INFORMATION)

        #Change to working directory
        os.chdir(dir)
        
        #Get Value of mpicores
        NumberofCores=self.MPICores.GetValue()
        
        
        #Run in parallel or single core
        if NumberofCores == 1:
            os.system("../felixsim") #single core
        else:
            os.system("mpirun -n "+str(NumberofCores)+" ../felixsim") #parallel

        self.Close(True)
    
    def InpCreate(self, event):
        
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
        
        for flag in self.main.notebook.page1.flagnames:
            index = self.main.notebook.page1.flagnames.index(flag)
            
            value = self.main.notebook.page1.flagObjectsChoices[index].GetCurrentSelection()
            inpfile.write(flag)
            inpfile.write("                = ")
            inpfile.write(str(value))
            inpfile.write("\n")
            
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

        optionsPanel.CIFtextbox.SetValue(self.User_CIFpath)

        
        
        #if not input_stream.IsOk():

        #   wx.LogError("Cannot open file '%s'."%openFileDialog.GetPath())
        #  return


        # def OnFileSelect(self,e):    

        
              
class Notebook(wx.Notebook):

    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, wx.ID_ANY)
        
        #Add pages to notebook for different tabs
        self.page1 = FlagPanel(self)
        self.page2 = RadiusPanel(self)
        self.page3 = BeamPanel(self)
        self.page4 = crystalPanel(self)
        self.page5 = microscopePanel(self)
        self.page6 = imagePanel(self)
        #page7 = optionPanel(self)
        
        self.AddPage(self.page1, "Flags")
        self.AddPage(self.page2, "Radius")
        self.AddPage(self.page3, "Beam")
        self.AddPage(self.page4, "Crystal")
        self.AddPage(self.page5, "Microscope")
        self.AddPage(self.page6, "Image")
        #self.AddPage(page7, "Options")

        
class MainFrame(wx.Frame):

    def __init__(self):
        
        wx.Frame.__init__(self, None, title="Felix")
        
		#Create a panel with a notebook on it
        panel = wx.Panel(self)
        self.notebook = Notebook(panel)
        option = optionPanel(panel, self)
                
        sizer = wx.BoxSizer(wx.VERTICAL)
 
        # add the widgets to the sizers
        sizer.Add(self.notebook, 0, wx.ALL, 5)
        sizer.Add(option, 0, wx.ALL|wx.CENTER, 5)
 
        panel.SetSizer(sizer)
        sizer.Fit(self)
 
        self.Show()
        
        
        
        
if __name__ == '__main__':
    app = wx.App()
    frame = MainFrame().Show()
    app.MainLoop()
        