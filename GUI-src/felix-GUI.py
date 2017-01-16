#!/usr/bin/env python2

import wx
import wx.lib.agw.floatspin as FS
import wx.lib.scrolledpanel
import math
import itertools
import shutil


class input_var(): #input variable
    
    def __init__(self, name):
        self.name = name
        #check name isinstance(name, str)
        
class choice_var(input_var): #choice variable - select option from list

    def __init__(self, name, choices, defaultRef=0):
        #all(isinstance(elem, str) for elem in choices)
        #isinstance(defaultRef)
        input_var.__init__(self, name)
        #preferably directly assign these to inputwidget
        self.choices = choices
        self.defaultRef = defaultRef
        
    def inputWidget(self, panel): 
        
        def Chosen(event):
            self.chosenRef = self.inputWidget.GetSelection()
            
        self.inputWidget = wx.Choice(panel, wx.ID_ANY, size=(100, -1),
                                     choices=self.choices, name=self.name)
        self.inputWidget.SetStringSelection(self.choices[self.defaultRef])
        self.inputWidget.Bind(wx.EVT_CHOICE, Chosen)    
              
        return self.inputWidget
        
    def writeInputLine(self):
        i = 30 - len(self.name)
        j = 35 - len(self.choices[self.chosenRef])
        return self.name + ':' + ' '*i + str(self.chosenRef) + ' '*j + \
                self.choices[self.chosenRef]

    def read(self, choiceRef):
        try: 
            choiceRef = int(choiceRef)
        except:
            return False
        else:
            if 0 <= choiceRef <= len(self.choices): 
                self.inputWidget.SetSelection(choiceRef)
                return True                
            else:
                return False

class value_var(input_var):
    
    def __init__(self, name):
        input_var.__init__(self, name)
        
    def inputWidget(self, panel):
        self.inputWidget = wx.TextCtrl(panel, wx.ID_ANY)
        
#==============================================================================
#         def checkChar(self, event):
#             keycode = event.GetKeyCode()
#             if keycode in list(range(48,58))+[46]:
#                 if(keycode != 46 or '.' not in self.box.GetValue()):
#                     event.Skip()                    
#         self.inputwidget.Bind(wx.EVT_CHAR, checkChar)
#==============================================================================
        return self.inputWidget
        #return wx.TextCtrl(panel, wx.ID_ANY)
                    
    def writeInputLine(self):
        i = 30 - len(self.name)
        return self.name + ':' + ' '*i + str(self.chosen)
        
    def read(self, value):
        try:
            value = float(value)
            self.inputWidget.SetValue(value)
            return True
        except ValueError:
            return False
            
class value_var2(input_var):
    
    def __init__(self, name, default, min_val, max_val, increment):
        input_var.__init__(self, name, default)
        self.increment = increment
        self.min_val = min_val
        self.max_val = max_val
        self.chosen = default
        
    def inputWidget(self, panel):
        
        def Chosen(event):
            self.chosen = self.inputWidget.GetValue()
        
        self.inputWidget = FS.FloatSpin(panel, size=(140, -1), value=self.default, 
                        min_val=self.min_val, max_val=self.max_val, 
                        increment=self.increment, agwStyle=FS.FS_RIGHT)        
        self.inputWidget.SetDigits(1)
        self.inputWidget.Bind(FS.EVT_FLOATSPIN, Chosen)        
        return self.inputWidget
        
    def writeInputLine(self):
        i = 30 - len(self.name)
        return self.name + ':' + ' '*i + str(self.chosen)
        
    def read(self, value):
        try:
            value = float(value)
            if self.min_val <= value <= self.max_val: 
                self.inputWidget.SetValue(value)
                return True
            else:
                return False
        except ValueError:
            return False

class combination_var(choice_var):
    
    def __init__(self, name, choices, defaultRef=0):
        #all(isinstance(elem, str) for elem in choices)
        comboChoices = []
        for i in range(1,len(choices)+1):
            combinations =  list(itertools.combinations(choices,i)) 
            for i in range(len(combinations)):
                temp_list = list(combinations[i])
                comboChoice = temp_list[0]
                for i in range(1, len(temp_list)):
                    comboChoice += ' & '
                    comboChoice += temp_list[i]
                comboChoices.append(comboChoice )
        choice_var.__init__(self, name, comboChoices, defaultRef)

        
iv = [] #INPUT VARIABLES      
iv.extend([
    choice_var('IWriteFLAG', ['0','5','11']),    
    combination_var('IImageFLAG', ['Diffractions','Reflections',
                    'Amplitude+Phase']),
    combination_var('IOutputFLAG', ['Nothing','Ug Matrix','Eigenspectra',
                    'Wavefunctions'], 0),
    choice_var('IBinorTextFLAG', ['Binary','Text'], 1),
    choice_var('IScatteringFactorMethod', ['Kirklands','Peng',
                    'Doyle+Turner'], 0),
    choice_var('IZolzFLAG', ['No','Yes'], 0),
    choice_var('ICentralBeamFlag', ['Not included','included'], 0),
    choice_var('IMaskFLAG', ['Circular','Square'], 0),
    choice_var('IAbsorbFLAG', ['No','Proprtional','Einstein Perturbative',
                    'Einstein Exact'], 0),
    choice_var('IAnisotropicDebyeWallerFLAG', ['No','Yes'], 0),
    choice_var('IPseudoCubicFLAG', ['Orthorhombic','PseudoCubic'], 0),
    choice_var('IXDirectionFLAG', ['Ignore','Use'], 0),
    value_var('IPixelCount'),
    value_var('IMinReflectionPool'),
    value_var('IMinStrongBeams'),
    value_var('IMinWeakBeams'),
    value_var('RBSBmax'),
    value_var('RBSPmax'),
    value_var('RDebyeWallerConstant'),
    value_var('RAbsorptionPer'),
    value_var('RConvergenceAngle'),
    value_var('IIncidentBeamDirectionX'),
    value_var('IIncidentBeamDirectionY'),
    value_var('IIncidentBeamDirectionZ'),
    value_var('IXDirectionX'),
    value_var('IXDirectionY'),
    value_var('IXDirectionZ'),
    value_var('INormalDirectionX'),
    value_var('INormalDirectionY'),
    value_var('INormalDirectionZ'),
    value_var('RAccelerationVoltage'),
    value_var('RInitialThickness'),
    value_var('RFinalThickness'),
    value_var('RDeltaThickness'),
    value_var('IReflectOut'),
    value_var('IImageOutputFLAG'),
    value_var('IRefineModeFLAG'),
    ])

  
print list(range(len(iv)))

class MainFrame(wx.Frame):

    def __init__(self):
      
        wx.Frame.__init__(self, None, title="Felix")
        #outerPanel = wx.Panel(self)
        #panel = wx.Panel(outerPanel)
        panel = wx.lib.scrolledpanel.ScrolledPanel(self, size=(1000,100))
        #with sizer this x scroll shouldn't be necessary
        panel.SetupScrolling(scroll_x=False)
        title = wx.StaticText(self, wx.ID_ANY, 'variables')
        
        numberRows = int(math.ceil(len(iv) / 3.0)) ## three on a row
        print(numberRows)
        
        
        # Makes a list of sizers for each row
        rowSizer = []
        for x in range(numberRows):
            rowSizer.append(wx.BoxSizer(wx.HORIZONTAL))
                      
        # Adds the objects from the lists to the respective rows in the sizer list
        for i in range(len(iv)):
            row = int(math.floor(i / 3.0))
            print(row)
            variableLabel = wx.StaticText(panel, wx.ID_ANY, iv[i].name)
            #box = wx.TextCtrl(panel, wx.ID_ANY)
            rowSizer[row].Add(variableLabel, 4, wx.ALL, 5)
            #rowSizer[row].Add(box, 2, wx.ALL, 5)
            rowSizer[row].Add(iv[i].inputWidget(panel), 2, wx.ALL, 5)

        
        # Adds spacers if necessary
        # Finds the number empty slots on the bottom row to add spacer later
        spacerNo = (3 - (len(iv) % 3)) % 3
        if spacerNo != 0:
            for x in range(spacerNo):
                rowSizer[-1].AddStretchSpacer(6)
        
        printButton = wx.Button(self, label = 'print')
        printButton.Bind(wx.EVT_BUTTON, self.onPrint)
        createINPButton = wx.Button(self, label = 'create .inp')
        createINPButton.Bind(wx.EVT_BUTTON, self.onCreateINP)
        loadINPButton = wx.Button(self, label = 'load from .inp')
        loadINPButton.Bind(wx.EVT_BUTTON, self.onLoadINP)
        tempButton = wx.Button(self, label = 'temp')
        tempButton.Bind(wx.EVT_BUTTON, self.onTemp)
        
          
        # Set up overall and title sizers
        panelSizer = wx.BoxSizer(wx.VERTICAL)
        #panelSizer2 = wx.BoxSizer(wx.VERTICAL)
        titleSizer = wx.BoxSizer(wx.HORIZONTAL)
        ButtonsSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        titleSizer.Add(title, 0, wx.ALL, 5)
        ButtonsSizer.Add(printButton, 0, wx.ALL, 5)
        ButtonsSizer.Add(createINPButton, 0, wx.ALL, 5)
        ButtonsSizer.Add(loadINPButton, 0, wx.ALL, 5)
        ButtonsSizer.Add(tempButton, 0, wx.ALL, 5)
        
        # Add list of row sizers (vertically) to panelSizer
        for i in range(len(rowSizer)):
            panelSizer.Add(rowSizer[i], 0, wx.ALL, 2)
            print 'row added, ', i
            #if i < 2:
                #panelSizer2.Add(rowSizer[i], 0, wx.ALL, 2)
        
        frameSizer = wx.BoxSizer(wx.VERTICAL)
        frameSizer.Add(titleSizer, 0, wx.CENTRE)            
        frameSizer.Add(panel, 1, wx.CENTER)
        frameSizer.Add(wx.StaticLine(self), 0, wx.ALL|wx.EXPAND, 2)
        frameSizer.Add(ButtonsSizer, 0, wx.CENTRE) 
        
        panel.SetSizer(panelSizer)
        panelSizer.Fit(panel)          
        #outerPanel.SetSizer(frameSizer)
        self.SetSizer(frameSizer)
        self.SetMinSize((1000,150))
        frameSizer.Fit(self)
        


        
    def onCreateINP(self, event):
        dlg = wx.DirDialog(self)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print "felix.inp was created in this path: %s" % path
            
            
            inpFilename = path + "/felix.inp"
            txtFilename = path + "/felix.txt"
        
            inpfile = open(inpFilename, "wb")
            
            inpfile.write(
                          "# Input file for felixsim/draw/refine "
                              "version:VERSION: Build:BUILD:"
                          "\n"
                          "# ------------------------------------\n"
                          "\n"
                          "# ------------------------------------\n"
                          "# felixsim input\n"
                          "\nhello world"
                          "\n"
                          )
            
            for i in range(len(iv)):
                inpfile.write(iv[i].writeInputLine() + "\n")
            
            inpfile.close()       
            shutil.copyfile(inpFilename, txtFilename)
            
        dlg.Destroy()
        
    def onLoadINP(self, event):
        dlg = wx.FileDialog(self, message="load a .inp file")
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print "Loading: %s" %path
            
            iList = list(range(len(iv)))
            for line in list(open(path, "rb")):
                for i in range(len(iv)):
                    if iv[i].name == line.split(":")[0]:
                        if iv[i].read(line.split()[1]) == True:
                            print iv[i].name + " recognised"
                            #iv[i].read()
                            iList.remove(i)
            for i in iList:
                print iv[i].name + " not recognised"
            
            
        dlg.Destroy()
        
    def onTemp(self, event):
        #print iv[1].inputWidget.SetSelection(7)
        print(None)
        
    def onPrint(self, event):
        print '-'*50
        for i in range(len(iv)):
            print iv[i].writeInputLine()
        print '-'*50
                    


if __name__ == '__main__':    
    app = wx.App(False)
    frame = MainFrame().Show()
    app.MainLoop()