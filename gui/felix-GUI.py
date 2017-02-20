#!/usr/bin/python2.7 

from __future__ import division
from inputvariables import iv
from inputvariables import seperator

import wx
import wx.lib.scrolledpanel
import math
import itertools
import shutil
import os 

### OVERALL FRAME SETUP
class MainFrame(wx.Frame):

    def __init__(self):
      
        wx.Frame.__init__(self, None, title="Felix")
        
        ### TITLE PANEL
        tpanel = wx.Panel(self)

        timage = wx.Image("image.png", wx.BITMAP_TYPE_PNG)
        timage.Rescale(660/6,300/6) #660*300
        timage = wx.BitmapFromImage(timage)
        timage = wx.StaticBitmap(tpanel, -1, timage)
        

        ttext = wx.StaticText(tpanel, wx.ID_ANY, 
            'FELIX\n'
            'GUI\n'
            'CRYSTAL...')
        ttext2 = wx.StaticText(tpanel, wx.ID_ANY, 
            'variables     #########################################\n'
            'felix limited...')
        

        tpanelSizer = wx.BoxSizer(wx.HORIZONTAL)
        tpanelSizer.Add(timage, 0, wx.ALL, 5)
        tpanelSizer.Add(ttext, 0, wx.ALL | wx.CENTRE, 5) 
        tpanelSizer.AddStretchSpacer(1)
        tpanelSizer.Add(ttext2, 2, wx.ALL | wx.CENTRE, 5)         

        tpanel.SetSizer(tpanelSizer)
        tpanelSizer.Fit(tpanel)  

        ### MAIN NOTEBOOK INITIALISE
        notebook = wx.Notebook(self, wx.ID_ANY)

        ## TAB - ABOUT
        taba = wx.Panel(notebook) #tab About (panel)
        text = wx.StaticText(taba, label= (
            "\n\n"
            "########################################\n"
            "first line - This is Felix\n"
            "\n"
            "3rd line - we assume this and have this functionaility"
            "\n\n\n"
            "Best Regards"
            "\n\n"))

        aimage = wx.Image("image.png", wx.BITMAP_TYPE_PNG)
        #aimage.Size((685,300))
        aimage = wx.BitmapFromImage(aimage)
        aimage = wx.StaticBitmap(taba, -1, aimage, size = (660,300))

        tabaSizer = wx.BoxSizer(wx.VERTICAL)
        tabaSizer.Add(text, 0, wx.CENTRE)
        tabaSizer.Add(aimage, 0, wx.CENTRE | wx.ALL,5)               
        taba.SetSizer(tabaSizer)
        tabaSizer.Fit(taba)           

        ## TAB - INPUT INITIALISE
        tabi = wx.Panel(notebook) #tab main (panel)
        
        #tab input panel 1
        mpanel1 = wx.lib.scrolledpanel.ScrolledPanel(tabi, size=(1100,400))
        mpanel1.SetupScrolling(scroll_x=False) #(tab) main panel 1

                #restruture as cycle through iv ready for group seperators
        
        mpanel1Sizer = wx.BoxSizer(wx.VERTICAL)
        perrow = 2
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        onrow = 0
        rmlist = []
        for i in range(len(iv)):
            if isinstance(iv[i],seperator): #or seperator
                for j in range(perrow - onrow):
                    rowSizer.AddStretchSpacer(17)
                mpanel1Sizer.Add(rowSizer, 0, wx.ALL | wx.EXPAND, 2)
                mpanel1Sizer.Add(wx.StaticLine(mpanel1), 0, wx.ALL | wx.EXPAND, 2)
                title = wx.StaticText(mpanel1, label = iv[i].name)
                title.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL, underline=True))
                mpanel1Sizer.Add(title, 0, wx.LEFT, 20)
                mpanel1Sizer.Add(wx.StaticLine(mpanel1), 0, wx.ALL | wx.EXPAND, 2)
                rmlist.append(i)
                rowSizer = wx.BoxSizer(wx.HORIZONTAL)
                onrow = 0
            else:
                print(i)
                if onrow == perrow:
                    mpanel1Sizer.Add(rowSizer, 0, wx.ALL | wx.EXPAND, 2)
                    rowSizer = wx.BoxSizer(wx.HORIZONTAL)
                    onrow = 0
                
                onrow += 1
                rowSizer.Add(wx.StaticText(mpanel1, wx.ID_ANY, iv[i].name), 8, wx.ALL, 5)
                iimage = wx.Image("info3.png", wx.BITMAP_TYPE_PNG)
                iimage.Rescale(25,25) #660*300
                iimage = wx.BitmapFromImage(iimage)
                iimage = wx.StaticBitmap(mpanel1, -1, iimage)
                iimage.Bind(wx.EVT_LEFT_DOWN, iv[i].onInfo)
                rowSizer.Add(iimage, 1, wx.CENTRE)
                #button = wx.Button(mpanel1, label = 'i', style = wx.BU_EXACTFIT)
                #button.Bind(wx.EVT_BUTTON, iv[i].onInfo)
                #rowSizer.Add(button, 1)
                rowSizer.Add(iv[i].inputWidget(mpanel1), 8)
                if onrow < perrow:                
                    rowSizer.AddStretchSpacer(4)    
               
        j = 0       #remove seperators from iv (input variables) 
        for i in rmlist:
            iv.pop(i-j)
            j += 1
        

        mpanel1.SetSizer(mpanel1Sizer)
        mpanel1Sizer.Fit(mpanel1)
       
        #tab input panel2
        mpanel2 = wx.Panel(tabi)
       
        printButton = wx.Button(mpanel2, label = 'print')
        printButton.Bind(wx.EVT_BUTTON, self.onPrint)
        printButton.Bind(wx.EVT_RIGHT_DOWN, self.infoPrint)
        createINPButton = wx.Button(mpanel2, label = 'create .inp')
        createINPButton.Bind(wx.EVT_BUTTON, self.onCreateINP)
        createINPButton.Bind(wx.EVT_RIGHT_DOWN, self.infoPrint)
        loadINPButton = wx.Button(mpanel2, label = 'load from .inp')
        loadINPButton.Bind(wx.EVT_BUTTON, self.onLoadINP)
        loadINPButton.Bind(wx.EVT_RIGHT_DOWN, self.infoPrint)
        
        mpanel2Sizer = wx.BoxSizer(wx.HORIZONTAL)     
        mpanel2Sizer.Add(printButton, 0, wx.ALL, 5)
        mpanel2Sizer.Add(createINPButton, 0, wx.ALL, 5)
        mpanel2Sizer.Add(loadINPButton, 0, wx.ALL, 5)
        
        mpanel2.SetSizer(mpanel2Sizer)
        mpanel2Sizer.Fit(mpanel2)
        
        #TAB - INPUT FINALISE
        tabiSizer = wx.BoxSizer(wx.VERTICAL)
        tabiSizer.Add(mpanel1, 1, wx.CENTRE | wx.EXPAND)
        tabiSizer.Add(wx.StaticLine(tabi), 0, wx.ALL | wx.EXPAND, 2)
        tabiSizer.Add(mpanel2, 0, wx.CENTRE)
        tabi.SetSizer(tabiSizer)
        tabiSizer.Fit(tabi)

        ##TAB - OUTPUT
        tabo = wx.Panel(notebook) #tab About (panel)
        text = wx.StaticText(tabo, label= (
            "From here, you can run\n"
            "the Felix simulation on a\n"
            "directory with the \n"
            "appropiate images and a\n"
            "felix.inp file"))

        
        findButton = wx.Button(tabo, label = 'findFelix')
        findButton.Bind(wx.EVT_BUTTON, self.onFind)
        findButton.Bind(wx.EVT_RIGHT_DOWN, self.infoPrint)
        self.coresInput = wx.TextCtrl(tabo, wx.ID_ANY)
        self.coresInput.Bind(wx.EVT_CHAR, self.checkChar)       
        runButton = wx.Button(tabo, label = 'run')
        runButton.Bind(wx.EVT_BUTTON, self.onRun)
        runButton.Bind(wx.EVT_RIGHT_DOWN, self.infoPrint)

        buttonoSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonoSizer.Add(findButton, 0, wx.LEFT, 40)
        buttonoSizer.Add(wx.StaticText(tabo, wx.ID_ANY, 'cores ='), 0, wx.CENTRE | wx.LEFT, 5)
        buttonoSizer.Add(self.coresInput, 0, wx.LEFT, 2)
        buttonoSizer.Add(runButton, 0, wx.LEFT, 5)
        
        taboSizer = wx.BoxSizer(wx.VERTICAL)
        taboSizer.Add(text, 0, wx.ALL, 20)
        taboSizer.Add(wx.StaticLine(tabo), 0, wx.ALL | wx.EXPAND, 2)
        taboSizer.Add(buttonoSizer)
        taboSizer.Add(wx.StaticLine(tabo), 0, wx.ALL | wx.EXPAND, 2)
        tabo.SetSizer(taboSizer)
        taboSizer.Fit(tabo) 

        ### MAIN NOTEBOOK FINALISE
        notebook.AddPage(taba, "about")
        notebook.AddPage(tabi, "input")
        notebook.AddPage(tabo, "output")
        
        ### OVERALL FRAME SIZERS
        frameSizer = wx.BoxSizer(wx.VERTICAL)
        frameSizer.Add(tpanel, 0, wx.LEFT)            
        frameSizer.Add(notebook, 1, wx.CENTER | wx.EXPAND | wx.ALL, 5)
        
        self.SetSizer(frameSizer)
        self.SetMinSize((1000,300))
        frameSizer.Fit(self)
        


    ### BUTTON FUNCTIONS (inside frame class)
    def onInfo(self, text):
        dlg = wx.MessageDialog(None, text,
            caption = "info", style = wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def checkChar(self, event):
            keycode = event.GetKeyCode()
            if keycode in list(range(48,58))+[46]+[8]+[37]+[39]:
                if(keycode != 46 or '.' not in self.coresInput.GetValue()):
                    event.Skip()
   
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
                                "version:VERSION: Build:BUILD:\n"
                          "# ------------------------------------\n"
                          "\n"
                          "# ------------------------------------\n"
                          "# felixsim input\n"
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
                if (line.split() and line.split()[0] not in ['','#']):
                    for i in range(len(iv)):
                        if iv[i].name == line.split()[0]:
                            value = line.split("=")[1].replace(',',' ').split()
                            if iv[i].read(value) == True:
                                print iv[i].name + " recognised"
                                #iv[i].read()
                                iList.remove(i)
            for i in iList:
                print iv[i].name + " not recognised"
               
        dlg.Destroy()
       
    def onFind(self, event):
        dlg = wx.FileDialog(self, message="select Felix location")
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            print "Loading: %s" %path
            print "current cores selected %s" %self.coresInput.GetValue()
               
        dlg.Destroy()   

    def onRun(self, event):
        dlg = wx.DirDialog(self)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            
            os.chdir(path)
            os.system('mpirun -n 2 ../../src/felixrefine')
            
        dlg.Destroy()
     
    def infoPrint(self, event):
        self.onInfo(
        "This will print to the terminal\n"
        "the current input variables\n"
        "states as shown in this GUI")

    def onPrint(self, event):
        print '-'*50
        for i in range(len(iv)):
            print iv[i].writeInputLine()
        print '-'*50
                    
if __name__ == '__main__':    
    app = wx.App(False)
    frame = MainFrame().Show()
    app.MainLoop()


