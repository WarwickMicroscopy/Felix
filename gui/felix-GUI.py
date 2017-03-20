#!/usr/bin/python2.7 


from __future__ import division

import time
t0 = time.time()

from inputvariables import iv
from inputvariables import seperator

print(time.time() - t0)
t0 = time.time()

import wx
import wx.lib.scrolledpanel
import math
import itertools
import shutil
import os

print(time.time() - t0) 


### OVERALL MAIN FRAME SETUP
class MainFrame(wx.Frame):

    def __init__(self):
      
        t0 = time.time()

        wx.Frame.__init__(self, None, title="Felix")
        

        ### TITLE PANEL
        tpanel = wx.Panel(self)

        timage = wx.Image("crystal.png", wx.BITMAP_TYPE_PNG)
        timage.Rescale(660/6,300/6) #660*300
        timage = wx.BitmapFromImage(timage)
        timage = wx.StaticBitmap(tpanel, -1, timage)
        

        ttext = wx.StaticText(tpanel, wx.ID_ANY, 
            'FELIX\n'
            'GUI\n'
            '+CRYSTALS')
        ttext2 = wx.StaticText(tpanel, wx.ID_ANY,
            '\nFelix: Bloch wave method diffraction pattern simulation software\n' 
            'Richard Beanland, Keith Evans, Rudolf A Roemer and Alex Hubert')
        

        tpanelSizer = wx.BoxSizer(wx.HORIZONTAL)
        tpanelSizer.Add(timage, 0, wx.ALL, 5)
        tpanelSizer.Add(ttext, 0, wx.ALL | wx.CENTRE, 5) 
        tpanelSizer.AddStretchSpacer(1)
        tpanelSizer.Add(ttext2, 2, wx.ALL | wx.CENTRE, 5)         

        tpanel.SetSizer(tpanelSizer)
        tpanelSizer.Fit(tpanel)

        print(time.time() - t0) 
        t0 = time.time()  

        ### MAIN NOTEBOOK INITIALISE
        notebook = wx.Notebook(self, wx.ID_ANY)

        ## TAB - ABOUT
        taba = wx.Panel(notebook) #tab About (panel)
        atext = wx.StaticText(taba, label= (
            "\n"          
            "(C) 2013/14, all rights reserved\n\n"
            "Free software under the GNU General Public License\n\n"
            "This GUI has seperate tabs, with capabilities to modify\n"
            "the input .inp file, run the simulation and display\n"
            "the images\n"
            "\n\n"
            "Right click on widgets or click on the info icon for\n"
            "info dialogs explaining their uses in more detail"
            "\n\n"))

        aimage = wx.Image("crystal.png", wx.BITMAP_TYPE_PNG)
        aimage = wx.BitmapFromImage(aimage)
        aimage = wx.StaticBitmap(taba, -1, aimage, size = (660,300))
        aimage2 = wx.Image("crystal2.jpg", wx.BITMAP_TYPE_JPEG)
        aimage2.Rescale(400,400)
        aimage2 = wx.BitmapFromImage(aimage2)
        aimage2 = wx.StaticBitmap(taba, -1, aimage2, size = (250,250))

        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        rowSizer.Add(aimage2, 0, wx.RIGHT | wx.CENTRE, 10)
        rowSizer.Add(atext, 0, wx.CENTRE)
        tabaSizer = wx.BoxSizer(wx.VERTICAL)
        tabaSizer.Add(rowSizer, 0, wx.CENTRE | wx.TOP, 5)
        tabaSizer.Add(aimage, 0, wx.CENTRE | wx.ALL,5)               
        taba.SetSizer(tabaSizer)
        tabaSizer.Fit(taba)

        print(time.time() - t0) 
        t0 = time.time()           

        ## TAB - INPUT INITIALISE
        tabi = wx.Panel(notebook) #tab main (panel)
        
        #tab input panel 1
        ipanel = wx.lib.scrolledpanel.ScrolledPanel(tabi, size=(1100,400))
        ipanel.SetupScrolling() #(tab) main panel 1

                #restruture as cycle through iv ready for group seperators
        
        ipanelSizer = wx.BoxSizer(wx.VERTICAL)
        perrow = 2
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        onrow = 0
        rmlist = []
        for i in range(len(iv)):
            if isinstance(iv[i],seperator): #or seperator
                for j in range(perrow - onrow):
                    rowSizer.AddStretchSpacer(17)
                ipanelSizer.Add(rowSizer, 0, wx.ALL | wx.EXPAND, 2)
                ipanelSizer.Add(wx.StaticLine(ipanel), 0, wx.ALL | wx.EXPAND, 2)
                title = wx.StaticText(ipanel, label = iv[i].name)
                title.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL, underline=True))
                ipanelSizer.Add(title, 0, wx.LEFT, 20)
                ipanelSizer.Add(wx.StaticLine(ipanel), 0, wx.ALL | wx.EXPAND, 2)
                rmlist.append(i)
                rowSizer = wx.BoxSizer(wx.HORIZONTAL)
                onrow = 0
            else:
                if onrow == perrow:
                    ipanelSizer.Add(rowSizer, 0, wx.ALL | wx.EXPAND, 2)
                    rowSizer = wx.BoxSizer(wx.HORIZONTAL)
                    onrow = 0
                
                onrow += 1
                rowSizer.Add(wx.StaticText(ipanel, wx.ID_ANY, iv[i].name), 8, wx.ALL, 5)
                iimage = wx.Image("info.png", wx.BITMAP_TYPE_PNG)
                #iimage.Rescale(25,25) #660*300
                iimage = wx.BitmapFromImage(iimage)
                iimage = wx.StaticBitmap(ipanel, -1, iimage)
                iimage.Bind(wx.EVT_LEFT_DOWN, lambda evt, text=iv[i].infoText: self.onInfo(evt,text))  
                rowSizer.Add(iimage, 1, wx.CENTRE)                
                rowSizer.Add(iv[i].inputWidget(ipanel), 8)
                #iv[i].widget.Bind(wx.EVT_RIGHT_DOWN, lambda evt, text=iv[i].infoText: self.onInfo(evt,text)) 
                if onrow < perrow:                
                    rowSizer.AddStretchSpacer(4)    
               
        j = 0       #remove seperators from iv (input variables) 
        for i in rmlist:
            iv.pop(i-j)
            j += 1
        
        print(time.time() - t0) 
        t0 = time.time()

        ipanel.SetSizer(ipanelSizer)
        ipanelSizer.Fit(ipanel)
       
        #tab input panel2
        ipanel2 = wx.Panel(tabi)
       
        printButton = wx.Button(ipanel2, label = 'print')
        printButton.Bind(wx.EVT_BUTTON, self.onPrint)
        printButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Print current variables state to terminal': self.onInfo(evt,text))
        createINPButton = wx.Button(ipanel2, label = 'create .inp')
        createINPButton.Bind(wx.EVT_BUTTON, self.onCreateINP)
        createINPButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Create .inp file with current variable states': self.onInfo(evt,text))
        loadINPButton = wx.Button(ipanel2, label = 'load from .inp')
        loadINPButton.Bind(wx.EVT_BUTTON, self.onLoadINP)
        loadINPButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Load variable states from .inp file': self.onInfo(evt,text))
        
        ipanel2Sizer = wx.BoxSizer(wx.HORIZONTAL)     
        ipanel2Sizer.Add(printButton, 0, wx.ALL, 5)
        ipanel2Sizer.Add(createINPButton, 0, wx.ALL, 5)
        ipanel2Sizer.Add(loadINPButton, 0, wx.ALL, 5)
        
        ipanel2.SetSizer(ipanel2Sizer)
        ipanel2Sizer.Fit(ipanel2)
        
        #TAB - INPUT FINALISE
        tabiSizer = wx.BoxSizer(wx.VERTICAL)
        tabiSizer.Add(ipanel, 1, wx.CENTRE | wx.EXPAND)
        tabiSizer.Add(wx.StaticLine(tabi), 0, wx.ALL | wx.EXPAND, 2)
        tabiSizer.Add(ipanel2, 0, wx.CENTRE)
        tabi.SetSizer(tabiSizer)
        tabiSizer.Fit(tabi)

        print(time.time() - t0) 
        t0 = time.time()

        ##TAB - OUTPUT
        tabo = wx.Panel(notebook) #tab About (panel)
        otext = wx.StaticText(tabo, label= (
            "From here, you can run the Felix simulation on a\n"
            "directory with the appropiate images and a\n"
            "felix.inp file"))
        
        findButton = wx.Button(tabo, label = 'findFelix')
        findButton.Bind(wx.EVT_BUTTON, self.onFind)
        findButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Find and select the felix run file': self.onInfo(evt,text))
        self.findText = wx.TextCtrl(tabo, wx.ID_ANY, size = (250,-1))
        self.findText.SetEditable(False)
        self.coresInput = wx.TextCtrl(tabo, wx.ID_ANY, value = '1', size = (40,-1))
        self.coresInput.Bind(wx.EVT_CHAR, self.checkChar)
        self.coresInput.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Type number of cores for mpirun': self.onInfo(evt,text))       
        runButton = wx.Button(tabo, label = 'run')
        runButton.Bind(wx.EVT_BUTTON, self.onRun)
        runButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Run Felix simulation on chosen directory assuming felix.inp file exists': self.onInfo(evt,text))
        viewerButton = wx.Button(tabo, label = 'viewer')
        viewerButton.Bind(wx.EVT_BUTTON, self.onViewer)
        viewerButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='View images in a folder': self.onInfo(evt,text))
        convertButton = wx.Button(tabo, label = 'convert')
        convertButton.Bind(wx.EVT_BUTTON, self.onConvert)
        convertButton.Bind(wx.EVT_RIGHT_DOWN, lambda evt, 
            text='Convert a .bin into a .gif': self.onInfo(evt,text))        

        buttonoSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonoSizer.Add(findButton, 0, wx.LEFT)
        buttonoSizer.Add(self.findText, 0, wx.LEFT, 2)
        buttonoSizer.Add(wx.StaticText(tabo, wx.ID_ANY, 'cores ='), 0, wx.CENTRE | wx.LEFT, 10)
        buttonoSizer.Add(self.coresInput, 0, wx.LEFT, 2)
        buttonoSizer2 = wx.BoxSizer(wx.HORIZONTAL)
        buttonoSizer2.Add(viewerButton, 0, wx.LEFT, 10)
        buttonoSizer2.Add(convertButton, 0, wx.LEFT, 10)
        
        taboSizer = wx.BoxSizer(wx.VERTICAL)
        taboSizer.Add(otext, 0, wx.CENTRE | wx.ALL, 20)
        taboSizer.Add(wx.StaticLine(tabo), 0, wx.ALL | wx.EXPAND, 2)
        taboSizer.Add(buttonoSizer, 0, wx.CENTRE)
        taboSizer.Add(wx.StaticLine(tabo), 0, wx.ALL | wx.EXPAND, 2)
        taboSizer.Add(runButton, 0, wx.CENTRE)
        taboSizer.Add(wx.StaticLine(tabo), 0, wx.ALL | wx.EXPAND, 2)
        taboSizer.Add(buttonoSizer2, 0, wx.CENTRE)
        taboSizer.Add(wx.StaticLine(tabo), 0, wx.ALL | wx.EXPAND, 2)
        tabo.SetSizer(taboSizer)
        taboSizer.Fit(tabo)

        print(time.time() - t0) 
        t0 = time.time() 

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
        
        print(time.time() - t0) 
        t0 = time.time() 

    ### BUTTON FUNCTIONS (inside frame class)
    def onInfo(self, event, text):
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
            self.findText.SetValue(path)
            print "Loading: %s" %path
            print "current cores selected %s" %self.coresInput.GetValue()
               
        dlg.Destroy()   

    def onRun(self, event):
        dlg = wx.DirDialog(self)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            
            os.chdir(path)
            os.system('mpirun -n ' + self.coresInput.GetValue() + ' ' + self.findText.GetValue())
            
        dlg.Destroy()
    
    def onViewer(self, event):
        dlg = wx.DirDialog(self)
        
        if dlg.ShowModal() == wx.ID_OK:
            viewerFrame = aviewerFrame(str(dlg.GetPath()))
            viewerFrame.Show()
            
        dlg.Destroy()        

    def onPrint(self, event):
        print '-'*50
        for i in range(len(iv)):
            print iv[i].writeInputLine()
        print '-'*50

    def onConvert(self, event):
        dlg = wx.FileDialog(self, message="select .bin to convert")
        
        if dlg.ShowModal() == wx.ID_OK:
            f = dlg.GetPath()
            print('converting following to .gif:' + f)
            os.system('convert -size 539x539 -depth 64 -define quantum:format=floating-point'
                ' -define quantum:scale=65535.0 -endian lsb GRAY:'+ f + ' ' + f.split('.')[0] + '.gif')
               
        dlg.Destroy() 

### VISUALISER FRAME
class aviewerFrame(wx.Frame):

    def __init__(self, path):

        wx.Frame.__init__(self, None, title="Visulise")
        panel = wx.lib.scrolledpanel.ScrolledPanel(self, size=(500,500))
        panel.SetupScrolling(scroll_x=False)
        print(path)
        print(os.listdir(path))

        sizer = wx.BoxSizer(wx.VERTICAL)
        
        for f in os.listdir(path):
            if os.path.splitext(f)[1].lower() in ('.jpg', '.png'):
                print(f)
                image = wx.Image(path + '/' + f, wx.BITMAP_TYPE_ANY)
                x,y = image.GetSize()
                print(x,y)
                c = 300/max(x,y)
                image.Rescale(c*x,c*y) 
                image = wx.BitmapFromImage(image)
                image = wx.StaticBitmap(panel, -1, image)
                sizer.Add(image, 0, wx.CENTRE | wx.ALL, 5)
                
        panel.SetSizer(sizer)

if __name__ == '__main__':    
    app = wx.App(False)
    frame = MainFrame().Show()
    app.MainLoop()


