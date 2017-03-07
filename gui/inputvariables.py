import itertools
import wx

### INPUT VARIABLES & WIDGET CLASSES
def combof(pool, maxn):
    combolist = []
    indexlist = []    
    for r in range(1, maxn + 1):
        n = len(pool)
        #if r > n:
        indices = range(0, r)
        tempstr = pool[0]
        temparray = [0]
        for i in indices[1:]:
            tempstr += ', ' + pool[i]
            temparray.append(i)

        combolist.append(tempstr)
        indexlist.append(temparray)

        while True:
            for i in reversed(range(r)):
                if indices[i] != i + n - r:
                    break
            else:
                break
            indices[i] += 1
            for j in range(i+1, r):
                indices[j] = indices[j-1] + 1
            tempstr = pool[indices[0]]
            temparray = [indices[0]]
            for i in indices[1:]:
                tempstr += ', ' + pool[i]
                temparray.append(i)
            
            combolist.append(tempstr)
            indexlist.append(temparray)

    return combolist, indexlist

class input_var(): #input variable
    
    def __init__(self, name, infoText=''):

        self.name = name
        if infoText == '':
            self.infoText = name
        else: 
            self.infoText = infoText

        #check name isinstance(name, str)
       
class combo_var(input_var): #combo variable - select option from list

    def __init__(self, name, choices, defaultRef=0, maxchoice = -1, infoText=''):
        if maxchoice == -1:
            maxchoice = 1

        #all(isinstance(elem, str) for elem in choices)
        input_var.__init__(self, name, infoText)
        self.choices, self.refList = combof(choices, maxchoice)
        self.defaultRef = defaultRef
        
    def inputWidget(self, panel): 
            
        self.widget = wx.Choice(panel, wx.ID_ANY, size=(100, -1),
                                     choices=self.choices, name=self.name)
        self.widget.SetStringSelection(self.choices[self.defaultRef])
        #self.widget.Bind(wx.EVT_RIGHT_DOWN, self.onInfo)
              
        return self.widget
        
    def writeInputLine(self):
        choiceRef = str(self.refList[self.widget.GetSelection()]).strip('[]')
        i = 30 - len(self.name)
        j = 15 - len(choiceRef)
        
        return self.name + ' '*i + '= ' + choiceRef

    def read(self, choiceRef):
        try: 
            choiceRef = [int(i) for i in choiceRef]
        except:
            return False
        else:
            if choiceRef in self.refList:
                self.widget.SetSelection(self.refList.index(choiceRef))
                return True
               
            else:
                return False


class value_var(input_var):
    
    def __init__(self, name, infoText=''):
        input_var.__init__(self, name, infoText)        

                        
    def inputWidget(self, panel):
        self.widget = wx.TextCtrl(panel, wx.ID_ANY)

        def checkChar(event):
            keycode = event.GetKeyCode()
            if keycode in list(range(48,58))+[46]+[8]+[37]+[39]:
                if(keycode != 46 or '.' not in self.widget.GetValue()):
                    event.Skip()
                    self.chosen = self.widget.GetValue()
                
                
        self.widget.Bind(wx.EVT_CHAR, checkChar)

        return self.widget
        #return wx.TextCtrl(panel, wx.ID_ANY)
                    
    def writeInputLine(self):
        i = 30 - len(self.name)
        if self.widget.GetValue() == '':
            value = 0
        else:
            value = self.widget.GetValue()
        return self.name + ' '*i + '= ' + str(value)
        
    def read(self, value):
        try:
            value = str(value[0])
            #print(value)
            self.widget.SetValue(value)
            return True
        except ValueError:
            return False


class seperator():

    def __init__(self, name):
        self.name = name


        
###INPUT VARIABLES      
iv =[
    seperator('Flag Variables'),
    combo_var('IWriteFLAG', ['0','5','11'], infoText = 'For detailed terminal info...'),    
    combo_var('IImageFLAG', ['Diffractions','Reflections',
                    'Amplitude+Phase'], 2, maxchoice = 3),
    combo_var('IOutputFLAG', ['Nothing','Ug Matrix','Eigenspectra',
                    'Wavefunctions'], 0),

    seperator('Beam Variables'),
    combo_var('IBinorTextFLAG', ['Binary','Text'], 1),
    combo_var('IScatteringFactorMethod', ['Kirklands','Peng',
                    'Doyle+Turner'], 0),
    combo_var('IZolzFLAG', ['No','Yes'], 0),
    combo_var('ICentralBeamFlag', ['Not included','included'], 0),
    combo_var('IMaskFLAG', ['Circular','Square'], 0),
    combo_var('IAbsorbFLAG', ['No','Proprtional','Einstein Perturbative',
                    'Einstein Exact'], 0),
    value_var('IMinReflectionPool'),
    combo_var('IAnisotropicDebyeWallerFLAG', ['No','Yes'], 0),
    combo_var('IPseudoCubicFLAG', ['Orthorhombic','PseudoCubic'], 0),
    combo_var('IXDirectionFLAG', ['Ignore','Use'], 0),
    value_var('IPixelCount'),
    value_var('IMinStrongBeams'),
    value_var('IMinWeakBeams'),
    value_var('RBSBmax'),

    seperator('Directional Vectors'),
    value_var('RBSPmax'),
    value_var('RDebyeWallerConstant'),
    value_var('RAbsorptionPer'),
    value_var('RConvergenceAngle'),
    value_var('IIncidentBeamDirectionX', infoText = 'IIncidentBeamDirectionX coordinate'),
    value_var('IIncidentBeamDirectionY'),
    value_var('IIncidentBeamDirectionZ'),
    value_var('IXDirectionX'),

    seperator('Miscellaneous Variables'),
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
    seperator('')
    ]

