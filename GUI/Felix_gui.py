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

class Notebook(wx.Notebook):

  def __init__(self, parent):
    wx.Notebook.__init__(self, parent, wx.ID_ANY)

    # Add pages to notebook for different tabs
    self.page1 = GuiPages.FlagPanel(self)
    self.page2 = GuiPages.RadiusPanel(self)
    self.page3 = GuiPages.BeamPanel(self)
    self.page4 = GuiPages.crystalPanel(self)
    self.page5 = GuiPages.microscopePanel(self)
    self.page6 = GuiPages.imagePanel(self)

    self.AddPage(self.page1, "Flags")
    self.AddPage(self.page2, "Radius")
    self.AddPage(self.page3, "Beam")
    self.AddPage(self.page4, "Crystal")
    self.AddPage(self.page5, "Microscope")
    self.AddPage(self.page6, "Image")


class MainFrame(wx.Frame):

  def __init__(self):

    wx.Frame.__init__(self, None, title="Felix")

    # Create a panel with a notebook on it
    panel = wx.Panel(self)
    self.notebook = Notebook(panel)
    self.option = GuiPages.optionPanel(panel, self)
    self.viewer = GuiPages.ViewerPanel(panel)

    sizerh = wx.BoxSizer(wx.HORIZONTAL)
    sizerv = wx.BoxSizer(wx.VERTICAL)

    # add the widgets to the sizers
    sizerv.Add(self.notebook, 0, wx.ALL, 5)
    sizerv.Add(self.option, 0, wx.ALL | wx.CENTER, 5)
    sizerh.Add(sizerv, 0)
    sizerh.Add(self.viewer, 0, wx.ALL | wx.EXPAND, 5)


    panel.SetSizer(sizerh)
    sizerh.Fit(self)

    self.Show()


if __name__ == '__main__':
  app = wx.App()
  frame = MainFrame().Show()
  app.MainLoop()
