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

  def __init__(self, parent, WikiObj):
    wx.Notebook.__init__(self, parent, wx.ID_ANY)

    # Add pages to notebook for different tabs
    self.page1 = GuiPages.FlagPanel(self, WikiObj)
    self.page2 = GuiPages.RadiusPanel(self, WikiObj)
    self.page3 = GuiPages.BeamPanel(self, WikiObj)
    self.page4 = GuiPages.crystalPanel(self, WikiObj)
    self.page5 = GuiPages.microscopePanel(self, WikiObj)
    self.page6 = GuiPages.imagePanel(self, WikiObj)
    #self.wiki = GuiPages.WikiPanel(self)


    self.AddPage(self.page1, "Flags")
    self.AddPage(self.page2, "Radius")
    self.AddPage(self.page3, "Beam")
    self.AddPage(self.page4, "Crystal")
    self.AddPage(self.page5, "Microscope")
    self.AddPage(self.page6, "Image")
    #self.AddPage(self.wiki, "wiki")

class WikiPreviewNotebook(wx.Notebook):
  def __init__(self, parent):
    wx.Notebook.__init__(self, parent, wx.ID_ANY)

    self.wikitext = GuiPages.WikiPanel(self)
    self.viewer = GuiPages.ViewerPanel(self)

    self.AddPage(self.wikitext, "Wiki")
    self.AddPage(self.viewer, "Image Preview")

class MainFrame(wx.Frame):

  def __init__(self):

    wx.Frame.__init__(self, None, title="Felix")
    # Create a panel with a notebook on it
    panel = wx.Panel(self)
    self.option = GuiPages.optionPanel(panel, self)

    self.wiki = WikiPreviewNotebook(panel)
    self.notebook = Notebook(panel, self.wiki)

    sizerh = wx.BoxSizer(wx.HORIZONTAL)
    sizerv = wx.BoxSizer(wx.VERTICAL)

    # add the widgets to the sizers
    sizerv.Add(self.notebook, 0, wx.ALL, 5)
    sizerv.Add(self.option, 0, wx.ALL | wx.CENTER, 5)
    sizerh.Add(sizerv, 2)
    sizerh.Add(self.wiki, 1, wx.ALL | wx.EXPAND, 5)


    panel.SetSizer(sizerh)
    sizerh.Fit(self)
    self.Center()

    self.Layout()


if __name__ == '__main__':
  app = wx.App()
  frame = MainFrame().Show()
  app.MainLoop()
