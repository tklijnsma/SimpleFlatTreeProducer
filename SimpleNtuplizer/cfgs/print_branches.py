#from random import shuffle
#from copy import copy

import ROOT

#sample_dir = '/pnfs/psi.ch/cms/trivcat/store/user/tklijnsm/DoublePhoton_FlatPt-5To300/crab_TK_Ntuplizer_A2/160512_144458/0000/'
#to_open_base = 'dcap://t3se01.psi.ch:22125/'

sample_dir   = ''
to_open_base = ''

physical_path = lambda output_file: to_open_base + sample_dir + output_file



output_file = 'DoubleElectron_AODSIM_example.root'




print 'Trying to read ' + output_file

root_fp = ROOT.TFile.Open( physical_path(output_file) )

EventTree = root_fp.Get( 'Events' )

for branch in EventTree.GetListOfBranches():
    print branch.GetName()
