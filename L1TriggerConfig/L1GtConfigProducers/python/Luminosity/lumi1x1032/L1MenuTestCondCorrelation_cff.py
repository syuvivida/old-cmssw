import FWCore.ParameterSet.Config as cms

# cff file for testing correlation conditions L1MenuTestCondCorrelation
# to be added after L1GtConfig.cff

# menu definition
from L1TriggerConfig.L1GtConfigProducers.l1GtTriggerMenuXml_cfi import *
l1GtTriggerMenuXml.TriggerMenuLuminosity = 'lumi1x1032'
l1GtTriggerMenuXml.DefXmlFile = 'L1MenuTestCondCorrelation.xml'
l1GtTriggerMenuXml.VmeXmlFile = ''

# prescale factors, trigger masks, trigger veto masks
# default: no prescale, no bit masked, no bit vetoed 

