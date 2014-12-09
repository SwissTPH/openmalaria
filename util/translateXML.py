#!/bin/env python2
# -*- coding: utf-8 -*-
# Note: ElementTree appears to have bugs in Python 3.3, but works in 2.7
# Most of this code works in both Python 2.x and 3.x, with a little compatibilty option:
from __future__ import print_function

# Tool for translating OpenMalaria XML files from one version to the next

import argparse, os, sys
# We should be able to use lxml or ElementTree. The former may have better formatting.
try:
    from lxml import etree as ET
    using_lxml=True
except ImportError:
    import xml.etree.ElementTree as ET
    using_lxml=False

LATEST_VERSION=33
KNOWN_VERSIONS=[32,33]    # all valid translation targets

class DocumentError(Exception):
    """Class for reporting errors in the input document."""
    def __init__(self,msg):
        Exception(self)
        self.lines = [msg]
    def __lshift__(self,line):
        self.lines.append(line)

# This is just so that comments also get added to the document tree when parsing
class MyTreeBuilder(ET.TreeBuilder):
   def comment(self, data):
       self.start(ET.Comment, {})
       self.data(data)
       self.end(ET.Comment)

# http://stackoverflow.com/a/4590052/314345
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

class Translator:
    def __init__(self,f):
        if using_lxml:
            self.tree = ET.parse(f)
        else:
            self.tree = ET.parse(f, ET.XMLParser(target=MyTreeBuilder()))
        self.root = self.tree.getroot()
        self.schema_ver = int(self.root.get('schemaVersion', 0))
        if self.schema_ver not in KNOWN_VERSIONS:
            err = DocumentError('Unknown document or document version')
            err << 'Found version '+version
            err << 'Known versions: '+str(KNOWN_VERSIONS)
            raise err
        self.namespace = 'http://openmalaria.org/schema/scenario_' + str(self.schema_ver)
        if self.root.tag != self.tag_name():
            err=DocumentError('Unexpected tag name or namespace')
            err << 'Expected <'+self.tag_name()+'>'
            err << 'Found <'+self.root.tag+'>'
            raise err
        self.interval=int(self.root.find('model').find('parameters').get('interval'))
        if self.interval not in (1,5):
            raise DocumentError('Unexpected time step: '+str(self.interval))
        self.steps_per_year = 365 / self.interval
    
    def write_to(self,f):
        indent(self.root)
        if using_lxml:
            self.tree.write(f, encoding='UTF-8', xml_declaration=True, pretty_print=False)
        else:
            ET.register_namespace('om',self.namespace)
            self.tree.write(f, encoding='UTF-8', xml_declaration=True)
    
    def update_to(self,target_version):
        while self.schema_ver < target_version:
            func_name = 'translate_' + str(self.schema_ver) + '_to_' + str(self.schema_ver+1)
            if not func_name in dir(self):      # do we have this function?
                raise NotImplementedError('Translation function '+func_name)
            getattr(self,func_name)()   # call it
            self.schema_ver += 1
            self.set_version_attrib(self.schema_ver)
    
    def tag_name(self):
        return '{' + self.namespace + '}scenario'
    
    def set_version_attrib(self,ver):
        oldVer = str(ver-1)
        newVer = str(ver)
        # This is the tag used by OpenMalaria to tell the version:
        self.root.set('schemaVersion', newVer)
        self.namespace = self.namespace.replace(oldVer, newVer)
        # This is the location of that namespace. It may include the version
        # twice (in the namespace name and the schema file name).
        schemaLocName = '{http://www.w3.org/2001/XMLSchema-instance}schemaLocation'
        self.root.set(schemaLocName, self.root.get(schemaLocName).replace(oldVer, newVer))
        if using_lxml:
            #lxml does not appear to be able to update the nsmap of an existing element
            old_root = self.root
            self.root = ET.Element(self.tag_name(), attrib=old_root.attrib, nsmap={'om':self.namespace})
            for child in old_root:
                self.root.append(child)
            self.tree = ET.ElementTree(self.root)
        else:
            self.root.tag = self.tag_name()
    
    def format_time(self, steps):
        # TODO: could use dates
        if steps % self.steps_per_year == 0:
            return str(steps / self.steps_per_year)+'y'
        return str(steps * self.interval)+'d'
    
    def get_screen_components(self, elt):
        components = []
        for deploy in elt.findall('deploy'):
            for depComp in deploy.findall('component'):
                if float(depComp.get('minAge',0.0)) != 0.0 or \
                    float(depComp.get('p',1.0)) != 1.0 or \
                    depComp.get('maxAge') != None:
                        # could do this via decisionTree instead
                        raise DocumentError('Lazily refusing to update <screen> with restricted deployments')
                components.append(depComp.get('id'))
        return components
    
    def update_survey_times(self, survey_times):
        """Do one pass through the list of times searching for patterns.
        Return true if any changes were made and should be re-run.
        Note that this search engine is sub-optimal in that if surveys are
        regularly taken at say 5 times every year over ten years it will not
        find the best set of sequences (hence why at least 5 matches are required)."""
        for i,start in enumerate(survey_times[0]):
            for y in survey_times[0][i+1:]:
                stride = y - start
                # start, start+stride are by definition in survey_times[0]
                # Consider a sequence if next two terms are present:
                if start+2*stride in survey_times[0] and \
                    start+3*stride in survey_times[0] and start+4*stride in survey_times[0]:
                    # match; calculate end and remove all matching elements:
                    j = 0
                    while True:
                        n = start + j*stride
                        if n in survey_times[0]:
                            survey_times[0].remove(n)
                            end = n+1
                            j += 1
                            continue
                        break
                    survey_times[1].append((start,stride,end))
                    # Now we need to restart our iteration
                    return True
        return False
    
    def translate_32_to_33(self):
        """Version 33 has many changes but is mostly backwards compatible.
        
        Translated:
        Intervention components: MDA → treatSimple
        Intervention component: syntax of 'screen' changed (requires manual fix)
        GARKI_DENSITY_BIAS option removed
        
        Not translated:
        Intervention components: MDA1D → treatPKPD or decisionTree
        Intervention components: MDA → decisionTree
        Syntax of decision trees in EventScheduler, MDA1D changed massively
        Time units and dates can be used instead of the old step numbers
        DecisionTree5Day is available as an optional replacement for ImmediateOutcomes
        Standard diagnostic may be specified in a library
        """
        for ev in self.root.iter('EventScheduler'):
            # i.e. if we have any elements by this name
            raise DocumentError('Cowardly refusing to update an XML using <EventScheduler>')
        
        # Compress survey time list if possible:
        surveys = self.root.find('monitoring').find('surveys')
        survey_times = ([],[]) # first list is individual times, second is for (start,stride,end) triples
        for survey in surveys.findall('surveyTime'):
            survey_times[0].append(int(survey.text))
            surveys.remove(survey)
        while self.update_survey_times(survey_times):
            pass
        # Now merge into a single list and sort:
        for time in survey_times[0]:
            survey_times[1].append((time,None,None))
        survey_times = sorted(survey_times[1], key = lambda tup: tup[0])
        for (time,stride,end) in survey_times:
            elt = ET.Element('surveyTime')
            elt.text = self.format_time(time)
            if stride != None:
                elt.set('repeatStep',self.format_time(stride))
                elt.set('repeatEnd',self.format_time(end))
            surveys.append(elt)
        
        # Translate diagnostics
        diagnostics = ET.Element('diagnostics')
        i = 0
        while i < len(self.root) and self.root[i] != 'model':
            i += 1
        self.root.insert(i-1, diagnostics)
        
        stdDiagName='standard'
        stdDiagUnits='Other'
        # Remove an obseleted option:
        modelOpt = self.root.find('model').find('ModelOptions')
        for opt in modelOpt:
            if opt.get('name') == 'GARKI_DENSITY_BIAS' and opt.get('value').lower() == 'true':
                stdDiagName='garki'
                stdDiagUnits='Garki'
                modelOpt.remove(opt)
                break
        
        # Move standard diagnostic description to new location:
        stdDiag = ET.Element('diagnostic', {'name':stdDiagName, 'units':stdDiagUnits})
        stdDiag.append(ET.Element('deterministic', {'minDensity':surveys.get('detectionLimit')}))
        diagnostics.append(stdDiag)
        surveys.set('diagnostic',stdDiagName)
        del surveys.attrib['detectionLimit']
        
        # Also set a neonatal mortality diagnostic (this is essentially a bug-fix)
        # Reason: the neonatal model was calibrated against this diagnostic
        neonatalDiag = ET.Element('diagnostic', {'name':'neonatal', 'units':'Other'})
        neonatalDiag.append(ET.Element('deterministic', {'minDensity':'40'}))
        diagnostics.append(neonatalDiag)
        self.root.find('model').find('clinical').insert(0,ET.Element('NeonatalMortality', {'diagnostic':'neonatal'}))
        
        intervs = self.root.find('interventions')
        for human in intervs.findall('human'):  # should only be one
            for comp in human.findall('component'):
                for child in comp:
                    if child.tag == 'MDA1D':
                        raise DocumentError('Cowardly refusing to update an XML using <MDA1D>')
                    elif child.tag == 'MDA':
                        for eff in child.findall('effects'):    # only one
                            opts = eff.findall('option')
                            if len(opts) > 1 or len(list(eff.iter('deploy'))) > 0:
                                # Could be translated to decisionTree without too much difficulty
                                raise DocumentError('Lazily refusing to update an <MDA> element not equivalent to some <treatSimple> element')
                            for opt in opts:
                                assert float(opt.get('pSelection')) == 1.0
                                name=opt.get('name')    #optional attribute
                                tsLiver = 0
                                tsBlood = 0
                                for clearInf in opt.findall('clearInfections'):
                                    ts = int(clearInf.get('timesteps'))
                                    stage = clearInf.get('stage')
                                    if stage == 'liver' or stage == 'both':
                                        tsLiver = ts
                                    if stage == 'blood' or stage == 'both':
                                        tsBlood = ts
                        comp.remove(child)
                        # don't use format_time here: -1 is special value
                        treatSimple = ET.Element('treatSimple',
                                                 {'durationLiver':str(tsLiver)+'t',
                                                  'durationBlood':str(tsBlood)+'t'})
                        comp.append(treatSimple)
                    elif child.tag =='screen':
                        diagElt = child.find('diagnostic')
                        posComponents = self.get_screen_components(child.find('positive'))
                        negComponents = self.get_screen_components(child.find('negative'))
                        
                        child.remove(diagElt)
                        child.remove(child.find('negative'))
                        child.remove(child.find('positive'))
                        
                        diagnostics.append(diagElt)
                        diagName = 'diagnostic' + str(len(diagnostics))
                        diagElt.set('name',diagName)
                        diagElt.set('units','?REPLACEME?')
                        print('Warning: units of diagnostic',diagName,'must be set manually')
                        print('Warning: replace ?REPLACEME? with Malariatherapy, Garki, or Other')
                        
                        child.set('diagnostic', diagName)
                        for comp_id in posComponents:
                            child.append(ET.Element('positive',{'id':comp_id}))
                        for comp_id in negComponents:
                            child.append(ET.Element('negative',{'id':comp_id}))

def die(*args):
    print(*args, file=sys.stderr)
    sys.exit(1)

def readArgs():
    parser = argparse.ArgumentParser(description="This tool translates one or more "+
            "OpenMalaria scenario files (XML) from one version to the next.")
    parser.add_argument('files', metavar='FILE', nargs='*',
                        help="A file to translate. Multiple files may be specified.")   
    parser.add_argument('-t', '--target-version', type=int, action="store",
                        metavar='VER', default=LATEST_VERSION, choices=KNOWN_VERSIONS,
                        help="Target version for translation. Default: "+str(LATEST_VERSION))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d','--dest', metavar='DIR', action='store', default=None,
                        help='Destination directory for translated scenarios')
    group.add_argument('-i', '--in-place', action='store_true', default=False,
                        help='If set, files will be updated in-place (incompatible with -d option')
    args = parser.parse_args()
    if args.dest and not os.path.isdir(args.dest):
        die('Not a directory:', args.dest)
    return args

def main():
    args = readArgs()
    # args.target_version : target version number
    # args.files : list of files to translate
    # args.dest : destination directory or None
    # args.in_place : overwrite input files (True) or not (False)
    for f_path in args.files:
        try:
            if not os.path.isfile(f_path):
                die('Not a file:', f_path)
            with open(f_path, 'r') as f:
                translator = Translator(f)
            print('Translating',f_path,'from',translator.schema_ver,'to',args.target_version)
            translator.update_to(args.target_version)
            if args.in_place:
                with open(f_path, 'w') as f:
                    translator.write_to(f)
            elif args.dest:
                path = os.path.join(args.dest, os.path.basename(f_path))
                with open(path, 'w') as f:
                    translator.write_to(f)
        except DocumentError,e:
            for line in e.lines:
                print('Error:',line)

if __name__ == "__main__":
    main()
