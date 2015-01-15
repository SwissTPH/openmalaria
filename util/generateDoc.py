#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of OpenMalaria.

Copyright (C) 2005-2014 Swiss Tropical Institute

OpenMalaria is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
"""

import argparse, sys, os.path, re
#try:
    #import lxml.etree as ET
#except ImportError:
import xml.etree.ElementTree as ET

linkbase='PLACEHOLDER' # note: if not specifying an output file, links will not work

all_links = set()
def make_link_from_name(name):
    a = 'elt-' + name
    global all_links
    i = 1
    link = a
    while True:
        if not link in all_links:
            all_links.add(link)
            return link
        i += 1
        link = a + '-' + str(i)

class DocWriter:
    """Writes formatted text to a file"""
    def __init__(self,f_out):
        self.f = f_out
    
    #"""wiki formatting functions: Google Code wiki"""
    #def heading(self, level, *title):
        #"""level: 1 is biggest heading, 2 next, 5 smallest"""
        #tag = "=" * level
        #args=[tag]
        #args.extend(title)
        #args.append(tag)
        #self.line(*args)
    #def bold(self, text):
        #return '*' + text + '*'
    #def link(self, target, text, internal=True):
        #"""return text for a link, to be injected into a line"""
        ## Format compatible with Google wiki:
        #if internal:
            #target = linkbase + '#' + target
        #else:
            #target = target
        #return '['+target+' '+text+']'
    #def bulleted(self, *args):
        #self.line('  *', *args)
    #def startcode(self):
        #self.line('{{{')
    #def endcode(self):
        #self.line('}}}')
    
    """wiki formatting functions: GitHub Markdown"""
    def heading(self, level, text):
        """level: 1 is biggest heading, 2 next, 5 smallest"""
        if level <= 2:
            self.line(text)
            char = '=' if level == 1 else '-'
            self.line(char * max(len(text), 3))
        else:
            self.line('#' * level, text)
        #self.line('<h'+str(level)+' id="'+text.replace(' ','_')+'">' + text + '</h'+str(level)+'>')
    def bold(self, text):
        return '**' + text + '**'
    def link(self, target, text, internal=True):
        """return text for a link, to be injected into a line"""
        # intended for internal links in HTML documents:
        if internal:
            # compatible with GitHub links:
            target = '#' + target
        else:
            target = target
        return '['+text+']('+target+')'
    def bulleted(self, *args):
        self.line('*  ', *args)
    def startcode(self, lang=''):
        self.line('```' + lang)
    def endcode(self):
        self.line('```')
    
    def head(self, schema_file, ver):
        # write header
        self.heading(1, 'Generated schema '+ver+' documentation')
        self.line('This page is automatically generated from the following schema file: `'+schema_file+'`.')
        self.line('I recommend against editing it because edits will likely be lost later.')
        self.line()
        self.line('Key:')
        self.startcode()
        self.line('    abc             required (one)')
        self.line('  [ def ]           optional (zero or one)')
        self.line('  ( ghi )*          any number (zero or more)')
        self.line('  ( jkl )+          at least one')
        self.line('  ( mno ){2,inf}    two or more occurrences')
        self.endcode()
        #self.line('<wiki:toc max_depth="3"/>')
    def finish(self):
        # write footer
        pass
    def line(self, *args):
        print(*args, file=self.f)
    def doc(self, name, doc):
        self.line()
        self.heading(4, 'Documentation ('+name+')')
        lines = doc.split('\n')
        for line in lines:
            self.line(line.lstrip())

def parse_appinfo(text):
    result = {}
    for pair in text.split(';'):
        parts = pair.split(':', maxsplit=1)
        if len(parts) == 2:
            result[parts[0].strip()] = parts[1].strip()
        else:
            assert len(pair.strip()) == 0 or die('bad appinfo:', text)
    return result

freq_symbs={
    (1,1) : (' ', ''),
    (0,1) : ('[', ']'),
    (0,'inf') : ('(', ')*'),
    (1,'inf') : ('(', ')+')
}
xsdns = 'http://www.w3.org/2001/XMLSchema'
xsdpre = '{' + xsdns + '}'
def replace_pre(name):
    if name is None:
        return None
    parts = name.split(':', 2)
    if len(parts) == 2 and parts[0] in nsmap:
        return '{' + nsmap[parts[0]] + '}' + parts[1]
    else:
        return name

class XSDType:
    """for built-in types (int, string, etc.)"""
    def __init__(self, name):
        self.name = name
        self.appinfo = {}
        self.doc = None
    def collect_elements(self, elements, stypes, parent):
        pass
    def type_spec(self):
        return self.name
    def has_elts(self):
        return 1        # we don't have child elts but we do have a type
    def has_attrs(self):
        return False
    def write_elt_spec(self, w):
        w.line('    ' + self.name)
    def write_doc(self, w, docname):
        pass

class Node:
    """Base schema node"""
    def __init__(self, node):
        self.doc = None
        self.appinfo = {}
        #print('Parsing', node.tag, node.get('name'))
        for child in node:
            self.read_child(child)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'annotation':
            doc_node = child.find(xsdpre + 'documentation')
            if doc_node is not None:
                self.doc = doc_node.text
            appinfo = child.find(xsdpre + 'appinfo')
            if appinfo is not None and appinfo.text is not None:
                self.appinfo = parse_appinfo(appinfo.text)
        else:
            die('unexpected child tag:', child.tag)

class ComplexType(Node):
    """Represent a type in the schema"""
    def __init__(self, node):
        self.attrs = []
        self.children = None # tree of child elements
        self.elements = [] # flat list of child elements
        self.base_name = None
        Node.__init__(self, node)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'all':
            assert self.children is None
            self.children = ('all', self.read_elts(child))
        elif child.tag == xsdpre + 'sequence':
            assert self.children is None
            self.children = ('seq', self.read_elts(child))
        elif child.tag == xsdpre + 'choice':
            assert self.children is None
            self.children = ('choice', self.read_elts(child))
        elif child.tag == xsdpre + 'complexContent':
            self.read_content(child)
        elif child.tag == xsdpre + 'simpleContent':
            assert self.children is None
            self.children = 'simple'    # this is only really to say we can't have another element
            self.read_content(child)
        elif child.tag == xsdpre + 'attribute':
            self.attrs.append(Attribute(child))
        else:
            Node.read_child(self, child)
    
    def read_elts(self, parent):
        result = []
        for child in parent:
            if child.tag == xsdpre + 'element':
                elt = Element(child)
                result.append(elt)
                self.elements.append(elt)
            elif child.tag == xsdpre + 'choice':
                result.append(('choice', self.read_elts(child)))
            else:
                die('unexpected child of <%s>: <%s>' %(parent.tag, child.tag))
        return result
    
    def read_content(self, parent):
        """Extension types"""
        for child in parent:
            if child.tag == xsdpre + 'extension':
                #Our extensions add attributes and sometimes a (seq) element
                # (where the base type is also seq or has no elements)
                assert self.base_name is None
                self.base_name = replace_pre(child.get('base'))
                for child2 in child:
                    self.read_child(child2)
    
    def collect_elements(self, elements, stypes, parent):
        if hasattr(self, 'base_name'):
            for attr in self.attrs:
                attr.set_type(stypes)
            if self.base_name is None:
                self.base_type = None
            else:
                self.base_type = stypes.get(self.base_name)
                assert self.base_type is not None or die('unknown type',self.base_name)
                self.base_type.collect_elements(elements, stypes, parent)
            del self.base_name
        for child in self.elements:
            child.collect_elements(elements, stypes, parent)
    
    def has_attrs(self):
        return len(self.attrs) > 0
    def has_elts(self):
        have_elts = 2 if len(self.elements) > 0 else 0
        if have_elts == 0 and self.base_type is not None:
            return self.base_type.has_elts()
        return have_elts
    def write_elt_links(self, w):
        #w.line()
        #w.heading(3, 'Elements')
        if self.base_type is not None:
            self.base_type.write_elt_links(w)
        for elt in self.elements:
            w.bulleted(w.link(elt.linkname, elt.name))
    def write_doc(self, w, docname):
        if self.doc is not None:
            w.doc(docname, self.doc)
        if self.base_type is not None:
            self.base_type.write_doc(w, 'base type')
    def write_attr_spec(self, w):
        for attr in self.attrs:
            if attr.use == 'optional' or attr.default:
                w.line('  [', attr.type_spec(),('] DEFAULT VALUE '+str(attr.default) if attr.default is not None else ']'))
            else:
                w.line('   ', attr.type_spec())
    def write_elt_spec(self, w):
        if self.base_type is not None:
            self.base_type.write_elt_spec(w) #TODO: is this correct?
        elif self.children is not None and self.children != 'simple':
            self.write_elt_rec(self.children, w, 0)
    def write_elt_rec(self, node, w, depth):
        if isinstance(node, Element):
            symbs = freq_symbs.get(node.occurs, None)
            if symbs is None:
                symbs = ('(', '){'+str(node.occurs[0])+','+str(node.occurs[1])+'}')
            w.line('| '*depth + symbs[0], '<' + node.name, '... />', symbs[1])
        else:
            assert len(node) == 2
            mode = node[0]
            if mode == 'all':
                w.line('| '*depth + 'IN ANY ORDER:')
            elif mode == 'seq':
                w.line('| '*depth + 'IN THIS ORDER:')
            elif mode == 'choice':
                w.line('| '*depth + 'EXACTLY ONE OF:')
            else:
                assert False
            for child in node[1]:
                self.write_elt_rec(child, w, depth+1)

class Attribute(Node):
    def __init__(self, node):
        self.name = node.get('name')
        self.mode = None
        self.default = node.get('default')
        self.use = node.get('use') # optional or required
        if self.use is not None and self.use not in ('optional', 'required'):
            die('bad use specification on attribute:', node.attrib)
        if (self.default is None) == (self.use is None):
            if self.use is None:
                print('Warning: attribute',node.get('name'),\
                    'does not specify use or default; assuming optional', file=sys.stderr)
                use = 'optional'
            else:
                die('attribute must not specify use and default:', node.attrib)
        self.type_name = replace_pre(node.get('type')) # None if no 'type' attr
        Node.__init__(self, node)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'simpleType':
            for child2 in child:
                assert child2.tag == xsdpre + 'restriction'
                # ignore 'base' attr (don't need)
                assert self.mode is None
                self.mode = 'enum'
                self.vals = []
                for child3 in child2:
                    assert child3.tag == xsdpre + 'enumeration'
                    assert len(child3) == 0
                    self.vals.append(child3.get('value'))
        else:
            Node.read_child(self, child)
    
    def set_type(self, stypes):
        if hasattr(self,'type_name'):
            if self.type_name is None:
                self.t = None
            else:
                self.t = stypes.get(self.type_name)
                assert self.t is not None or die('unknown type',self.type_name)
            del self.type_name
    
    def type_spec(self):
        rst=self.name + '='
        if self.mode == 'enum':
            sep='('
            for val in self.vals:
                rst += sep + '"' + val + '"'
                sep = ' or '
            rst += ')'
        else:
            rst += self.t.type_spec()
        return rst
    
    def writedoc(self, w):
        w.line()
        name = self.appinfo.get('name', self.name)
        w.heading(4, name)
        w.startcode('xml')
        w.line(self.type_spec())
        w.endcode()
        if 'units' in self.appinfo:
            w.line(w.bold('Units:'), self.appinfo['units'])
        if 'min' in self.appinfo:
            w.line(w.bold('Min:'), self.appinfo['min'])
        if 'max' in self.appinfo:
            w.line(w.bold('Max:'), self.appinfo['max'])
        if self.default:
            w.line(w.bold('Default value:'), self.default)
            
        if self.doc is not None:
            w.line()
            lines = self.doc.split('\n')
            for line in lines:
                w.line(line.lstrip())

class Element(Node):
    """Represent an element type in the schema"""
    def __init__(self, node):
        self.name = node.get('name')
        self.elt_type = None
        self.type_name = replace_pre(node.get('type')) # None if no 'type' attr
        maxO = node.get('maxOccurs', 1)
        maxO = 'inf' if maxO == 'unbounded' else int(maxO)
        self.occurs = (int(node.get('minOccurs', 1)), maxO)
        Node.__init__(self, node)
    
    def read_child(self, child):
        if child.tag == xsdpre + 'complexType':
            assert self.elt_type is None and self.type_name is None or die('redundant type specification')
            self.elt_type = ComplexType(child)
        else:
            Node.read_child(self, child)
    
    def depth(self):
        return 0 if self.parent is None else 1 + self.parent.depth()
    
    def collect_elements(self, elements, stypes, parent):
        if self in elements:
            if parent.depth() < self.parent.depth():
                self.parent = parent
        else:
            self.parent = parent
            self.linkname = make_link_from_name(self.name)
            elements.append(self)
            
            if self.elt_type is None:
                self.elt_type = stypes.get(self.type_name)
            assert self.elt_type is not None or die('type not found:', self.type_name)
            self.elt_type.collect_elements(elements, stypes, self)
    
    def writedoc(self, w):
        w.line()
        w.line()
        w.heading(1, '<a name="'+self.linkname+'"></a> ' + self.appinfo.get('name', self.name))
        w.line(self.breadcrumb(w, None))
        w.line()
        #w.heading(5, 'specification')
        w.startcode('xml')
        have_attrs = self.elt_type.has_attrs()
        have_elts = self.elt_type.has_elts()
        w.line('<'+self.name+('' if have_attrs else ('>' if have_elts>0 else '/>')))
        if have_attrs:
            self.elt_type.write_attr_spec(w)
            w.line('  >' if have_elts>0 else '  />')
        if have_elts>0:
            self.elt_type.write_elt_spec(w)
            w.line('</'+self.name+'>')
        w.endcode()
        
        if have_elts>1:
            self.elt_type.write_elt_links(w)
        
        if self.doc is not None:
            w.doc('element', self.doc)
        self.elt_type.write_doc(w, 'type')
        
        if have_attrs:
            w.line()
            w.heading(3, 'Attributes')
            for attr in self.elt_type.attrs:
                attr.writedoc(w)
    
    def breadcrumb(self, w, parent):
        parent = self.parent if parent is None else parent
        r = '→ ' if parent is None else parent.breadcrumb(w, None) + ' → '
        return r + w.link(self.linkname, self.name)

class FixedAttribute:
    def __init__(self, string):
        self.string = string
        self.use = 'required'
        self.default = None
    def set_type(self, stypes):
        pass
    def type_spec(self):
        return self.string
    def writedoc(self, w):
        pass

def translate(f_in, f_out, schema_file):
    global nsmap
    global all_links
    all_links=set() # reset
    
    m = re.match('scenario[^0-9]*([0-9_]+)', schema_file)
    ver = str(m.group(1)) if m is not None else '??'
    w = DocWriter(f_out)
    w.head(schema_file, ver)
    
    tree = ET.parse(f_in)
    root = tree.getroot()
    assert root.tag == xsdpre + 'schema'
    targetns = root.get('targetNamespace')
    ompre = '' if targetns is None else '{' + targetns + '}'
    try:
        nsmap = root.nsmap
    except AttributeError:
        # default ElementTree implementation doesn't let us access this, so guess:
        nsmap = { 'xs' : xsdns }
        if targetns is not None:
            nsmap['om'] = targetns
    
    stypes = {
        xsdpre + 'boolean' : XSDType('boolean'),
        xsdpre + 'int' : XSDType('int'),
        xsdpre + 'integer' : XSDType('integer'),
        xsdpre + 'decimal' : XSDType('decimal'),
        xsdpre + 'double' : XSDType('double'),
        xsdpre + 'string' : XSDType('string')
    }
    omroot = None
    for child in root:
        if child.tag == xsdpre + 'complexType':
            stypes[ompre + child.get('name')] = ComplexType(child)
        elif child.tag == xsdpre + 'element':
            assert omroot is None or die('schema contains multiple root elements')
            omroot = Element(child)
        else:
            die('unexpected tag in schema:', child.tag)
    
    # handle namespace crap
    omroot.elt_type.attrs.append(FixedAttribute(
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'))
    if targetns is None: # pre 32
        omroot.elt_type.attrs.append(FixedAttribute(
            'xsi:noNamespaceSchemaLocation="'+schema_file+'"'))
    else:
        omroot.elt_type.attrs.append(FixedAttribute(
            'xmlns:om="http://openmalaria.org/schema/scenario_' + ver + '"'))
        omroot.elt_type.attrs.append(FixedAttribute(
            'xsi:schemaLocation="http://openmalaria.org/schema/scenario_'+ver+' '+schema_file+'"'))
    
    elements=[]
    omroot.collect_elements(elements, stypes, None)
    
    for elt in elements:
        elt.writedoc(w)
    
    w.finish()

def die(*args):
    print('Error: ', *args, file=sys.stderr)
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="""This tool converts an
                        OpenMalaria schema file to a more readable documentation format.""")
    parser.add_argument('schema', metavar='SCHEMA', nargs='+',
                        help="Schema file to be translated")
    parser.add_argument('-o','--output', metavar='FILE', action='store',
                        help='File to output to (if not given, output is to stdout')
    args = parser.parse_args()
    # args.schema : file to translate
    # args.output : output file
    global linkbase
    if len(args.schema) == 1:
        schema_name = os.path.basename(args.schema[0])
        with open(args.schema[0], 'r') as f_in:
            if args.output is not None:
                linkbase = os.path.splitext(os.path.basename(args.output))[0]
                with open(args.output, 'w') as f_out:
                    translate(f_in, f_out, schema_name)
            else:
                translate(f_in, sys.stdout, schema_name)
    else:
        assert len(args.schema) > 1
        assert args.output is None or die('cannot use -o / --output option with more than one input file')
        generated=[]
        for sch_file in args.schema:
            print('Reading',sch_file,'...', file=sys.stderr)
            schema = os.path.basename(sch_file)
            m = re.match('scenario_([0-9_]+)\.xsd', schema)
            assert m is not None or die('Schema file name does not match regular expression:',schema)
            linkbase = 'GeneratedSchema'+str(m.group(1))+'Doc'
            generated.append((linkbase,schema))
            with open(sch_file, 'r') as f_in:
                with open(linkbase + '.md', 'w') as f_out:
                    translate(f_in, f_out, schema)
        with open('GeneratedSchemaDocIndex.md', 'w') as f_out:
            w = DocWriter(f_out)
            w.heading(1, 'Generated Schema Documentation')
            w.line('This documentation was generated with the following command.')
            w.line('Comments welcome but be warned edits will most likely be lost.')
            w.startcode('sh')
            w.line(' '.join(sys.argv))
            w.endcode()
            w.line()
            w.heading(2, 'Index')
            for link, schema in generated:
                w.bulleted(w.link(link, schema, internal=False))
            w.finish()

if __name__ == "__main__":
    main()
