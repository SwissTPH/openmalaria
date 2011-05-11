#!/usr/bin/python
# -*- coding: utf-8 -*-

# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# 
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import sys
import os.path
from optparse import OptionParser
import xml.dom.minidom as xd
import codecs


def findIncludes(node,inlined,toInline):
    nextElt=node
    while nextElt is not None:
        node=nextElt
        nextElt=node.nextSibling
        if node.nodeType==xd.Node.ELEMENT_NODE and node.tagName=="xs:include":
            incl=node.getAttribute("schemaLocation")
            if not os.path.isfile(incl):
                raise Exception,"schema '{0}' not found".format(incl)
            if not incl in inlined:
                if not incl in toInline:
                    #print "new include: {0}".format(incl)
                    toInline.add(incl)
            node.parentNode.removeChild(node)
            node.unlink()

def inlineSchema(sourcePath,destinationPath):
    """Open sourcePath. For every included XSD, look for the document in the
    current directory, and inline it. Write result to destinationPath."""
    doc=xd.parse(sourcePath)
    docElt=doc.documentElement
    if not (docElt.nodeType==xd.Node.ELEMENT_NODE and docElt.tagName=="xs:schema"):
        raise Exception, "file {0} doesn't start with an elt of name xs:schema".format(sourcePath)
    
    inlined=set()
    toInline=set()
    
    findIncludes(docElt.firstChild,inlined,toInline)
    
    while len(toInline):
        incl=toInline.pop()
        inclDoc=xd.parse(incl)
        iDElt=inclDoc.documentElement
        if not iDElt.nodeType==xd.Node.ELEMENT_NODE and iDElt.tagName=="xs:schema":
            raise Exception, "file {0} doesn't start with an elt of name xs:schema".format(incl)
        findIncludes(iDElt.firstChild,inlined,toInline)
        for child in iDElt.childNodes:
            docElt.appendChild(child.cloneNode(True))
        inclDoc.unlink()
        inlined.add(incl)
        #print "added: {0}".format(incl)
    
    destFile=codecs.open(destinationPath,'w','utf-8')
    docElt.writexml(destFile)
    destFile.close()


def main(args):
    #try:
        parser = OptionParser(usage="Usage: %prog [options] SCHEMA OUT_PATH",
                            description="""Opens SCHEMA, inlines any included schemas (looking for them in the current working directory) and writes to OUT_PATH.""")
        
        (options, others) = parser.parse_args(args=args[1:])
        if not len(others)==2:
            raise Exception,"Expected 2 arguments"
        
        inlineSchema(others[0],others[1])
    #except Exception,e:
    #    print str(e)
    #    return -1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
