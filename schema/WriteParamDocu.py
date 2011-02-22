#Write parameter documenation to an Excel File. You may need to
#install http://www.python-excel.org/ xlwt, e.g.:
#apt-get install python-xlwt

from xlwt import Workbook
from xlwt import easyxf
from xml.etree.ElementTree import ElementTree

import string
import re

appinfoOrder=['name','units','min','max','exposed','sweepable']
docu_xf = easyxf('align: wrap on, vert centre')

def normalizeNewlines(string):
    string=" ".join(string.split())
    return re.sub(r'(\r\n|\r|\n)', ' ', string)

def formatElement(el,path):
    global row,sheet, tree
    if (el.find("{http://www.w3.org/2001/XMLSchema}annotation")):
        sheet.write(row,0,path)
        annotation=el.find("{http://www.w3.org/2001/XMLSchema}annotation")
        docu=annotation.find("{http://www.w3.org/2001/XMLSchema}documentation").text
        content=string.strip(docu)
        sheet.write(row,1,normalizeNewlines(content),docu_xf)
        appInfo=string.strip(el.find("{http://www.w3.org/2001/XMLSchema}annotation").find("{http://www.w3.org/2001/XMLSchema}appinfo").text)
        appInfoList=string.split(appInfo,";")[0:-1]
        for keyValue in appInfoList:
            splitPair=string.split(keyValue,":")
            colIndex=appinfoOrder.index(string.strip(str(splitPair[0])))+2
            sheet.write(row,colIndex,splitPair[1],docu_xf)
        row=row+1
        
def drillDown(el,path):
    name=el.get("name")
    if (name):
        path=path+"/"+name
        formatElement(el,path)
    for elem in el.getchildren():
        drillDown(elem,path)
    elType=el.get("type")
    if (elType):
        for typeDefinition in tree.findall("{http://www.w3.org/2001/XMLSchema}complexType"):
            if (typeDefinition.get("name")==elType):
                drillDown(typeDefinition,path)
                break      
    
def main():
    global row, sheet, tree
    row=1
    book = Workbook()
    sheet = book.add_sheet('Params')
    sheet.write(0,0,"Parameter name and path (in XML document)")
    sheet.write(0,1,"Parameter documentation       ")
    sheet.write(0,2,"Parameter name (in GUI)")
    sheet.write(0,3,"Parameter units")
    sheet.write(0,4,"Parameter min value")
    sheet.write(0,5,"Parameter max value")
    sheet.write(0,6,"Parameter exposed in GUI")
    sheet.write(0,7,"Parameter is sweepable")
    tree = ElementTree()
    tree.parse("scenario.xsd")
    #we know that the first element in the schema defines scenario
    scenarioElement=tree.find("{http://www.w3.org/2001/XMLSchema}element")
    drillDown(scenarioElement,"")
    book.save('Documentation.xls')

if __name__ == '__main__':
    main()
