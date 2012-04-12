#!/usr/bin/python
#Write parameter documenation to an Excel File. You may need to
#install http://www.python-excel.org/ xlwt, e.g.:
#apt-get install python-xlwt
#or download the windows installer

from xlwt import Workbook
from xlwt import easyxf
from xml.etree.ElementTree import ElementTree

import string
import re

appinfoOrder=['name','units','min','max','exposed','sweepable']
docu_xf = easyxf('align: wrap on, vert centre')
worksheets = {}

def normalizeNewlines(string):
    string=" ".join(string.split())
    return re.sub(r'(\r\n|\r|\n)', ' ', string)

def formatElement(el,path):
    if (el.find("{http://www.w3.org/2001/XMLSchema}annotation") is not None):
        splitPath=string.split(path,'/')
        #set default component to "scenario", i.e. write to the "Scenario" worksheet below
        component="scenario"
        printPath=string.replace(path,"/"+component,"")
        #update component if a know child element of scenario, i.e. write to another worksheet below
        if (len(splitPath)>2):
            if (splitPath[2] in worksheets):
                component=splitPath[2]
                printPath=string.replace(printPath,"/"+component,"")
        sheet= worksheets[component][0]
        row=worksheets[component][1]
        sheet.write(row,0,string.lstrip(printPath,"/"))
        annotation=el.find("{http://www.w3.org/2001/XMLSchema}annotation")
        docuElem=annotation.find("{http://www.w3.org/2001/XMLSchema}documentation")
        if (docuElem is not None):
            docu=docuElem.text
        else:
            docu="TODO"
        content=string.strip(docu)
        sheet.write(row,1,normalizeNewlines(content),docu_xf)
        appInfoElem=el.find("{http://www.w3.org/2001/XMLSchema}annotation").find("{http://www.w3.org/2001/XMLSchema}appinfo")
        if (appInfoElem is not None):
            appInfo=string.strip(appInfoElem.text)
        else:
            appInfo="name:TODO"
        appInfoList=string.split(appInfo,";")[0:-1]
        for keyValue in appInfoList:
            splitPair=string.split(keyValue,":")
            colIndex=appinfoOrder.index(string.strip(str(splitPair[0])))+2
            sheet.write(row,colIndex,splitPair[1],docu_xf)
        #update next row to be written in that sheet
        worksheets[component][1]=worksheets[component][1]+1
        
def drillDown(el,path,isExtType):
    name=el.get("name")
    if (name and not isExtType):
        path=path+"/"+name
        formatElement(el,path)
    for elem in el.getchildren():
        drillDown(elem,path,False)
    elType=el.get("type")
    if (elType):
        for typeDefinition in tree.findall("{http://www.w3.org/2001/XMLSchema}complexType"):
            if (typeDefinition.get("name")==elType):
                drillDown(typeDefinition,path,True)
                break      

def initWorksheets(elementName,sheetName):
    sheet=book.add_sheet(sheetName)
    worksheets[elementName]=[sheet,1]
    sheet.col(0).width =10000
    sheet.col(1).width =15000
    sheet.col(2).width =6500
    sheet.col(3).width =4300
    sheet.col(4).width =4500
    sheet.col(5).width =4800
    sheet.col(6).width =6000
    sheet.col(7).width =5500
    sheet.write(0,0,"Parameter name and path (in XML document)")
    sheet.write(0,1,"Parameter documentation       ")
    sheet.write(0,2,"Parameter name (in GUI)")
    sheet.write(0,3,"Parameter units")
    sheet.write(0,4,"Parameter min value")
    sheet.write(0,5,"Parameter max value")
    sheet.write(0,6,"Parameter exposed in GUI")
    sheet.write(0,7,"Parameter is sweepable")
    
def main():
    global tree,book
    book = Workbook()
    initWorksheets("scenario", "Scenario")
    initWorksheets("demography", "Demography")
    initWorksheets("monitoring", "Monitoring")
    initWorksheets("interventions", "Interventions")
    initWorksheets("healthSystem", "HealthSystem")
    initWorksheets("entomology", "Entomology")
    initWorksheets("pharmacology", "Pharmacology")
    initWorksheets("model", "Model")
    tree = ElementTree()
    tree.parse("../schema/scenario_30.xsd")
    #we know that the first element in the schema defines scenario
    scenarioElement=tree.find("{http://www.w3.org/2001/XMLSchema}element")
    drillDown(scenarioElement,"",False)
    book.save('Documentation.xls')

if __name__ == '__main__':
    main()
