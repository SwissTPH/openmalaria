These programs are for converting an XML-based decision tree into a compressed list of outcomes, each associated with a probability.

decisionTree.xml is the source tree, and compressedTree.xml the latest output version (produced by CaseManagementTree.d).

There are two programs here:
 *  CaseManagementTree.d: the latest converter, written in D
 *  CM.java: an older java-based converter
They don't work in the same way and output is not the same; however it has roughly the same format and in both cases is intended to be copied into scenario.xml files for the openmalaria simulator to represent the case-management system.
