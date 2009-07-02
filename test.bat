mkdir test\sandbox
cd test\sandbox
del checkpoint
copy ..\original\densities.csv .
copy ..\original\scenario_5.xsd .
copy ..\original\scenario* .
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario1.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario1.xml
fc ..\original\original1.txt output.txt
rename output.txt output1.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario2.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario2.xml
fc ..\original\original2.txt output.txt
rename output.txt output2.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario3.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario3.xml
fc ..\original\original3.txt output.txt
rename output.txt output3.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario4.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario4.xml
fc ..\original\original4.txt output.txt
rename output.txt output4.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario5.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario5.xml
fc ..\original\original5.txt output.txt
rename output.txt output5.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario6.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario6.xml
fc ..\original\original6.txt output.txt
rename output.txt output6.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario7.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario7.xml
fc ..\original\original7.txt output.txt
rename output.txt output7.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario9.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario9.xml
fc ..\original\original9.txt output.txt
rename output.txt output9.txt
del checkpoint
rem ..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario10.xml
rem ..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario10.xml
rem fc ..\original\original10.txt output.txt
rem rename output.txt output10.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario11.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario11.xml
fc ..\original\original11.txt output.txt
rename output.txt output11.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario12.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario12.xml
fc ..\original\original12.txt output.txt
rename output.txt output12.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenarioIPT.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenarioIPT.xml
fc ..\original\originalIPT.txt output.txt
rename output.txt outputIPT.txt
del checkpoint
rem Checkpointing is broken:
..\..\malariacontrol_6.17_windows_intelx86 --scenario scenarioDummyPKPD.xml
fc ..\original\originalDummyPKPD.txt output.txt
rename output.txt outputDummyPKPD.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --scenario scenarioCevCq.xml
fc ..\original\originalCevCq.txt output.txt
rename output.txt outputDummyCevCq.txt
del checkpoint
pause
