mkdir test\sandbox
cd test\sandbox
del checkpoint
copy ..\original\densities.csv .
copy ..\original\scenario_5.xsd .
copy ..\original\scenario* .
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario1.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario1.xml
fc ..\original\original1.txt output.txt
copy output.txt output1.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario2.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario2.xml
fc ..\original\original2.txt output.txt
copy output.txt output2.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario3.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario3.xml
fc ..\original\original3.txt output.txt
copy output.txt output3.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario4.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario4.xml
fc ..\original\original4.txt output.txt
copy output.txt output4.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario5.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario5.xml
fc ..\original\original5.txt output.txt
copy output.txt output5.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario6.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario6.xml
fc ..\original\original6.txt output.txt
copy output.txt output6.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario7.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario7.xml
fc ..\original\original7.txt output.txt
copy output.txt output7.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario8.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario8.xml
fc ..\original\original8.txt output.txt
copy output.txt output8.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario9.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario9.xml
fc ..\original\original9.txt output.txt
copy output.txt output9.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario10.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario10.xml
fc ..\original\original10.txt output.txt
copy output.txt output10.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario11.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario11.xml
fc ..\original\original11.txt output.txt
copy output.txt output11.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario12.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenario12.xml
fc ..\original\original12.txt output.txt
copy output.txt output12.txt
del checkpoint
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenarioIPT.xml
..\..\malariacontrol_6.17_windows_intelx86 --checkpoint --scenario scenarioIPT.xml
fc ..\original\originalIPT.txt output.txt
copy output.txt outputIPT.txt
del checkpoint
pause
