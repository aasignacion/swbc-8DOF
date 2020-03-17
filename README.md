# swbc-8DOF

To run this program, download the SWBC library in http://stanford-wbc.sourceforge.net/doxy/tutorials-page.html. Moreover, the dependecies must be downloaded and then built properly.

After this, copy the four files: Control_code.cpp, tutsim1.cpp, tutsim1.hpp, tutrob1.xml and CMakelist.txt.

Details on SWBC: Stanford whole-body control contains library for whole body control of multi-link robot manipulator. Tutorial C++ program are provided for both joint-space and task-space controller.

Details on this code: This program is to verify the controller proposed in https://www.researchgate.net/publication/316335592_Robust_Impedance_Control_of_High-DOF_Robot_Based_on_ISMC_and_DOB.

This code is based on Tut05 code for task-space control or end-effector control. Modes are created to show the efficacy of the proposed controller.

Mode 0: No control Mode 1: PD control with no disturbance Mode 2: PD control with with disturbance Mode 3: ISMC with disturbance Mode 4: Proposed Control: ISMC + DOB with disturbance

Modes can switched by clicking the toggl3 button on the FLTK GUI.
