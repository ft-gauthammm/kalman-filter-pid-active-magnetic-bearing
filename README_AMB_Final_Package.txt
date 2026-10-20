AMB FINAL PACKAGE
1) AMB_PID_KF_Project_V7_AUDITED.m — complete script based on uploaded V6 logic; not executed here.
2) AMB_Mathematical_Model_and_PID_Demo_FINAL.m — companion baseline PID demo.
3) AMB_Final_Mathematical_Model_and_Results_Guide.docx — equations, parameters, results, limitations.
4) AMB_Final_System_Diagram.png / .svg — final system diagram.
5) AMB_V6_Comparison(1).csv, AMB_V6_NoiseSweep(1).csv, AMB_V6_Mismatch(1).csv, AMB_V6_RotorSpeed(1).csv — original uploaded V6 outputs.

MODEL:
m*x_ddot = ks*x + ki*i + Fd
m=46 kg; ki=129.7 N/A; ks=+2.4e6 N/m (positive destabilizing magnitude)
G(s)=X(s)/I(s)=ki/(m*s^2-ks)
poles=+/-sqrt(ks/m)=+/-228.416 rad/s
z=[x;v], zdot=A*z+B*i+E_F*Fd; y=x+noise
A=[[0,1],[ks/m,0]], B=[[0],[ki/m]], E_F=[[0],[1/m]]
dt=0.2 ms, 5 kHz controller; saturation +/-15 A.

IMPORTANT: Single-axis linearized simulation only. The 20 kHz switching frequency is not explicitly modeled. The optimizer is grid search plus local refinement, not UCB. Re-run MATLAB before treating outputs as regenerated.
