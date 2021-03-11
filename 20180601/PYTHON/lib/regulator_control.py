import opendssdirect as dss
import numpy as np
from lib.helper import voltage_regulator_index_dict

def regulator_control(XNR, vr_lines, tf_lines, nnode, nline):
# returns associated indices (values as list) of a bus's voltage regulators (keys)
    vr_idx_dict = voltage_regulator_index_dict() 
    vr_line_idx = range(0, vr_lines)
   
    # flag if need to rerun NR3
    flag = 0
    vr_line_counter = 0
    XNR_final = XNR
   
    for k in vr_idx_dict.keys():
        for vridx in vr_idx_dict[k]: #{bus: [indices in dss.RegControls.AllNames(), ...]}
            
            dss.RegControls.Name(dss.RegControls.AllNames()[vridx])
            dss.Circuit.SetActiveBus(dss.CktElement.BusNames()[0].split(".")[0])
            winding = dss.RegControls.Winding()
        
            Vbase = dss.Bus.kVBase() *1000
            Sbase = 10**6
            Ibase = Sbase / Vbase
            band = dss.RegControls.ForwardBand()
            target_voltage = dss.RegControls.ForwardVreg()

            idxbs = dss.Circuit.AllBusNames().index((dss.CktElement.BusNames()[0].split('.')[0]))   
        
            ph = dss.CktElement.BusNames()[0].split('.')[1:] 
            ph_arr = [0, 0, 0]
            for i in ph:
                ph_arr[int(i) - 1] = 1
            if len(ph) == 0:
                ph_arr = [1, 1, 1]
            for ph in range(len(ph_arr)):
                if ph_arr[ph] == 1: #loop over existing phases of voltage regulator
                    
                    NR_voltage = np.abs((XNR[2*nnode*ph + 2*idxbs]  + (1j * XNR[2*nnode*ph + 2*idxbs + 1])) * Vbase / dss.RegControls.PTRatio())

                    if dss.RegControls.ForwardR() and dss.RegControls.ForwardX() and dss.RegControls.CTPrimary():
                        # if LDC exists

                        #vr_line_counter - counts the number of lines passed; two lines for every phase
                        #vridx - index of current voltage regulator in dss.RegControls.AllNames()
                        #tf_lines - number of transformers
                        
                        line_idx =  2*vr_line_idx[vr_line_counter] + 2*(winding - 1)

                        I_reg = XNR[2*3*(nnode+nline) + 2*tf_lines + line_idx] + \
                            1j * XNR[2*3*(nnode+nline) + 2*tf_lines +  line_idx + 1]

                        V_drop = (dss.RegControls.ForwardR() + 1j*dss.RegControls.ForwardX()) / 0.2 * (I_reg * Ibase / dss.RegControls.CTPrimary())     
                    
                        V_drop = (dss.RegControls.ForwardR() + 1j*dss.RegControls.ForwardX()) / 0.2 * (I_reg * Ibase / (dss.RegControls.CTPrimary()/0.2))
                        V_R = np.abs(NR_voltage - V_drop)

                        abs_diff = np.abs(V_R - target_voltage)
                        V_compare = V_R
                        print('vcompare', dss.RegControls.Name(), V_compare)
                                
                    else:
                        # if LDC term does not exist
                        print('**** LDC missing term ***** ')
                        abs_diff = np.abs(NR_voltage - target_voltage)
                        V_compare = NR_voltage
  
                    print('absolute difference: ', abs_diff, "\n")
                    vr_line_counter += 1
                  
                    # compare NR3 voltage to forward Vreg voltage +- band
                    if abs_diff <= band: #converges
                        XNR_final = XNR
                        continue
                            
                    elif abs_diff > band:
                        if V_compare > (target_voltage + band): #NR3 voltage above forward-Vreg
                            if dss.RegControls.TapNumber() <= -16 :
                                print('Tap Number Out of Bounds' )
                                XNR_final = XNR
                        
                            else:
                                print('Decrease Tap Number')
                                dss.RegControls.TapNumber(dss.RegControls.TapNumber() - 1)
                                print('New tap number ', dss.RegControls.TapNumber())
                                flag = 1 #run NR3 again
                        else: #NR3 voltage below forward-Vreg
                            if dss.RegControls.TapNumber() >= 16:
                                print('Tap Number Out of Bounds' )
                                print('New tap number ', dss.RegControls.TapNumber())
                                XNR_final = XNR
                        
                            else:
                                print('Increase tap number')
                                dss.RegControls.TapNumber(dss.RegControls.TapNumber() + 1)
                                flag = 1 #run NR3 again
        if flag == 1:
            return flag, XNR_final
    return flag, XNR_final