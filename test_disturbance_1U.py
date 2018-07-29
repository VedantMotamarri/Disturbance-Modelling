import numpy as np
import Disturbance_1U_02 as dist
import unittest
from satellite import Satellite
from constants_1U import m_INERTIA
from ddt import data, ddt , unpack

@ddt
class TestGravityGradient(unittest.TestCase):
    #Test the data type of output
    def test_gg_data_type(self):
        sat = Satellite(np.array([1.,0.,0.,0.,0.,0.,0.]),13.)
        sat.setPos(np.array([7070e3,0.,0.]))  
        result = dist.vTGgB(sat)
        self.assertEqual(type(result),np.ndarray)
    
    #Take position vector to be eigen vector of m_INERTIA along with q =[1,0,0,0] thus torque will be zero
    l, w = np.linalg.eig(m_INERTIA)
    @data(7070e3*w[:,0],7070e3*w[:,1],7070e3*w[:,2])
    def test_inertia_eigenvec(self,value):
        state = np.array([1.,0.,0.,0.,0.1,-0.02,-0.2])
        mySat = Satellite(state,128.05)
        mySat.setPos(value)
        mySat.setVel(np.array([5.60,-5.0,0.0]))
        result = dist.vTGgB(mySat)
        print result
        self.assertTrue(np.allclose(result,[0.,0.,0.]))

class TestAeroDrag(unittest.TestCase):
    #test the data type of the output
    
    def test_aero_type(self):
        sat = Satellite(np.array([np.sqrt(0.5),-np.sqrt(0.5),0.,0.,0.1,0.23,0.]),13.)
        sat.setPos(np.array([0.,0.,7e6]))
        
        sat.setVel(np.array([0,2e3,6e3]))
        result = dist.vTAdB(sat)
       
        self.assertEqual(type(result),np.ndarray)
    def test_aero_value(self):
        sat = Satellite(np.array([np.sqrt(0.5),-np.sqrt(0.5),0.,0.,0.1,0.23,0.]),13.)
        sat.setPos(np.array([0.,0.,7e6]))
        sat.setVel(np.array([0,2e3,6e3]))
        result = dist.vTAdB(sat)
        print result
        self.assertTrue(np.allclose(result, [2.99654080e-10,-2.57065600e-11,-7.71196800e-11]))
       
         
class TestSolarTorque(unittest.TestCase):
        #test the data type of the output
    def test_solar_torque_type(self):
        sat = Satellite(np.array([1.,0.,0.,0.,0.,0.,0.]),13.)
        sat.setPos(np.array([7070e3,0.,0.]))
        v_sv_i=np.array([1.0,0.0,0.0])
        sat.setSun_i(v_sv_i)
        sat.setLight(0)  
        result = dist.vTSdB(sat)
        self.assertEqual(type(result),np.ndarray)
    
    def test_solar_torque_value(self):
        state = np.array([1.,0.,0.,0.,0.1,-0.02,-0.2])
        mySat = Satellite(state,128.05)
        v_sv_i=np.array([1.0,0.0,0.0])            #sun vector in eci frame
        mySat.setSun_i(v_sv_i)
        mySat.setLight(0)
        result = dist.vTSdB(mySat)
        print result
        self.assertTrue(np.allclose(result,[  0.00000000e+00,-3.66624000e-11,3.17376000e-10]))

if __name__=="__main__":
    unittest.main(verbosity=2)
