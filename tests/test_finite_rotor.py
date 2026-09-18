import math
import unittest
import numpy as np
from gravcomm import FiniteRotor,TwoMassRotor,Actuator,simulate_actuator


class FiniteRotorTest(unittest.TestCase):
    def setUp(self):
        self.rotor=FiniteRotor(np.array([[.2,0,0],[-.2,0,0]]),np.array([2.,2.]),.2)

    def test_point_mass_limit_and_harmonic_sign(self):
        point=TwoMassRotor.balanced_quadrupole(4.,.2)
        phi=np.linspace(0,2*np.pi,137)
        np.testing.assert_allclose(self.rotor.acceleration_component(.5,phi),point.acceleration_component(.5,phi),rtol=1e-14)
        np.testing.assert_allclose(self.rotor.harmonic_coefficient(.5),point.harmonic_coefficient(.5),rtol=1e-13)
        self.assertAlmostEqual(self.rotor.polar_inertia,.16)
        self.assertAlmostEqual(self.rotor.quadrupole_anisotropy.real,.16)

    def test_axisymmetric_ring_has_negligible_ac(self):
        phi=np.arange(256)*2*np.pi/256
        ring=FiniteRotor(np.column_stack([.2*np.cos(phi),.2*np.sin(phi),np.zeros(256)]),np.full(256,4/256),.2)
        self.assertLess(ring.harmonic_amplitude(1),1e-24)
        self.assertGreater(ring.polar_inertia,.15)

    def test_uniform_sphere_against_external_point_mass_field(self):
        # Independent spherical-coordinate quadrature and Newton's sphere theorem.
        x,w=np.polynomial.legendre.leggauss(8);radius=.04
        r=(x+1)*radius/2;wr=w*radius/2
        mu=x;az=np.arange(32)*2*np.pi/32
        points=[];weights=[]
        for ri,wi in zip(r,wr):
            for z,wz in zip(mu,w):
                for phi in az:
                    points.append([.2+ri*np.sqrt(1-z*z)*np.cos(phi),ri*np.sqrt(1-z*z)*np.sin(phi),ri*z])
                    weights.append(wi*ri*ri*wz*2*np.pi/32)
        masses=np.array(weights)*2/(4*np.pi*radius**3/3)
        sphere=FiniteRotor(np.array(points),masses,.24)
        point=TwoMassRotor.asymmetric_dipole(2,.2)
        phases=np.array([0,.3,1.2,2.4])
        np.testing.assert_allclose(sphere.acceleration_component(.5,phases),point.acceleration_component(.5,phases),rtol=1e-11)

    def test_orientation_phase_and_axial_symmetry(self):
        a=.37;p=self.rotor.positions_m.copy()
        p[:,0]=self.rotor.positions_m[:,0]*np.cos(a);p[:,1]=self.rotor.positions_m[:,0]*np.sin(a)
        shifted=FiniteRotor(p,self.rotor.masses_kg,.2)
        np.testing.assert_allclose(shifted.harmonic_coefficient(.5),self.rotor.harmonic_coefficient(.5)*np.exp(2j*a),rtol=1e-13)
        np.testing.assert_allclose(shifted.acceleration_component(.5,[0,1,2],'z'),0,atol=1e-30)

    def test_actuator_uses_achieved_phase_for_finite_source(self):
        drive=Actuator(1,10,5,.02,100,100,10)
        t=np.linspace(0,1,101)
        trace=simulate_actuator(drive,t,np.ones(101),1)
        np.testing.assert_allclose(trace.acceleration(self.rotor,.5,harmonic=None),
            self.rotor.acceleration_component(.5,trace.phase_rad),rtol=1e-14)

    def test_invalid_and_immutable_quadrature(self):
        with self.assertRaises(ValueError):FiniteRotor([[0,0,0]],[-1],1)
        with self.assertRaises(ValueError):self.rotor.acceleration_component(.1,0)
        with self.assertRaises(ValueError):self.rotor.positions_m[0,0]=0


if __name__=='__main__':unittest.main()
