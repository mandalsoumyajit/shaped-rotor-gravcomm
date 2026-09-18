"""Run under WSL: python3 -m unittest discover -s .../fea -p test_fea.py."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from run_calculix import parse_dat
from postprocess_carrier import shape, GAUSS, read_displacements, tensor


class DatParserTest(unittest.TestCase):
    def test_static_fields_are_not_mixed_with_normalized_modal_fields(self):
        sample = '''
 forces (fx,fy,fz) for set FIXED and time 0.1000000E+01
 1 -42 0 0
 displacements (vx,vy,vz) for set TIP and time 0.1000000E+01
 2 0.001 0 0
 stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set EALL and time 1
 1 1 100 0 0 0 0 0
 MODE NO EIGENVALUE FREQUENCY
 1 394784.176 628.31853 100.0 0
 P A R T I C I P A T I O N   F A C T O R S
 1 2 3 4 5 6 7
 forces (fx,fy,fz) for set FIXED and time 2
 1 999999 0 0
 stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set EALL and time 2
 1 1 999999 0 0 0 0 0
'''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'test.dat'
            path.write_text(sample)
            blocks,freq = parse_dat(path)
        self.assertEqual(blocks['forces'],[[1.,-42.,0.,0.]])
        self.assertEqual(blocks['stresses'],[[1.,1.,100.,0.,0.,0.,0.,0.]])
        self.assertEqual(freq,[100.])


class PostprocessTest(unittest.TestCase):
    def test_tetra_shape_partition_and_affine_strain(self):
        vertices=np.array([[0.,0,0],[2,0,0],[0,3,0],[0,0,4]])
        xyz=np.vstack([vertices,[(vertices[i]+vertices[j])/2 for i,j in
                                 ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))]])
        A=np.array([[.01,.02,.03],[.04,.05,.06],[.07,.08,.09]])
        u=xyz@A.T
        volume=0.
        for L in GAUSS:
            N,D=shape(L);J=xyz.T@D
            self.assertAlmostEqual(sum(N),1.)
            np.testing.assert_allclose(np.sum(D,axis=0),0,atol=1e-14)
            np.testing.assert_allclose(N@xyz,L@vertices,atol=1e-14)
            np.testing.assert_allclose(u.T@D@np.linalg.inv(J),A,atol=1e-14)
            volume+=np.linalg.det(J)/24
        self.assertAlmostEqual(volume,4.)

    def test_frd_adjacent_signed_numbers_and_modes(self):
        def line(n,values):return ' -1'+f'{n:10d}'+''.join(f'{v:12.5E}' for v in values)+'\n'
        data=' -4  DISP        4    1\n'+line(1,[1.,-2.,3.])+' -3\n'
        data+='    1PMODE                         1\n -4  DISP        4    1\n'+line(1,[-4.,5.,-6.])+' -3\n'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.frd';path.write_text(data);result=read_displacements(path)
        np.testing.assert_array_equal(result[0][1],[1,-2,3])
        np.testing.assert_array_equal(result[1][1],[-4,5,-6])

    def test_stress_tensor_order(self):
        np.testing.assert_array_equal(tensor([1,2,3,4,5,6]),[[1,4,5],[4,2,6],[5,6,3]])

    def test_volume_stress_averages_do_not_cancel_opposing_stress(self):
        from static_regions import stress_averages
        tension=np.diag([100.,0.,0.])
        mean_vm,vm_mean=stress_averages([tension,-tension],[1.,3.])
        self.assertAlmostEqual(mean_vm,100.)
        self.assertAlmostEqual(vm_mean,50.)
        np.testing.assert_allclose(stress_averages([tension+20*np.eye(3)],[2.]),[100.,100.])


if __name__=='__main__':
    unittest.main()
