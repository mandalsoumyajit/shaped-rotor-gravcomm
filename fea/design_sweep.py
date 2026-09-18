"""Fixed-mass geometry screening followed by selected structural evaluations.

Maximize radial second-harmonic gravity at a sensor in the rotor plane 1 m
from the axis. Constraints: 7.9135094916 kg total, 250 mm radius, 90 mm axial
envelope, 1800 rpm. This is a sampled trade study, not a global optimum.
"""
import csv
import itertools
import json
from pathlib import Path

import gmsh
import numpy as np
from scipy.optimize import brentq
from shaped_rotor import solve
from postprocess_carrier import quadrature, read_input, write_csv

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'/'design_sweep'
MASS=7.913509491619677
G=6.67430e-11


def build(spec,insert_size):
    """Populate active Gmsh model and return material-tag sets and CAD integrals."""
    gmsh.model.add('candidate')
    occ=gmsh.model.occ
    c=spec['center'];w=spec['spoke'];rim=spec['rim'];t=.02
    sphere=spec['shape']=='sphere'
    r=insert_size if sphere else spec['radius']
    h=2*r if sphere else insert_size
    def cyl(x,r,z=-t/2,height=t):return (3,occ.addCylinder(x,0,z,0,0,height,r))
    outer,inner=cyl(0,.25),cyl(0,.25-rim)
    ring,_=occ.cut([outer],[inner])
    R=.25-rim+.001
    parts=[cyl(0,.03),cyl(c,r+.015),cyl(-c,r+.015),
           (3,occ.addBox(-R,-w/2,-.01,2*R,w,.02)),
           (3,occ.addBox(-w/2,-R,-.01,w,2*R,.02))]
    carrier,_=occ.fuse(ring,parts)
    fillet=spec.get('fillet',0.)
    if fillet:
        occ.synchronize()
        adjacent={}
        for dim,surface in gmsh.model.getBoundary(carrier,oriented=False):
            for _,edge in gmsh.model.getBoundary([(dim,surface)],oriented=False):
                adjacent.setdefault(edge,set()).add(surface)
        edges=[]
        for edge,surfaces in adjacent.items():
            bounds=gmsh.model.getBoundingBox(1,edge)
            if (len(surfaces)==2 and gmsh.model.getType(1,edge)=='Line'
                and bounds[3]-bounds[0]<1e-6 and bounds[4]-bounds[1]<1e-6
                and abs(bounds[5]-bounds[2]-.02)<1e-6):
                edges.append(edge)
        if not edges:raise ValueError('No carrier junctions selected for filleting')
        carrier=occ.fillet([tag for dim,tag in carrier],edges,[fillet])
        assert all(dim==3 for dim,tag in carrier)
    carrier,_=occ.cut(carrier,[cyl(0,.006,z=-.011,height=.022)])
    def insert(x):
        return (3,occ.addSphere(x,0,0,r)) if sphere else cyl(x,r,z=-h/2,height=h)
    # Subtract geometrically exact insert volume, then replace with steel.
    carrier,_=occ.cut(carrier,[insert(c),insert(-c)])
    steel=[insert(c),insert(-c)]
    _,mapping=occ.fragment(carrier,steel)
    alu={tag for row in mapping[:len(carrier)] for dim,tag in row if dim==3}
    steel={tag for row in mapping[len(carrier):] for dim,tag in row if dim==3}
    assert alu.isdisjoint(steel)
    occ.synchronize()
    stats={};quad=0.
    for name,tags,rho in [('aluminum',alu,2700.),('steel',steel,7850.)]:
        mass=inertia=0.
        for tag in tags:
            v=occ.getMass(3,tag);com=occ.getCenterOfMass(3,tag)
            J=np.asarray(occ.getMatrixOfInertia(3,tag)).reshape(3,3)
            mass+=rho*v;inertia+=rho*(J[2,2]+v*(com[0]**2+com[1]**2))
            quad+=rho*(J[1,1]-J[0,0]+v*(com[0]**2-com[1]**2))
        stats[name+'_mass_kg']=mass;stats[name+'_polar_inertia_kg_m2']=inertia
    stats.update(quadrupole_x2_minus_y2_kg_m2=quad,insert_radius_m=r,insert_height_m=h,
                 outer_extent_m=max(.25,c+r+.015))
    return alu,steel,stats


def design(spec,enforce_envelope=True):
    gmsh.initialize();gmsh.option.setNumber('General.Terminal',0)
    try:
        def residual(size):
            gmsh.clear();_,_,stats=build(spec,size)
            return stats['aluminum_mass_kg']+stats['steel_mass_kg']-MASS
        if spec['shape']=='sphere':
            size=brentq(residual,.02,.06,xtol=1e-10)
        else:
            # Carrier geometry is independent of height once insert spans plate.
            r=spec['radius'];residual(.09)
            _,_,stats=build_after_clear(spec,.09)
            size=(MASS-stats['aluminum_mass_kg'])/(2*7850*np.pi*r*r)
        gmsh.clear();_,_,stats=build(spec,size)
        if enforce_envelope and (stats['outer_extent_m']>.25000001 or stats['insert_height_m']>.09000001 or stats['insert_height_m']<.020001):
            return None
        assert abs(stats['aluminum_mass_kg']+stats['steel_mass_kg']-MASS)<1e-6
        return dict(**spec,insert_size=size,**stats)
    finally:gmsh.finalize()


def build_after_clear(spec,size):
    gmsh.clear();return build(spec,size)


def mesh_builder(candidate):
    def mesh(h,folder):
        gmsh.initialize()
        try:
            gmsh.option.setNumber('General.Terminal',0);gmsh.option.setNumber('General.NumThreads',1)
            alu,steel,stats=build(candidate,candidate['insert_size'])
            gmsh.option.setNumber('Mesh.MeshSizeMax',h);gmsh.option.setNumber('Mesh.MeshSizeMin',.003)
            for dim,tag in gmsh.model.getEntities(0):
                p=gmsh.model.getValue(dim,tag,[])
                gmsh.model.mesh.setSize([(dim,tag)],.003 if np.hypot(p[0],p[1])<.007 else h)
            if candidate.get('fillet',0.):
                gmsh.option.setNumber('Mesh.MeshSizeMin',min(.003,h/12))
                for dim,tag in gmsh.model.getEntities(1):
                    box=gmsh.model.getBoundingBox(dim,tag)
                    if gmsh.model.getType(dim,tag)=='Circle' and max(np.array(box[3:])-box[:3])<.02:
                        pts=gmsh.model.getBoundary([(dim,tag)],oriented=False)
                        gmsh.model.mesh.setSize(pts,h/12)
            if candidate.get('faceted_geometry',False):
                gmsh.option.setNumber('Mesh.SecondOrderLinear',1)
                # Resolve finite insert/boss circles without imposing tiny
                # curvature sizes over the entire solid volume.
                for dim,tag in gmsh.model.getEntities(1):
                    box=gmsh.model.getBoundingBox(dim,tag)
                    span=max(box[3]-box[0],box[4]-box[1])
                    if gmsh.model.getType(dim,tag)=='Circle' and .02<span<.12:
                        length=gmsh.model.occ.getMass(dim,tag)
                        gmsh.model.mesh.setTransfiniteCurve(tag,max(3,int(np.ceil(length/.004))+1))
            gmsh.model.mesh.generate(3);gmsh.model.mesh.setOrder(2)
            # Curved small-radius fillets can invert midside-node tetrahedra.
            # Check the solver quadrature before accepting/exporting the mesh.
            def min_jacobian():
                values=[]
                for typ in gmsh.model.mesh.getElementTypes(3):
                    points,_=gmsh.model.mesh.getIntegrationPoints(typ,'Gauss2')
                    _,det,_=gmsh.model.mesh.getJacobians(typ,points)
                    values.append(float(min(det)))
                return min(values)
            if min_jacobian()<=0:
                gmsh.model.mesh.optimize('HighOrder')
            if min_jacobian()<=0:raise RuntimeError('Invalid curved tetrahedra remain after high-order optimization')
            gmsh.write(str(folder/'mesh.msh'));gmsh.write(str(folder/'geometry.brep'))
            tags,xyz,_=gmsh.model.mesh.getNodes();nodes=dict(zip(map(int,tags),xyz.reshape(-1,3)))
            fixed=[]
            for dim,tag in gmsh.model.getEntities(2):
                box=gmsh.model.getBoundingBox(dim,tag)
                if (gmsh.model.getType(dim,tag)=='Cylinder' and max(abs(box[0]),abs(box[1]),abs(box[3]),abs(box[4]))<.006001):
                    fixed.extend(map(int,gmsh.model.mesh.getNodes(dim,tag,includeBoundary=True)[0]))
            stats['_fixed_nodes']=sorted(set(fixed))
            assert len(stats['_fixed_nodes'])>12
            elements={'ALUMINUM':[],'STEEL':[]}
            for name,volumes in [('ALUMINUM',alu),('STEEL',steel)]:
                for vol in volumes:
                    types,ids,conns=gmsh.model.mesh.getElements(3,vol)
                    for typ,eids,conn in zip(types,ids,conns):
                        _,dim,_,n,local,primary=gmsh.model.mesh.getElementProperties(typ)
                        assert n==10 and primary==4
                        local=np.asarray(local).reshape(n,dim);order=list(range(4))
                        for i,j in ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3)):
                            d=np.linalg.norm(local-(local[i]+local[j])/2,axis=1);assert min(d)<1e-12
                            order.append(int(np.argmin(d)))
                        elements[name]+=[(int(e),list(map(int,c[order]))) for e,c in zip(eids,conn.reshape(-1,n))]
            return nodes,elements,stats
        finally:gmsh.finalize()
    return mesh


def gravity(folder,nphase=128):
    nodes,elements,materials=read_input(folder/'model.inp')
    q=quadrature(nodes,elements,materials)
    p=np.array([v[3] for v in q]);mass=np.array([v[4]*v[5] for v in q])
    phi=np.arange(nphase)*2*np.pi/nphase
    result={}
    for d in (.5,1.,2.):
        values=[]
        for angle in phi:
            x=p[:,0]*np.cos(angle)-p[:,1]*np.sin(angle)
            y=p[:,0]*np.sin(angle)+p[:,1]*np.cos(angle)
            dx=x-d;r2=dx*dx+y*y+p[:,2]**2
            values.append(G*np.sum(mass*dx/r2**1.5))
        c=2*np.mean(np.array(values)*np.exp(-2j*phi))
        result[f'a2_at_{d:g}m_m_s2']=float(abs(c))
    result['quadrature_mass_kg']=float(sum(mass))
    return result


def main():
    OUT.mkdir(exist_ok=True)
    baseline=design(dict(shape='cylinder',radius=.03,center=.2,spoke=.025,rim=.035))
    candidates=[]
    # Balanced coarse grid: all candidates keep thickness, hub, boss wall and
    # bore fixed. Shape/radius/placement/width are explicit bounded variables.
    for shape,radius,center,spoke,rim in itertools.product(['cylinder','sphere'],[.03,.035],[.18,.2,.21],[.02,.03,.04],[.02,.035]):
        if shape=='sphere' and radius!=.03:continue
        spec=dict(shape=shape,radius=radius,center=center,spoke=spoke,rim=rim)
        if shape=='cylinder' and center+radius+.015>.25000001:continue
        item=design(spec)
        if item:candidates.append(item)
    candidates.sort(key=lambda r:r['quadrupole_x2_minus_y2_kg_m2'],reverse=True)
    write_csv(OUT/'cad_screen.csv',candidates)
    # Select distinct useful directions, not every near-duplicate of the best Q.
    chosen={'baseline':baseline,'best_quadrupole':candidates[0],
            'best_sphere':next(r for r in candidates if r['shape']=='sphere'),
            'wide_spokes':next(r for r in candidates if r['spoke']==.04)}
    rows=[]
    for name,candidate in chosen.items():
        print('Evaluating '+name+': '+json.dumps(candidate),flush=True)
        row=solve(.018,1800.,OUT/name,mesh_builder(candidate),spec=candidate)
        folder=OUT/name/'h0.018_rpm1800'
        row.update(gravity(folder));row.update(case=name,shape=candidate['shape'],center_m=candidate['center'],
                spoke_width_m=candidate['spoke'],rim_width_m=candidate['rim'])
        rows.append(row);write_csv(OUT/'structural_screen.csv',rows)
    (OUT/'assumptions.json').write_text(json.dumps(dict(total_mass_kg=MASS,radius_limit_m=.25,
        axial_envelope_m=.09,rpm=1800,objective='max radial second harmonic at 1 m; CAD quadrupole pre-screen',
        structural_constraints='screen only: stress and modal margins need final design requirements',
        baseline=baseline,selected=chosen),indent=2))


if __name__=='__main__':main()
