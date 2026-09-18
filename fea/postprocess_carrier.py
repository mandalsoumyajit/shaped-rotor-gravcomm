"""Location-matched static stresses and mass-weighted modal diagnostics.

Uses element integration-point stresses, avoiding nodal averaging across the
steel/aluminium interface. C3D10 Gauss order follows CalculiX gauss.f gauss3d5.
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

from run_calculix import parse_dat

ROOT=Path(__file__).resolve().parent
EDGES=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
B=.138196601125011
GAUSS=np.full((4,4),B)
np.fill_diagonal(GAUSS,1-3*B)
DL=np.array([[-1.,-1.,-1.],[1,0,0],[0,1,0],[0,0,1]])


def shape(L):
    N=np.r_[L*(2*L-1),[4*L[i]*L[j] for i,j in EDGES]]
    D=np.vstack([(4*L-1)[:,None]*DL,
                 [4*(L[i]*DL[j]+L[j]*DL[i]) for i,j in EDGES]])
    return N,D


def read_input(path):
    nodes={};elements={};materials={};mode=None;material=None
    for line in path.read_text().splitlines():
        if line.startswith('*'):
            mode=None
            if line.startswith('*NODE,'):mode='node'
            if line.startswith('*ELEMENT,'):
                mode='element';material=line.split('ELSET=')[1].strip()
            continue
        if mode=='node':
            p=line.split(',');nodes[int(p[0])]=np.array(list(map(float,p[1:4])))
        elif mode=='element':
            p=list(map(int,line.split(',')));elements[p[0]]=p[1:];materials[p[0]]=material
    return nodes,elements,materials


def read_displacements(path):
    """FRD fixed-width scientific numbers need not have whitespace separators."""
    result={};mode=0;active=False
    with path.open() as f:
        for line in f:
            if '1PMODE' in line:mode=int(line.split()[-1])
            if line.startswith(' -4'):
                active=line.split()[1]=='DISP'
                if active:result[mode]={}
            elif line.startswith(' -3'):active=False
            elif active and line.startswith(' -1'):
                node=int(line[3:13]);result[mode][node]=np.array([float(line[13+12*k:25+12*k]) for k in range(3)])
    return result


def tensor(s):
    xx,yy,zz,xy,xz,yz=s
    return np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])


def quadrature(nodes,elements,materials):
    records=[]
    for eid,conn in elements.items():
        xyz=np.array([nodes[n] for n in conn])
        rho=2700. if materials[eid]=='ALUMINUM' else 7850.
        for ip,L in enumerate(GAUSS,1):
            N,D=shape(L);J=xyz.T@D;det=np.linalg.det(J)
            if det<=0:raise ValueError(f'Nonpositive Jacobian in {eid}')
            records.append((eid,ip,N,N@xyz,det/24,rho))
    return records


def interface_forces(nodes,elements,materials,stress):
    """Reconstruct linear stress from 4 Gauss samples, integrate on curved faces.

    Seven-point triangular quadrature; reported forces are a postprocessing
    equilibrium check, not a nonlinear contact solution. Normals point OUT of
    aluminium. Force on steel is minus the corresponding aluminium traction.
    """
    faces={}
    for eid,c in elements.items():
        for face in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
            key=tuple(sorted(c[i] for i in face))
            faces.setdefault(key,[]).append((eid,face))
    tri=[([1/3]*3,.225)]
    for a,b,w in ((.059715871789770,.470142064105115,.132394152788506),
                  (.797426985353087,.101286507323456,.125939180544827)):
        for k in range(3):
            v=[b,b,b];v[k]=a;tri.append((v,w))
    rows={sign:dict(force=np.zeros(3),area=0.,normal_min=float('inf'),normal_max=-float('inf')) for sign in (-1,1)}
    for neighbors in faces.values():
        if len(neighbors)!=2 or materials[neighbors[0][0]]==materials[neighbors[1][0]]:continue
        eid,face=next(v for v in neighbors if materials[v[0]]=='ALUMINUM')
        c=elements[eid];xyz=np.array([nodes[n] for n in c])
        stress_corner=np.linalg.solve(GAUSS,np.array([stress[eid,ip] for ip in range(1,5)]))
        opposite=next(i for i in range(4) if i not in face)
        sign=1 if np.mean(xyz[list(face),0])>0 else -1
        r=rows[sign]
        for bary,w in tri:
            L=np.zeros(4);L[list(face)]=bary
            N,D=shape(L)
            # Natural coordinates xi,eta,zeta = L[1:].
            d1=(np.eye(4)[face[1]]-np.eye(4)[face[0]])[1:]
            d2=(np.eye(4)[face[2]]-np.eye(4)[face[0]])[1:]
            normal=np.cross(xyz.T@D@d1,xyz.T@D@d2)
            if normal@(xyz[opposite]-N@xyz)>0:normal=-normal
            jac=np.linalg.norm(normal);unit=normal/jac
            S=tensor(L@stress_corner)
            r['force']-=S@normal*w/2
            r['area']+=jac*w/2
            sn=unit@S@unit
            r['normal_min']=min(r['normal_min'],sn);r['normal_max']=max(r['normal_max'],sn)
    return rows


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def consistent_insert_forces(nodes,elements,materials,displacements):
    """Recover consistent interface nodal force from steel stress and body load.

    Integrates B^T sigma and N^T rho*a using the solver's C3D10 quadrature.
    Avoids extrapolating discontinuous elemental stresses onto a curved boundary.
    Force is on the insert from its bonded carrier, not a contact-pressure field.
    """
    alu_nodes={n for eid,c in elements.items() if materials[eid]=='ALUMINUM' for n in c}
    result={-1:np.zeros(3),1:np.zeros(3)}
    E=200e9;nu=.3;mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu));omega=1800*2*np.pi/60
    for eid,c in elements.items():
        if materials[eid]!='STEEL':continue
        boundary=np.array([n in alu_nodes for n in c])
        if not any(boundary):continue
        xyz=np.array([nodes[n] for n in c]);u=np.array([displacements[n] for n in c])
        sign=1 if np.mean(xyz[:,0])>0 else -1
        for L in GAUSS:
            N,D=shape(L);J=xyz.T@D;grad=D@np.linalg.inv(J);vol=np.linalg.det(J)/24
            du=u.T@grad;strain=(du+du.T)/2;S=2*mu*strain+lam*np.trace(strain)*np.eye(3)
            body=7850*omega**2*(N@xyz)*[1,1,0]
            result[sign]+=vol*np.sum((grad@S-N[:,None]*body)[boundary],axis=0)
    return {str(sign):dict(force_n=v.tolist(),inward_force_n=-sign*v[0]) for sign,v in result.items()}


def analyze(folder,output,plots=False):
    nodes,elements,materials=read_input(folder/'model.inp')
    blocks,freqs=parse_dat(folder/'model.dat')
    stress={(int(v[0]),int(v[1])):np.array(v[2:8]) for v in blocks['stresses']}
    q=quadrature(nodes,elements,materials)
    samples=[]
    for eid,ip,N,p,volume,rho in q:
        if materials[eid]!='ALUMINUM':continue
        S=tensor(stress[eid,ip]);rad=np.hypot(p[0],p[1]);angle=np.arctan2(p[1],p[0])
        tangent=np.array([-np.sin(angle),np.cos(angle),0.])
        vm=np.sqrt(((S[0,0]-S[1,1])**2+(S[1,1]-S[2,2])**2+(S[2,2]-S[0,0])**2)/2
                   +3*(S[0,1]**2+S[0,2]**2+S[1,2]**2))
        samples.append([*p,volume,rad,np.degrees(angle)%360,tangent@S@tangent,vm,S[0,0]])
    samples=np.array(samples)
    sectors=[]
    for center in (45,135,225,315):
        select=(samples[:,4]>.22)&(samples[:,4]<.245)&(abs(samples[:,5]-center)<15)
        hoop=samples[select,6];weights=samples[select,3]
        sectors.append(dict(mesh=folder.name,region=f'rim_{center}deg_pm15',points=int(sum(select)),
                            hoop_mean_pa=float(np.average(hoop,weights=weights)),
                            hoop_min_pa=float(min(hoop)),hoop_max_pa=float(max(hoop))))
    interface=interface_forces(nodes,elements,materials,stress)
    m=7850*np.pi*.03**2*.09;omega=1800*2*np.pi/60;force=m*omega**2*.2
    for sign,r in interface.items():
        r['force']=r['force'].tolist()
        r['analytical_force_n']=force
        r['inward_force_n']=-sign*r['force'][0]
        r['force_error_percent']=100*(r['inward_force_n']/force-1)
        r['projected_nominal_stress_pa']=r['inward_force_n']/(.06*.02)
    disp=read_displacements(folder/'model.frd')
    recovered=consistent_insert_forces(nodes,elements,materials,
        {int(v[0]):np.array(v[1:4]) for v in blocks['displacements']})
    for r in recovered.values():
        r['analytical_force_n']=force;r['force_error_percent']=100*(r['inward_force_n']/force-1)
        r['projected_nominal_stress_pa']=r['inward_force_n']/(.06*.02)
    modal=[]
    polar_inertia=sum(v*rho*(p[0]**2+p[1]**2) for _,_,_,p,v,rho in q)
    for mode in range(1,9):
        u=disp[mode];energies=np.zeros(3);total=0.;insert=0.;angular_participation=0.
        coms={-1:np.zeros(3),1:np.zeros(3)};masses={-1:0.,1:0.}
        for eid,ip,N,p,volume,rho in q:
            value=N@np.array([u[n] for n in elements[eid]])
            mass=rho*volume;norm=value@value
            energies+=mass*value**2;total+=mass*norm
            angular_participation+=mass*value@np.array([-p[1],p[0],0.])
            if materials[eid]=='STEEL':
                insert+=mass*norm;sign=1 if p[0]>0 else -1
                coms[sign]+=mass*value;masses[sign]+=mass
        row=dict(mode=mode,frequency_hz=freqs[mode-1],x_fraction=energies[0]/total,
                 y_fraction=energies[1]/total,z_fraction=energies[2]/total,
                 insert_kinetic_fraction=insert/total,generalized_mass_kg=total)
        row['angular_acceleration_effective_inertia_fraction']=angular_participation**2/(total*polar_inertia)
        for sign,label in ((-1,'left'),(1,'right')):
            for axis,value in zip('xyz',coms[sign]/masses[sign]):row[f'{label}_com_{axis}']=value
        modal.append(row)
    peak=samples[np.argmax(samples[:,7])]
    details=dict(mesh=folder.name,rim=sectors,extrapolated_interface_tractions=interface,
                 consistent_interface_forces=recovered,
                 peak_vm=dict(x_m=peak[0],y_m=peak[1],z_m=peak[2],value_pa=peak[7]),
                 modes=modal,integrated_mass_kg=sum(v*rho for _,_,_,_,v,rho in q))
    (output/f'{folder.name}.json').write_text(json.dumps(details,indent=2))
    if plots:
        write_csv(output/'mode_summary.csv',modal)
        # Top faces from tetra connectivity, preserving nodes for displacement.
        boundary={}
        for eid,c in elements.items():
            for face in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
                key=tuple(sorted(c[i] for i in face));boundary.setdefault(key,[]).append(eid)
        faces=[key for key,v in boundary.items() if len(v)==1 and
               (all(abs(nodes[n][2]-.01)<1e-7 for n in key) or all(abs(nodes[n][2]-.045)<1e-7 for n in key))]
        fig,axes=plt.subplots(2,4,figsize=(16,8))
        for mode,ax in enumerate(axes.flat,1):
            maxu=max(np.linalg.norm(v) for v in disp[mode].values())
            polys=[];values=[]
            for face in faces:
                polys.append(np.array([nodes[n][:2] for n in face])*1000)
                values.append(np.mean([disp[mode][n][2] for n in face])/maxu)
            coll=PolyCollection(polys,array=np.array(values),cmap='coolwarm',clim=(-1,1),edgecolors='none')
            ax.add_collection(coll);ax.autoscale();ax.set_aspect('equal')
            # In-plane displacement arrows; normalized shape, not forced response.
            subset=list(nodes)[::max(1,len(nodes)//180)]
            p=np.array([nodes[n] for n in subset]);u=np.array([disp[mode][n] for n in subset])/maxu
            ax.quiver(p[:,0]*1000,p[:,1]*1000,u[:,0],u[:,1],scale=12,width=.004,color='black')
            ax.set_title(f'Mode {mode}: {freqs[mode-1]:.2f} Hz\nz-motion fraction {modal[mode-1]["z_fraction"]:.1%}',fontsize=10)
            ax.set_xticks([]);ax.set_yticks([])
        fig.suptitle('Prestressed carrier modes at 1800 rpm: color = normalized axial displacement; arrows = in-plane motion')
        fig.colorbar(coll,ax=axes.ravel().tolist(),shrink=.6,label='Axial displacement / maximum vector displacement')
        fig.savefig(output/'carrier_modes.png',dpi=160,bbox_inches='tight');plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(11,4.5))
        for ax,index,label in ((axes[0],6,'Hoop stress (MPa)'),(axes[1],7,'von Mises stress (MPa)')):
            sc=ax.scatter(samples[:,0]*1000,samples[:,1]*1000,c=samples[:,index]/1e6,s=3,cmap='viridis')
            ax.set_aspect('equal');ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)');fig.colorbar(sc,ax=ax,label=label)
        fig.suptitle('Aluminium integration-point stresses, projected through thickness')
        fig.tight_layout();fig.savefig(output/'stress_locations.png',dpi=160);plt.close(fig)
    return details


def main():
    output=ROOT/'results'/'stress_modes';output.mkdir(exist_ok=True)
    rows=[]
    for h in (.025,.018,.013):
        d=analyze(ROOT/'results'/'shaped_rotor'/f'h{h:g}_rpm1800',output,plots=h==.013)
        rows+=d['rim']
        print(json.dumps({k:v for k,v in d.items() if k!='modes'}),flush=True)
    write_csv(output/'rim_comparison.csv',rows)


if __name__=='__main__':main()
